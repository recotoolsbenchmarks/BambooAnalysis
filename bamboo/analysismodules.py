"""
Minimally, ``bambooRun`` needs a class with a constructor that takes a single argument
(the list of command-line arguments that it does not recognize as its own), and a
``run`` method  that takes no arguments.
:py:mod:`bamboo.analysismodules` provides more interesting base classes, starting from
:py:class:`~bamboo.analysismodules.AnalysisModule`, which implements a large part of
the common functionality for loading samples and distributing worker tasks.
:py:class:`~bamboo.analysismodules.HistogramsModule` specializes this further
for modules that output stack histograms, and
:py:class:`~bamboo.analysismodules.NanoAODHistoModule` supplements this
with loading the decorations for NanoAOD, and merging of the counters for generator weights etc.
"""
import argparse
from collections import defaultdict
from itertools import chain
from functools import partial
import logging
logger = logging.getLogger(__name__)
import os.path
import re
import datetime
from timeit import default_timer as timer
import resource
import urllib.parse
from .analysisutils import addLumiMask, downloadCertifiedLumiFiles, parseAnalysisConfig, getAFileFromAnySample, readEnvConfig

def reproduceArgv(args, group):
    # Reconstruct the module-specific arguments (to pass them to the worker processes later on)
    assert isinstance(group, argparse._ArgumentGroup)
    argv = []
    for action in group._group_actions:
        if isinstance(action, argparse._StoreTrueAction):
            if getattr(args, action.dest):
                argv.append(action.option_strings[0])
        elif isinstance(action, argparse._StoreAction):
            argv.append(action.option_strings[0])
            argv.append(getattr(args, action.dest))
        elif isinstance(action, argparse._AppendAction):
            items = getattr(args, action.dest)
            if items:
                for item in items:
                    argv.append(action.option_strings[0])
                    argv.append(item)
        else:
            raise RuntimeError("Reconstruction of action {0} not supported".format(action))
    return argv

def modAbsPath(modArg):
    # Put absolute path if module is specified by file
    mod_clName = None
    if ":" in modArg:
        modArg, mod_clName = modArg.split(":")
    if os.path.isfile(modArg):
        modArg = os.path.abspath(modArg)
    if mod_clName:
        modArg = ":".join((modArg, mod_clName))
    return modArg

def parseRunRange(rrStr):
    return tuple(int(t.strip()) for t in rrStr.split(","))

def parseEras(eraStr):
    if ":" in eraStr:
        eraMode, eras = eraStr.split(":")
        return (eraMode, list(eras.split(",")))
    else:
        if eraStr in ("all", "combined", "split"):
            return (eraStr, None)
        else:
            return ("all", list(eraStr.split(",")))

class SampleTask:
    """ Information about work to be done for single sample (input and output file, parameters etc.) """
    def __init__(self, name, inputFiles=None, outputFile=None, kwargs=None, config=None, resolver=None):
        self.name = name
        self._inputFiles = inputFiles
        self.outputFile = outputFile
        self.kwargs = kwargs
        self.config = config
        self._resolver = resolver
    @property
    def inputFiles(self):
        if self._inputFiles is None:
            self._inputFiles = self._resolver(self.name, self.config)
        return self._inputFiles

def getBatchDefaults(backend):
    if backend == "slurm":
        return {
            "sbatch_time"     : "0-01:00",
            "sbatch_mem"      : "2048",
            "stageoutFiles"   : ["*.root"],
            "sbatch_workdir"  : os.getcwd(),
            "sbatch_additionalOptions" : [ "--export=ALL" ],
            }
    elif backend == "htcondor":
        return {"cmd" : [
            "universe     = vanilla",
            "+MaxRuntime  = {0:d}".format(20*60), # 20 minutes
            "getenv       = True"
            ]}

class AnalysisModule:
    """ Base analysis module
    
    Adds common infrastructure for parsing analysis config files
    and running on a batch system, with customization points for
    concrete classes to implement (most importantly
    :py:meth:`~bamboo.analysismodules.AnalysisModule.processTrees`
    and :py:meth:`~bamboo.analysismodules.AnalysisModule.postProcess`)
    """
    CustomSampleAttributes = ["db", "split", "files", "run_range", "certified_lumi_file"]
    def __init__(self, args):
        """ Constructor

        set up argument parsing, calling :py:meth:`~bamboo.analysismodules.AnalysisModule.addArgs`
        and :py:meth:`~bamboo.analysismodules.AnalysisModule.initialize`

        :param args: list of command-line arguments that are not parsed by ``bambooRun``
        """
        parser = argparse.ArgumentParser(description=(
            "Run an analysis, i.e. process the samples in an analysis description file with a module (subclass of bamboo.analysismodules.AnalysisModule). "
            "There are three modes, specified by the --distributed option: if unspecified, one program processes all samples and collects the results; "
            "--distributed=driver does the same, but launches worker tasks "
            "(the same program with --distributed=worker, therefore some of the options only apply to 'driver' or 'worker' mode, or have a different interpretation) "
            "to process the samples, for instance on a batch system (depending on the settings in the --envConfig file)."))
        parser.add_argument("-m", "--module", type=str, default="bamboo.analysismodules:AnalysisModule", help="Module to run (format: modulenameOrPath[:classname])")
        parser.add_argument("-v", "--verbose", action="store_true", help="Run in verbose mode")
        parser.add_argument("--distributed", type=str, help="Role in distributed mode (sequential mode if not specified)", choices=["worker", "driver", "finalize"])
        parser.add_argument("input", nargs="*", help="Input: analysis description yml file (driver mode) or files to process (worker mode)")
        parser.add_argument("-o", "--output", type=str, default=".", help="Output directory (driver mode) or file (worker mode) name")
        parser.add_argument("--interactive", "-i", action="store_true", help="Interactive mode (initialize to an IPython shell for exploration)")
        parser.add_argument("--maxFiles", type=int, default=-1, help="Maximum number of files to process per sample (all by default, 1 may be useful for tests)")
        parser.add_argument("-t", "--threads", type=int, help="Enable implicit multithreading, specify number of threads to launch")
        driver = parser.add_argument_group("non-worker mode only (--distributed=driver, finalize, or unspecified) optional arguments")
        driver.add_argument("--redodbqueries", action="store_true", help="Redo all DAS/SAMADhi queries even if results can be read from cache files")
        driver.add_argument("--overwritesamplefilelists", action="store_true", help="Write DAS/SAMADhi results to files even if files exist (meaningless without --redodbqueries)")
        driver.add_argument("--envConfig", type=str, help="Config file to read computing environment configuration from (batch system, storage site etc.)")
        driver.add_argument("--onlypost", action="store_true", help="Only run postprocessing step on previous results")
        driver.add_argument("--eras", type=parseEras, default=("all", None), help="Select eras to consider, and which plots to make (format: '[all|split|combined]:[era1,era2...]')")
        worker = parser.add_argument_group("worker mode only (--distributed=worker) arguments")
        worker.add_argument("--treeName", type=str, default="Events", help="Tree name (default: Events)")
        worker.add_argument("--runRange", type=parseRunRange, help="Run range (format: 'firstRun,lastRun')")
        worker.add_argument("--certifiedLumiFile", type=str, help="(local) path of a certified lumi JSON file")
        worker.add_argument("--sample", type=str, help="Sample name (as in the samples section of the analysis configuration)")
        worker.add_argument("--input", action="append", help="File with the list of files to read", dest="filelists")
        worker.add_argument("--anaConfig", type=str, default=None, help="Analysis description yml file provided to the driver")
        specific = parser.add_argument_group("module-specific arguments")
        self.addArgs(specific)
        self.args = parser.parse_args(args)
        self.specificArgv = reproduceArgv(self.args, specific)
        self._envConfig = None
        self._analysisConfig = None
        self.initialize()
    def addArgs(self, parser):
        """ Hook for adding module-specific argument parsing (receives an argument group), parsed arguments are available in ``self.args`` afterwards """
        pass
    def initialize(self):
        """ Hook for module-specific initialization (called from the constructor after parsing arguments) """
        pass
    @property
    def inputs(self):
        inputs = list(self.args.input)
        if self.args.distributed == "worker" and self.args.filelists:
            for ifl in self.args.filelists:
                with open(ifl) as iflf:
                    inputs += [ ln.strip() for ln in iflf if len(ln.strip()) > 0 ]
        return inputs
    @property
    def envConfig(self):
        if self._envConfig is None:
            self._envConfig = readEnvConfig(self.args.envConfig)
        return self._envConfig
    @property
    def analysisConfig(self):
        if self._analysisConfig is None:
            if self.args.distributed == "worker":
                anaCfgName = self.args.anaConfig
            elif ( not self.args.distributed ) or self.args.distributed in ("driver", "finalize"):
                anaCfgName = self.args.input[0]
            else:
                raise RuntimeError("--distributed should be either worker, driver, finalize, or unspecified (for sequential mode)")
            self._analysisConfig = parseAnalysisConfig(anaCfgName)
        return self._analysisConfig
    def getSampleFilesResolver(self, cfgDir="."):
        from .analysisutils import sample_resolveFiles
        resolver = partial(sample_resolveFiles, redodbqueries=self.args.redodbqueries, overwritesamplefilelists=self.args.overwritesamplefilelists, envConfig=self.envConfig, cfgDir=cfgDir)
        if self.args.maxFiles <= 0:
            return resolver
        else:
            def maxFilesResolver(smpName, smpConfig, maxFiles=-1, resolve=None):
                inputFiles = resolve(smpName, smpConfig)
                if maxFiles > 0 and maxFiles < len(inputFiles):
                    logger.warning("Only processing first {0:d} of {1:d} files for sample {2}".format(maxFiles, len(inputFiles), smpName))
                    inputFiles = inputFiles[:maxFiles]
                return inputFiles
            return partial(maxFilesResolver, maxFiles=self.args.maxFiles, resolve=resolver)
    def getATree(self, fileName=None, sampleName=None):
        """ Retrieve a representative TTree, e.g. for defining the plots or interactive inspection, and a dictionary with metadata """
        if self.args.distributed == "worker":
            from .root import gbl
            tup = gbl.TChain(self.args.treeName)
            if not tup.Add(self.inputs[0], 0):
                raise IOError("Could not open file {}".format(self.inputs[0]))
            return tup, "", {}
        elif ( not self.args.distributed ) or self.args.distributed in ("driver", "finalize"):
            if len(self.args.input) != 1:
                raise RuntimeError("Main process (driver, finalize, or non-distributed) needs exactly one argument (analysis description YAML file)")
            analysisCfg = self.analysisConfig
            if fileName and sampleName:
                sampleCfg = analysisCfg["samples"][sampleName]
            else:
                filesResolver = self.getSampleFilesResolver(cfgDir=os.path.dirname(os.path.abspath(anaCfgName)))
                sampleName,sampleCfg,fileName = getAFileFromAnySample(analysisCfg["samples"], resolveFiles=filesResolver)
            logger.debug("getATree: using a file from sample {0} ({1})".format(sampleName, fileName))
            from .root import gbl
            tup = gbl.TChain(analysisCfg.get("tree", "Events"))
            if not tup.Add(fileName, 0):
                raise IOError("Could not open file {}".format(fileName))
            return tup, sampleName, sampleCfg
        else:
            raise RuntimeError("--distributed should be either worker, driver, finalize, or unspecified (for sequential mode)")
    def customizeAnalysisCfg(self, analysisCfg):
        """ Hook to modify the analysis configuration before jobs are created (only called in driver or non-distributed mode) """
        pass
    def run(self):
        """ Main method

        Depending on the arguments passed, this will:

        * if ``-i`` or ``--interactive``, call :py:meth:`~bamboo.analysismodules.AnalysisModule.interact`
          (which could do some initialization and start an IPython shell)
        * if ``--distributed=worker`` call :py:meth:`~bamboo.analysismodules.AnalysisModule.processTrees`
          with the appropriate input, output, treename, lumi mask and run range
        * if ``--distributed=driver`` or not given (sequential mode): parse the analysis configuration file,
          construct the tasks with :py:meth:`~bamboo.analysismodules.AnalysisModule.getTasks`, run them
          (on a batch cluster or in the same process with :py:meth:`~bamboo.analysismodules.AnalysisModule.processTrees`),
          and finally call :py:meth:`~bamboo.analysismodules.AnalysisModule.postProcess` with the results.
        """
        if self.args.interactive:
            self.interact()
        else:
            if self.args.distributed == "worker":
                if ( not self.args.output.endswith(".root") ) or os.path.isdir(self.args.output):
                    raise RuntimeError("Output for worker processes needs to be a ROOT file")
                inputFiles = self.inputs
                if self.args.maxFiles > 0 and self.args.maxFiles < len(inputFiles):
                    logger.warning("Only processing first {0:d} of {1:d} files".format(self.args.maxFiles, len(inputFiles)))
                    inputFiles = inputFiles[:self.args.maxFiles]
                logger.info("Worker process: calling processTrees for {mod} with ({0}, {1}, treeName={treeName}, certifiedLumiFile={certifiedLumiFile}, runRange={runRange}, sample={sample})".format(inputFiles, self.args.output, mod=self.args.module, treeName=self.args.treeName, certifiedLumiFile=self.args.certifiedLumiFile, runRange=self.args.runRange, sample=self.args.sample))
                if self.args.anaConfig:
                    sampleCfg = self.analysisConfig["samples"][self.args.sample]
                else:
                    sampleCfg = None
                if self.args.threads:
                    from .root import gbl
                    logger.info(f"Enabling implicit MT for {self.args.threads} threads")
                    gbl.ROOT.EnableImplicitMT(self.args.threads)
                self.processTrees(inputFiles, self.args.output, tree=self.args.treeName, certifiedLumiFile=self.args.certifiedLumiFile, runRange=self.args.runRange, sample=self.args.sample, sampleCfg=sampleCfg)
            elif ( not self.args.distributed ) or self.args.distributed in ("driver", "finalize"):
                if len(self.args.input) != 1:
                    raise RuntimeError("Main process (driver or non-distributed) needs exactly one argument (analysis description YAML file)")
                anaCfgName = self.args.input[0]
                workdir = self.args.output
                analysisCfg = self.analysisConfig
                self.customizeAnalysisCfg(analysisCfg)
                filesResolver = self.getSampleFilesResolver(cfgDir=os.path.dirname(os.path.abspath(anaCfgName)))
                tasks = self.getTasks(analysisCfg, tree=analysisCfg.get("tree", "Events"), resolveFiles=(filesResolver if not self.args.onlypost else None))
                resultsdir = os.path.join(workdir, "results")
                if self.args.onlypost:
                    if os.path.exists(resultsdir):
                        aProblem = False
                        for tsk in tasks:
                            fullOutName = os.path.join(resultsdir, tsk.outputFile)
                            if not os.path.exists(fullOutName):
                                logger.error(f"Output file for {tsk.name} not found ({fullOutName})")
                                aProblem = True
                        if aProblem:
                            logger.error(f"Not all output files were found, cannot perform post-processing")
                            return
                    else:
                        logger.error(f"Results directory {resultsdir} does not exist, cannot perform post-processing")
                        return
                elif self.args.distributed == "finalize":
                    tasks_notfinalized = [ tsk for tsk in tasks if not os.path.exists(os.path.join(resultsdir, tsk.outputFile)) ]
                    if not tasks_notfinalized:
                        logger.info(f"All output files were found, so no finalization was redone. If you merged the outputs of partially-done MC samples please remove them from {resultsdir} and rerun to pick up the rest.")
                    else:
                        def cmdMatch(ln, smpNm):
                            return f" --sample={smpNm} " in ln or ln.endswith(f" --sample={smpNm}")
                        from .batch import getBackend
                        batchBackend = getBackend(self.envConfig["batch"]["backend"])
                        outputs = batchBackend.findOutputsForCommands(os.path.join(workdir, "batch"),
                                { tsk.name: partial(cmdMatch, smpNm=tsk.name) for tsk in tasks_notfinalized })
                        aProblem = False
                        for tsk in tasks_notfinalized:
                            nExpected, tskOut = outputs[tsk.name]
                            if not tskOut:
                                logger.error(f"No output files for sample {tsk.name}")
                                aProblem = True
                            tskOut_by_name = defaultdict(list)
                            for fn in tskOut:
                                tskOut_by_name[os.path.basename(fn)].append(fn)
                            for outFileName, outFiles in tskOut_by_name.items():
                                if nExpected != len(outFiles):
                                    logger.error(f"Not all jobs for {tsk.name} produced an output file {outFileName} ({len(outFiles):d}/{nExpected:d} found), cannot finalize")
                                    aProblem = True
                                else:
                                    haddCmd = ["hadd", "-f", os.path.join(resultsdir, outFileName)]+outFiles
                                    import subprocess
                                    try:
                                        logger.debug("Merging outputs for sample {0} with {1}".format(tsk.name, " ".join(haddCmd)))
                                        subprocess.check_call(haddCmd, stdout=subprocess.DEVNULL)
                                    except subprocess.CalledProcessError:
                                        loger.error("Failed to run {0}".format(" ".join(haddCmd)))
                                        aProblem = True
                        if aProblem:
                            logger.error("Could not finalize all tasks so no post-processing will be run (rerun in verbose mode for the full list of directories and commands)")
                            return
                        else:
                            logger.info("All tasks finalized")
                elif not tasks:
                    logger.warning("No tasks defined, skipping to postprocess")
                else: ## need to run tasks
                    downloadCertifiedLumiFiles(tasks, workdir=workdir)
                    if os.path.exists(resultsdir):
                        logger.warning("Output directory {0} exists, previous results may be overwritten".format(resultsdir))
                    else:
                        os.makedirs(resultsdir)
                    ## store one "skeleton" tree (for more efficient "onlypost" later on
                    aTask = tasks[0]
                    aFileName = aTask.inputFiles[0]
                    from .root import gbl
                    aFile = gbl.TFile.Open(aFileName)
                    if not aFile:
                        logger.warning(f"Could not open file {aFileName}, no skeleton tree will be saved")
                    else:
                        treeName = analysisCfg.get("tree", "Events")
                        aTree = aFile.Get(treeName)
                        if not aTree:
                            logger.warning(f"Could not get {treeName} from file {aFileName}, no skeleton tree will be saved")
                        else:
                            outfName = os.path.join(resultsdir, "__skeleton__{0}.root".format(aTask.kwargs["sample"]))
                            outf = gbl.TFile.Open(outfName, "RECREATE")
                            if "/" in treeName:
                                outf.mkdir("/".join(treeName.split("/")[:-1])).cd()
                            skeletonTree = aTree.CloneTree(1) ## copy header and a single event
                            outf.Write()
                            outf.Close()
                            logger.debug(f"Skeleton tree written to {outfName}")
                    ## run all tasks
                    if not self.args.distributed: ## sequential mode
                        for task in tasks:
                            output = os.path.join(resultsdir, task.outputFile)
                            logger.info("Sequential mode: calling processTrees for {mod} with ({0}, {1}, {2}".format(task.inputFiles, output, ", ".join("{0}={1}".format(k,v) for k,v in task.kwargs.items()), mod=self.args.module))
                            if "runRange" in task.kwargs:
                                task.kwargs["runRange"] = parseRunRange(task.kwargs["runRange"])
                            if self.args.threads:
                                from .root import gbl
                                logger.info(f"Enabling implicit MT for {self.args.threads} threads")
                                gbl.ROOT.EnableImplicitMT(self.args.threads)
                            self.processTrees(task.inputFiles, output, sampleCfg=task.config, **task.kwargs)
                    elif self.args.distributed == "driver":
                        ## construct the list of tasks
                        from .batch import splitInChunks, writeFileList, SplitAggregationTask, HaddAction, format_runtime, getBackend
                        commArgs = [
                              "bambooRun"
                            , "--module={0}".format(modAbsPath(self.args.module))
                            , "--distributed=worker"
                            , "--anaConfig={0}".format(os.path.abspath(anaCfgName))
                        ] + self.specificArgv + (["--verbose"] if self.args.verbose else []) + ([f"-t {self.args.threads}"] if self.args.threads else [])
                        beTasks = []
                        for tsk in tasks:
                            split = 1
                            if tsk.config and "split" in tsk.config:
                                split = tsk.config["split"]
                            if split >= 0: ## at least 1 (no splitting), at most the numer of arguments (one job per input)
                                chunks = splitInChunks(tsk.inputFiles, nChunks=max(1, min(split, len(tsk.inputFiles))))
                            else: ## at least 1 (one job per input), at most the number of arguments (no splitting)
                                chunks = splitInChunks(tsk.inputFiles, chunkLength=max(1, min(-split, len(tsk.inputFiles))))
                            cmds = []
                            os.makedirs(os.path.join(workdir, "infiles"), exist_ok=True)
                            for i,chunk in enumerate(chunks):
                                cfn = os.path.join(workdir, "infiles", "{0}_in_{1:d}.txt".format(tsk.kwargs["sample"], i))
                                writeFileList(chunk, cfn)
                                cmds.append(" ".join([str(a) for a in commArgs] + [
                                    "--input={0}".format(os.path.abspath(cfn)), "--output={0}".format(tsk.outputFile)]+
                                    ["--{0}={1}".format(key, value) for key, value in tsk.kwargs.items()]
                                    ))
                            beTasks.append(SplitAggregationTask(cmds, finalizeAction=HaddAction(cmds, outDir=resultsdir, options=["-f"])))
                        ## submit to backend
                        backend = self.envConfig["batch"]["backend"]
                        batchBackend = getBackend(backend)
                        clusJobs = batchBackend.jobsFromTasks(beTasks, workdir=os.path.join(workdir, "batch"), batchConfig=self.envConfig.get(backend), configOpts=getBatchDefaults(backend))
                        for j in clusJobs:
                            j.submit()
                        logger.info("The status of the batch jobs will be periodically checked, and the outputs merged if necessary. "
                                "If only few jobs (or the monitoring loop) fail, it may be more efficient to resubmit or rerun them manually "
                                "and (if necessary) merge the outputs and produce the final results either by rerunning with --driver=finalize, "
                                "or manually, using the commands that will be printed if any jobs fail, and by rerunning with the --onlypost option.")
                        clusMon = batchBackend.makeTasksMonitor(clusJobs, beTasks, interval=int(self.envConfig["batch"].get("update", 120)))
                        collectResult = clusMon.collect() ## wait for batch jobs to finish and finalize

                        if any(not tsk.failedCommands for tsk in beTasks):
                            ## Print time report (possibly more later)
                            logger.info("Average runtime for successful tasks (to further tune job splitting):")
                            for tsk in beTasks:
                                if not tsk.failedCommands:
                                    totTime = sum((next(clus for clus in clusJobs if cmd in clus.commandList).getRuntime(cmd) for cmd in tsk.commandList), datetime.timedelta())
                                    nTasks = len(tsk.commandList)
                                    smpName = next(arg for arg in tsk.commandList[0].split() if arg.startswith("--sample=")).split("=")[1]
                                    logger.info(" - {0}: {1} ({2:d} jobs, total: {3})".format(smpName, format_runtime(totTime/nTasks), nTasks, format_runtime(totTime)))

                        if not collectResult["success"]:
                            # Print missing hadd actions to be done when (if) those recovery jobs succeed
                            haddCmds = list(chain.from_iterable(tsk.finalizeAction.getActions() for tsk in beTasks if tsk.failedCommands and tsk.finalizeAction))
                            logger.error("Finalization commands to be run are:\n{0}".format("\n".join(" ".join(cmd) for cmd in haddCmds)))
                            logger.error("Since there were failed jobs, I'm not doing the postprocessing step. "
                                    "After performing recovery actions (see above) you may run me again with the --distributed=finalize (to merge) or --onlypost (if merged manually) option instead.")
                            return
                try:
                    self.postProcess(tasks, config=analysisCfg, workdir=workdir, resultsdir=resultsdir)
                except Exception as ex:
                    logger.exception(ex)
                    logger.error("Exception in postprocessing. If the worker job results (e.g. histograms) were correctly produced, you do not need to resubmit them, and may recover by running with the --onlypost option instead.")
            else:
                raise RuntimeError("--distributed should be either worker, driver, or be unspecified (for sequential mode)")

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, sample=None, sampleCfg=None):
        """ worker method: produce results (e.g. histograms or trees) from the input files

        should be implemented by concrete modules

        :param inputFiles: input file names
        :param outputFile: output file name
        :param tree: key name of the tree inside the files
        :param certifiedLumiFile: lumi mask json file name
        :param runRange: run range to consider (for efficiency of the lumi mask)
        :param sample: sample name (key in the samples block of the configuration file)
        :param sampleCfg: that sample's entry in the configuration file
        """
        pass
    def getTasks(self, analysisCfg, resolveFiles=None, **extraOpts):
        """ Get tasks from analysis configs (and args), called in for driver or sequential mode

        :returns: a list of :py:class:`~bamboo.analysismodules.SampleTask` instances
        """
        tasks = []
        sel_eras = analysisCfg["eras"].keys()
        if self.args.eras[1]:
            sel_eras = list(era for era in sel_eras if era in self.args.eras[1])
        for sName, sConfig in analysisCfg["samples"].items():
            opts = dict(extraOpts)
            if "certified_lumi_file" in sConfig:
                opts["certifiedLumiFile"] = sConfig.get("certified_lumi_file")
            if "run_range" in sConfig:
                opts["runRange"] = ",".join(str(rn) for rn in sConfig.get("run_range"))
            opts["sample"] = sName
            if "era" not in sConfig or sConfig["era"] in sel_eras:
                tasks.append(SampleTask(sName, outputFile=f"{sName}.root", kwargs=opts, config=sConfig, resolver=resolveFiles))
        return tasks

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Do postprocessing on the results of the tasks, if needed

        should be implemented by concrete modules

        :param taskList: ``(inputs, output), kwargs`` for the tasks (list, string, and dictionary)
        :param config: parsed analysis configuration file
        :param workdir: working directory for the current run
        :param resultsdir: path with the results files
        """
        pass
    def interact(self):
        """ Interactive mode (load some things and embed IPython)

        should be implemented by concrete modules
        """
        pass ## define things and embed IPython

class HistogramsModule(AnalysisModule):
    """ Base histogram analysis module """
    def __init__(self, args):
        """ Constructor

        Defines a ``plotList`` member variable, which will store a list of plots,
        (:py:meth:`~bamboo.analysismodules.HistogramsModule.processTrees` will call
        :py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots`, store the result there,
        and retrieve their results afterwards to trigger the actual processing of the files)
        """
        super(HistogramsModule, self).__init__(args)
        self.plotList = []
        self.plotDefaults = {}

    def addArgs(self, parser):
        super(HistogramsModule, self).addArgs(parser)
        parser.add_argument("--plotIt", type=str, default="plotIt", help="plotIt executable to use (default is taken from $PATH)")

    def initialize(self):
        """ initialize """
        if self.args.distributed == "worker" and len(self.inputs) == 0:
            raise RuntimeError("Worker task needs at least one input file")

    def interact(self):
        """ Interactively inspect a decorated input tree

        Available variables: ``tree`` (decorated tree), ``tup`` (raw tree),
        ``noSel`` (root selection), ``backend``, ``runAndLS`` (e.g. ``(runExpr, lumiBlockExpr)``)
        (the inputs for the lumi mask), and ``op`` (:py:mod:`bamboo.treefunctions`).
        """
        tup, smpName, smpCfg = self.getATree()
        tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=smpName, sampleCfg=smpCfg)
        import bamboo.treefunctions as op
        import IPython
        IPython.embed()

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, sample=None, sampleCfg=None):
        """ Worker sequence: produce histograms from the input files

        More in detail, this will load the inputs, call :py:meth:`~bamboo.analysismodules.HistogramsModule.prepareTree`,
        add a lumi mask if requested, call :py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots`,
        run over all files, and write the produced histograms to the output file.
        """
        from .root import gbl
        tup = gbl.TChain(tree)
        for fName in inputFiles:
            if not tup.Add(fName, 0):
                raise IOError("Could not open file {}".format(fName))
        tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=sample, sampleCfg=sampleCfg)
        if certifiedLumiFile:
            noSel = addLumiMask(noSel, certifiedLumiFile, runRange=runRange, runAndLS=runAndLS)

        outF = gbl.TFile.Open(outputFile, "RECREATE")
        logger.info("Starting to define plots")
        start = timer()
        self.plotList = self.definePlots(tree, noSel, sample=sample, sampleCfg=sampleCfg)
        end = timer()
        maxrssmb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        logger.info(f"{len(self.plotList):d} plots defined in {end - start:.2f}s, max RSS: {maxrssmb:.2f}MB")
        if hasattr(backend, "buildGraph"):
            logger.info("Starting to build RDataFrame graph")
            start = timer()
            backend.buildGraph(self.plotList)
            end = timer()
            maxrssmb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
            logger.info(f"RDataFrame graph finished in {end - start:.2f}s, max RSS: {maxrssmb:.2f}MB")
        from .dataframebackend import _RDFNodeStats, _RDFHistoNDStats
        logger.info(f"Number of uses per node type: {_RDFNodeStats!s}")
        logger.info(f"HistoND calls per column type: {_RDFHistoNDStats!s}")
        from .dataframebackend import _RDFHistoND_methods
        logger.debug(f"HistoND helper instantiations: {_RDFHistoND_methods!s}")
        ## make a list of suggested nuisance parameters
        systNuis = []
        for systN, systVars in backend.allSysts.items():
            for varn in systVars:
                for suff in ("up", "down"):
                    if varn.endswith(suff):
                        varn = varn[:-len(suff)]
                if varn not in systNuis:
                    systNuis.append(varn)
        if systNuis:
            logger.info("Systematic shape variations impacting any plots: {0}".format(", ".join(systNuis)))

        outF.cd()
        logger.info("Starting to fill plots")
        start = timer()
        numHistos = 0
        for p in self.plotList:
            for h in backend.getResults(p):
                numHistos += 1
                h.Write()
        end = timer()
        maxrssmb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        logger.info(f"Plots finished in {end - start:.2f}s, max RSS: {maxrssmb:.2f}MB ({numHistos:d} histograms)")
        self.mergeCounters(outF, inputFiles, sample=sample)
        outF.Close()
        return backend
    # processTrees customisation points
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        """ Create decorated tree, selection root (noSel), backend, and (run,LS) expressions

        should be implemented by concrete modules

        :param tree: decorated tree
        :param sample: sample name (as in the samples section of the analysis configuration file)
        :param sampleCfg: that sample's entry in the configuration file
        """
        return tree, None, None, None
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        """ Main method: define plots on the trees (for a give systematic variation)

        should be implemented by concrete modules, and return a list of
        :py:class:`bamboo.plots.Plot` objects.
        The structure (name, binning) of the histograms should not depend on the sample, era,
        and the list should be the same for all values
        (the weights and systematic variations associated with weights or collections
        may differ for data and different MC samples, so the actual set of histograms
        will not be identical).

        :param tree: decorated tree
        :param noSel: base selection
        :param sample: sample name (as in the samples section of the analysis configuration file)
        :param sampleCfg: that sample's entry in the configuration file
        """
        return [] ## plot list
    def mergeCounters(self, outF, infileNames, sample=None):
        """ Merge counters

        should be implemented by concrete modules

        :param outF: output file (TFile pointer)
        :param infileNames: input file names
        :param sample: sample name
        """
        pass
    def readCounters(self, resultsFile):
        """ Read counters from results file

        should be implemented by concrete modules, and return a dictionary with
        counter names and the corresponding sums

        :param resultsFile: TFile pointer to the results file
        """
        return dict()

    def getPlotList(self, fileHint=None, sampleHint=None, resultsdir=None):
        """
        Helper method for postprocessing: construct the plot list

        The path (and sample name) of an input file can be specified,
        otherwise the results directory is searched for a skeleton tree.
        Please note that in the latter case, the skeleton file is arbitrary
        (in practice it probably corresponds to the first sample encountered
        when running in sequential or ``--distributed=driver`` mode), so if
        the postprocessing depends on things that are different between
        samples, one needs to be extra careful to avoid surprises.

        :param fileHint: name of an input file for one of the samples
        :param sampleHint: sample name for the input file passed in ``fileHint``
        :param resultsdir: directory with the produced results files (mandatory if no ``fileHint`` and ``sampleHint`` are passed)
        """
        if fileHint is not None and sampleHint is not None:
            pass
        elif resultsdir is not None:
            try:
                import os
                prefix = "__skeleton__"
                suffix = ".root"
                skelFn = next(fn for fn in os.listdir(resultsdir) if fn.startswith(prefix) and fn.endswith(suffix))
                fileHint = os.path.join(resultsdir, skelFn)
                sampleHint = skelFn[len(prefix):-len(suffix)]
            except StopIteration:
                raise RuntimeError(f"No skeleton file found in {resultsdir}")
        else:
            raise RuntimeError("Either the results directory, or an input file and corresponding sample name, needs to be specified")
        tup, smpName, smpCfg = self.getATree(fileName=fileHint, sampleName=sampleHint)
        tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=smpName, sampleCfg=smpCfg)
        if "certified_lumi_file" in smpCfg:
            lumiFile = os.path.join(workdir, urllib.parse.urlparse(smpCfg["certified_lumi_file"]).path.split("/")[-1])
            noSel = addLumiMask(noSel, lumiFile, runRange=smpCfg.get("run_range"), runAndLS=runAndLS)
        return self.definePlots(tree, noSel, sample=smpName, sampleCfg=smpCfg)

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: run plotIt

        The list of plots is created if needed (from a representative file,
        this enables rerunning the postprocessing step on the results files),
        and then plotIt is executed
        """
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [ ap for ap in self.plotList if isinstance(ap, CutFlowReport) ]
        plotList_plotIt = [ ap for ap in self.plotList if ( isinstance(ap, Plot) or isinstance(ap, DerivedPlot) ) and len(ap.binnings) == 1 ]
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(config, plotList_cutflowreport, resultsdir=resultsdir, readCounters=self.readCounters, eras=self.args.eras, verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            eraMode, eras = self.args.eras
            if eras is None:
                eras = list(config["eras"].keys())
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=(eraMode, eras), workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt, eras=(eraMode, eras), verbose=self.args.verbose)

class NanoAODModule(AnalysisModule):
    """ A :py:class:`~bamboo.analysismodules.AnalysisModule` extension for NanoAOD, adding decorations and merging of the counters """
    def isMC(self, sampleName):
        return not any(sampleName.startswith(pd) for pd in ("BTagCSV", "BTagMu", "Charmonium", "DisplacedJet", "DoubleEG", "DoubleMuon", "DoubleMuonLowMass", "EGamma", "FSQJet1", "FSQJet2", "FSQJets", "HTMHT", "HeavyFlavour", "HighEGJet", "HighMultiplicity", "HighPtLowerPhotons", "IsolatedBunch", "JetHT", "MET", "MinimumBias", "MuOnia", "MuonEG", "NoBPTX", "SingleElectron", "SingleMuon", "SinglePhoton", "Tau", "ZeroBias"))
    def prepareTree(self, tree, sample=None, sampleCfg=None, description=None, lazyBackend=False):
        """ Add NanoAOD decorations, and create an RDataFrame backend

        In addition to the arguments needed for the base class
        :py:meth:`~bamboo.analysismodules.AnalysisModule.prepareTree`` method,
        a description of the tree, and settings for reading systematic variations
        or corrections from alternative branches, or calculating these on the fly,
        should be passed, such that the decorations can be constructed accordingly.

        :param description: description of the tree format, and configuration for reading or calculating systematic variations and corrections, a :py:class:`~bamboo.treedecorators.NanoAODDescription` instance (see also :py:meth:`bamboo.treedecorators.NanoAODDescription.get`)
        """
        from bamboo.treedecorators import decorateNanoAOD
        from bamboo.dataframebackend import DataframeBackend, LazyDataframeBackend
        backendCls = (LazyDataframeBackend if lazyBackend else DataframeBackend)
        t = decorateNanoAOD(tree, description=description)
        be, noSel = backendCls.create(t)
        return t, noSel, be, (t.run, t.luminosityBlock)
    def mergeCounters(self, outF, infileNames, sample=None):
        """ Merge the ``Runs`` trees """
        from .root import gbl
        cruns = gbl.TChain("Runs")
        for fn in infileNames:
            cruns.Add(fn)
        outF.cd()
        runs = cruns.CloneTree()
        runs.Write("Runs")
    def readCounters(self, resultsFile):
        """ Sum over each leaf of the (merged) ``Runs`` tree (except ``run``) """
        runs = resultsFile.Get("Runs")
        from .root import gbl
        if ( not runs ) or ( not isinstance(runs, gbl.TTree) ):
            raise RuntimeError("No tree with name 'Runs' found in {0}".format(resultsFile.GetName()))
        sums = dict()
        runs.GetEntry(0)
        for lv in runs.GetListOfLeaves():
            lvn = lv.GetName()
            if lvn != "run":
                if lv.GetLeafCount():
                    lvcn = lv.GetLeafCount().GetName()
                    if lvcn in sums:
                        del sums[lvcn]
                    sums[lvn] = [ lv.GetValue(i) for i in range(lv.GetLeafCount().GetValueLong64()) ]
                else:
                    sums[lvn] = lv.GetValue()
        for entry in range(1, runs.GetEntries()):
            runs.GetEntry(entry)
            for cn, vals in sums.items():
                if hasattr(vals, "__iter__"):
                    entryvals = getattr(runs, cn)
                    ## warning and workaround (these should be consistent for all NanoAODs in a sample)
                    if len(vals) != len(entryvals):
                        logger.error("Runs tree: array of sums {0} has a different length in entry {1:d}: {2:d} (expected {3:d})".format(cn, entry, len(entryvals), len(vals)))
                    for i in range(min(len(vals), len(entryvals))):
                        vals[i] += entryvals[i]
                else:
                    sums[cn] += getattr(runs, cn)
        return sums

class NanoAODHistoModule(NanoAODModule, HistogramsModule):
    """ A :py:class:`~bamboo.analysismodules.HistogramsModule` implementation for NanoAOD, adding decorations and merging of the counters """
    def __init__(self, args):
        super(NanoAODHistoModule, self).__init__(args)

class SkimmerModule(AnalysisModule):
    """ Base skimmer module """
    def __init__(self, args):
        """ Constructor """
        super(SkimmerModule, self).__init__(args)

    def addArgs(self, parser):
        super(SkimmerModule, self).addArgs(parser)
        parser.add_argument("--keepOriginalBranches", action="store_true", help="Keep all original branches (in addition to those defined by the module)")
        parser.add_argument("--maxSelected", type=int, default=-1, help="Maximum number of accepted events (default: -1 for all)")
        parser.add_argument("--outputTreeName", type=str, default="Events", help="Name of the output tree")

    def initialize(self):
        """ initialize """
        if self.args.distributed == "worker" and len(self.inputs) == 0:
            raise RuntimeError("Worker task needs at least one input file")

    def interact(self):
        """ Interactively inspect a decorated input tree

        Available variables: ``tree`` (decorated tree), ``tup`` (raw tree),
        ``noSel`` (root selection), ``backend``, and ``runAndLS`` (e.g. ``(runExpr, lumiBlockExpr)``)
        (the inputs for the lumi mask), and ``op`` (:py:mod:`bamboo.treefunctions`).
        """
        tup, smpName, smpCfg = self.getATree()
        tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=smpName, sampleCfg=smpCfg)
        import bamboo.treefunctions as op
        import IPython
        IPython.embed()

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, sample=None, sampleCfg=None):
        """ Worker sequence: produce histograms from the input files

        More in detail, this will load the inputs, call :py:meth:`~bamboo.analysismodules.SkimmerModule.prepareTree`,
        add a lumi mask if requested, call :py:meth:`~bamboo.analysismodules.SkimmerModule.defineSkimSelection`,
        run over all files, and write the skimmed trees to the output file.
        """
        treeName = tree
        from .root import gbl
        tup = gbl.TChain(treeName)
        for fName in inputFiles:
            if not tup.Add(fName, 0):
                raise IOError("Could not open file {}".format(fName))
        tree, noSel, backend, runAndLS = self.prepareTree(tup, sample=sample, sampleCfg=sampleCfg)
        if certifiedLumiFile:
            noSel = addLumiMask(noSel, certifiedLumiFile, runRange=runRange, runAndLS=runAndLS)

        finalSel, brToKeep = self.defineSkimSelection(tree, noSel, sample=sample, sampleCfg=sampleCfg)
        defBr = dict((k,v) for k,v in brToKeep.items() if v is not None)
        origBr = list(k for k,v in brToKeep.items() if v is None)

        backend.writeSkim(finalSel, outputFile, self.args.outputTreeName, definedBranches=defBr, origBranchesToKeep=(None if self.args.keepOriginalBranches else origBr), maxSelected=self.args.maxSelected)

        outF = gbl.TFile.Open(outputFile, "UPDATE")
        self.mergeCounters(outF, inputFiles, sample=sample)
        outF.Close()
        return backend

    # processTrees customisation points
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        """ Create decorated tree, selection root (noSel), backend, and (run,LS) expressions

        should be implemented by concrete modules

        :param tree: decorated tree
        :param sample: sample name (as in the samples section of the analysis configuration file)
        :param sampleCfg: that sample's entry in the configuration file
        """
        return tree, None, None, None
    def defineSkimSelection(self, tree, noSel, sample=None, sampleCfg=None):
        """ Main method: define a selection for the skim

        should be implemented by concrete modules, and return a
        :py:class:`bamboo.plots.Selection` object

        :param tree: decorated tree
        :param noSel: base selection
        :param sample: sample name (as in the samples section of the analysis configuration file)
        :param sampleCfg: that sample's entry in the configuration file

        :returns: the skim :py:class:`bamboo.plots.Selection`, and a map ``{ name: expression }`` of branches to store (to store all the branches of the original tree in addition, pass --keepOriginalBranches to bambooRun; individual branches can be added by with an entry ``name: None`` entry)
        """
        return noSel, {}

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: write the equivalent analysis.yml file """
        pass ## TODO implement

class NanoAODSkimmerModule(NanoAODModule, SkimmerModule):
    """ A :py:class:`~bamboo.analysismodules.SkimmerModule` implementation for NanoAOD, adding decorations and merging of the counters """
    def __init__(self, args):
        super(NanoAODSkimmerModule, self).__init__(args)

class DataDrivenContribution:
    """
    Configuration helper class for data-driven contributions

    An instance is constructed for each contribution in any of the scenarios by
    the :py:meth:`bamboo.analysismodules.DataDrivenBackgroundAnalysisModule.initialize`
    method, with the name and configuration dictionary found in YAML file.
    The :py:meth:`~bamboo.analysismodules.DataDrivenContribution.usesSample`:,
    :py:meth:`~bamboo.analysismodules.DataDrivenContribution.replacesSample`: and
    :py:meth:`~bamboo.analysismodules.DataDrivenContribution.modifiedSampleConfig`:
    methods can be customised for other things than using the data samples
    to estimate a background contribution.
    """
    def __init__(self, name, config):
        self.name = name
        self.config = config
        self.uses = config.get("uses", [])
        self.replaces = config.get("replaces", [])
    def usesSample(self, sampleName, sampleConfig):
        """ Check if this contribution uses a sample (name or group in 'uses') """
        return sampleName in self.uses or sampleConfig.get("group") in self.uses
    def replacesSample(self, sampleName, sampleConfig):
        """ Check if this contribution replaces a sample (name or group in 'replaces') """
        return sampleName in self.replaces or sampleConfig.get("group") in self.replaces
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        """
        Construct the sample configuration for the reweighted counterpart of a sample

        The default implementation assumes a data sample and turns it into a MC sample
        (the luminosity is set as ``generated-events`` to avoid changing the normalisation).
        """
        modCfg = dict(sampleConfig)
        modCfg.update({
            "type": "mc",
            "generated-events": lumi,
            "cross-section": 1.,
            "group": self.name
            })
        return modCfg

class DataDrivenBackgroundAnalysisModule(AnalysisModule):
    """
    :py:class:`~bamboo.analysismoduldes.AnalysisModule` with support for data-driven backgrounds

    A number of contributions can be defined, each based on a list of samples or
    groups needed to evaluate the contribution (typically just data) and a list
    of samples or groups that should be left out when making the plot with
    data-driven contributions.
    The contributions should be defined in the analysis YAML file, with a block
    ``datadriven`` (at the top level) that could look as follows:

    .. code-block:: yaml

     datadriven:
       chargeMisID:
         uses: [ data ]
         replaces: [ DY ]
       nonprompt:
         uses: [ data ]
         replaces: [ TTbar ]

    The ``--datadriven`` command-line switch then allows to specify a scenario
    for data-driven backgrounds, i.e. a list of data-driven contributions to
    include (``all`` and ``none`` are also possible, the latter is the default
    setting).
    The parsed contributions are available as ``self.datadrivenContributions``,
    and the scenarios (each list is a list of contributions) as
    ``self.datadrivenScenarios``.
    """
    def addArgs(self, parser):
        super(DataDrivenBackgroundAnalysisModule, self).addArgs(parser)
        parser.add_argument("--datadriven", action="append", help="Scenarios for data-driven backgrounds ('all' for all available in yaml, 'none' for everything from MC (default), or a comma-separated list of contributions; several can be specified)")
    def initialize(self):
        super(DataDrivenBackgroundAnalysisModule, self).initialize()
        ddConfig = self.analysisConfig.get("datadriven", dict())
        if not self.args.datadriven:
            scenarios = [ [] ]
        elif not ddConfig:
            raise RuntimeError("--datadriven argument passed, but no 'datadriven' top-level block was found in the analysis YAML configuration file")
        else:
            scenarios = []
            for arg in self.args.datadriven:
                if arg.lower() == "all":
                    sc = tuple(sorted(ddConfig.keys()))
                    logger.info("--datadriven=all selected contributions {0}".format(", ".join(sc)))
                elif arg.lower() == "none":
                    sc = tuple()
                else:
                    sc = tuple(sorted(arg.split(",")))
                    if any(st not in ddConfig for st in sc):
                        raise RuntimeError("Unknown data-driven contribution(s): {0}".format(", ".join(st for st in sc if st not in ddConfig)))
                if sc not in scenarios:
                    scenarios.append(sc)
            if not scenarios:
                raise RuntimeError("No data-driven scenarios, please check the arguments ({0}) and config ({1})".format(self.args.datadriven, str(config)))
        self.datadrivenScenarios = scenarios
        self.datadrivenContributions = {
                contribName: DataDrivenContribution(contribName, config)
                for contribName, config in ddConfig.items()
                if any(contribName in scenario for scenario in scenarios)
                }

class DataDrivenBackgroundHistogramsModule(DataDrivenBackgroundAnalysisModule, HistogramsModule):
    """
    :py:class:`~bamboo.analysismoduldes.HistogramsModule` with support for data-driven backgrounds

    see the :py:class:`~bamboo.analysismoduldes.DataDrivenBackgroundAnalysisModule`
    class for more details about configuring data-driven backgrounds, and the
    :py:class:`~bamboo.plots.SelectionWithDataDriven` class for ensuring the
    necessary histograms are filled correctly.
    This class makes sure the histograms for the data-driven contributions are
    written to different files, and runs ``plotIt`` for the different scenarios.
    """
    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, sample=None, sampleCfg=None):
        backend = super(DataDrivenBackgroundHistogramsModule, self).processTrees(inputFiles, outputFile, tree=tree, certifiedLumiFile=certifiedLumiFile, runRange=runRange, sample=sample, sampleCfg=sampleCfg)
        ddFiles = {}
        from .root import gbl
        from .plots import SelectionWithDataDriven
        for p in self.plotList:
            if isinstance(p.selection, SelectionWithDataDriven):
                for ddSuffix in p.selection.dd:
                    if p.selection.dd[ddSuffix] is not None:
                        if ddSuffix not in ddFiles:
                            ddFiles[ddSuffix] = gbl.TFile.Open("{1}{0}{2}".format(ddSuffix, *os.path.splitext(outputFile)), "RECREATE")
                        ddFiles[ddSuffix].cd()
                        for h in backend.getResults(p, key=(p.name, ddSuffix)):
                            h.Write()
        for ddF in ddFiles.values():
            ddF.Close()
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir)
        from .plots import Plot, DerivedPlot, CutFlowReport, SelectionWithDataDriven
        plotList_cutflowreport = [ ap for ap in self.plotList if isinstance(ap, CutFlowReport) ]
        plotList_plotIt = [ ap for ap in self.plotList if ( isinstance(ap, Plot) or isinstance(ap, DerivedPlot) ) and len(ap.binnings) == 1 ]
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(config, plotList_cutflowreport, resultsdir=resultsdir, readCounters=self.readCounters, eras=self.args.eras, verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            eraMode, eras = self.args.eras
            if eras is None:
                eras = list(config["eras"].keys())
            # - step 1: group per max-available scenario (combination of data-driven contributions - not the same as configured scenarios, since this one is per-plot)
            plotList_per_availscenario = defaultdict(list)
            for p in plotList_plotIt:
                if isinstance(p.selection, SelectionWithDataDriven):
                    avSc = tuple(p.selection.dd.keys())
                    plotList_per_availscenario[avSc].append(p)
                else:
                    plotList_per_availscenario[tuple()].append(p)
            usedButNotUsedConfigContribs = set(contrib for avSc in plotList_per_availscenario.keys() for contrib in avSc if avSc not in self.datadrivenContributions)
            if len(usedButNotUsedConfigContribs) > 0:
                logger.error("Data-driven contributions {0} were used for defining selections and plots, but are not defined in the configuration file (present there are {1}); the plots will not be produced".format(
                    ", ".join(usedButNotUsedConfigContribs),
                    ", ".join(self.datadrivenContributions.keys())
                    ))
            # - step 2: collect plots for actual scenarios (avoiding adding them twice)
            plots_per_scenario = defaultdict(dict)
            for sc in self.datadrivenScenarios:
                for avSc, avPlots in plotList_per_availscenario.items():
                    plots_per_scenario[tuple(dd for dd in sc if dd in avSc)].update({p.name: p for p in avPlots})
            for sc, scPlots in plots_per_scenario.items():
                scName = "_".join(chain(["plots"], sc))
                cfgName = os.path.join(workdir, f"{scName}.yml")
                ## Modified samples: first remove replaced ones
                modSamples = { smp: smpCfg for smp, smpCfg in config["samples"].items()
                        if not any(self.datadrivenContributions[contrib].replacesSample(smp, smpCfg) for contrib in sc) }
                ## Then add data-driven contributions
                for contribName in sc:
                    contrib = self.datadrivenContributions[contribName]
                    for smp, smpCfg in config["samples"].items():
                        if contrib.usesSample(smp, smpCfg):
                            if "era" in smpCfg:
                                lumi = config["eras"][smpCfg["era"]]["luminosity"]
                            else:
                                lumi = sum(eraCfg["luminosity"] for eraCfg in config["eras"].values())
                            modSamples[f"{smp}{contribName}"] = contrib.modifiedSampleConfig(smp, smpCfg, lumi=lumi)
                config_sc = dict(config) # shallow copy
                config_sc["samples"] = modSamples
                ## TODO possible improvement: cache counters (e.g. by file name) - needs a change in readCounters interface
                writePlotIt(config_sc, list(scPlots.values()), cfgName, eras=(eraMode, eras), workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
                runPlotIt(cfgName, workdir=workdir, plotsdir=scName, plotIt=self.args.plotIt, eras=(eraMode, eras), verbose=self.args.verbose)
