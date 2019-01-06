"""
Analysis module base classes
"""
import argparse
import logging
logger = logging.getLogger(__name__)
import os.path
from .analysisutils import addLumiMask, downloadCertifiedLumiFiles, parseAnalysisConfig, readEnvConfig, runPlotIt

def reproduceArgv(args, group):
    """ Reconstruct the module-specific arguments (to pass them to the worker processes later on) """
    assert isinstance(group, argparse._ArgumentGroup)
    argv = []
    for action in group._group_actions:
        if isinstance(action, argparse._StoreTrueAction):
            if getattr(args, action.dest):
                argv.append(action.option_strings[0])
        elif isinstance(action, argparse._StoreAction):
            argv.append(action.option_strings[0])
            argv.append(getattr(args, action.dest))
        else:
            raise RuntimeError("Reconstruction of action {0} not supported".format(action))
    return argv

class AnalysisModule(object):
    """
    Base analysis module
    
    construct from an arguments list and provide a run method (called by bambooRun)
    """
    def __init__(self, args):
        parser = argparse.ArgumentParser(description=(
            "Run an analysis, i.e. process the samples in an analysis description file with a module (subclass of bamboo.analysismodules.AnalysisModule). "
            "There are three modes, specified by the --distributed option: if unspecified, one program processes all samples and collects the results; "
            "--distributed=driver does the same, but launches worker tasks "
            "(the same program with --distributed=worker, therefore some of the options only apply to 'driver' or 'worker' mode, or have a different interpretation) "
            "to process the samples, for instance on a batch system (depending on the settings in the --envConfig file)."))
        parser.add_argument("-m", "--module", type=str, default="bamboo.analysismodules:AnalysisModule", help="Module to run (format: modulenameOrPath[:classname])")
        parser.add_argument("-v", "--verbose", action="store_true", help="Run in verbose mode")
        parser.add_argument("--distributed", type=str, help="Role in distributed mode (sequential mode if not specified)", choices=["worker", "driver"])
        parser.add_argument("input", nargs="*", help="Input: analysis description yml file (driver mode) or files to process (worker mode)")
        parser.add_argument("-o", "--output", type=str, default=".", help="Output directory (driver mode) or file (worker mode) name")
        parser.add_argument("--interactive", "-i", action="store_true", help="Interactive mode (initialize to an IPython shell for exploration)")
        driver = parser.add_argument_group("driver mode only (--distributed=driver or unspecified) optional arguments")
        driver.add_argument("--redodbqueries", action="store_true", help="Redo all DAS/SAMADhi queries even if results can be read from cache files")
        driver.add_argument("--overwritesamplefilelists", action="store_true", help="Write DAS/SAMADhi results to files even if files exist (meaningless without --redodbqueries)")
        driver.add_argument("--envConfig", type=str, help="Config file to read computing environment configuration from (batch system, storage site etc.)")
        driver.add_argument("--plotIt", type=str, default="plotIt", help="plotIt executable to use (default is taken from $PATH)")
        worker = parser.add_argument_group("worker mode only (--distributed=worker) arguments")
        worker.add_argument("--treeName", type=str, default="Events", help="Tree name (default: Events)")
        worker.add_argument("--runRange", type=(lambda x : tuple(int(t.strip()) for t in x.split(","))), help="Run range (format: 'firstRun,lastRun')")
        worker.add_argument("--certifiedLumiFile", type=str, help="(local) path of a certified lumi JSON file")
        specific = parser.add_argument_group("module-specific arguments")
        self.addArgs(specific)
        self.args = parser.parse_args(args)
        self.specificArgv = reproduceArgv(self.args, specific)
        self.initialize()
    def addArgs(self, parser):
        """ Customize the ArgumentParser (set description and add arguments, if needed) """
        pass
    def initialize(self):
        """ Do more initialization that depends on command-line arguments """
        pass
    def getATree(self):
        if self.args.distributed == "worker":
            import ROOT
            tup = ROOT.TChain(self.args.treeName)
            tup.Add(self.args.input[0])
            return tup
        elif ( not self.args.distributed ) or self.args.distributed == "driver":
            if len(self.args.input) != 1:
                raise RuntimeError("Main process (driver or non-distributed) needs exactly one argument (analysis description YAML file)")
            anaCfgName = self.args.input[0]
            analysisCfg = parseAnalysisConfig(anaCfgName, redodbqueries=self.args.redodbqueries, overwritesamplefilelists=self.args.overwritesamplefilelists)
            import ROOT
            tup = ROOT.TChain(analysisCfg.get("tree", "Events"))
            tup.Add(next(analysisCfg["samples"].values())["files"][0])
            return tup
        else:
            raise RuntimeError("--distributed should be either worker, driver, or be unspecified (for sequential mode)")
    def run(self):
        """ Main method """
        if self.args.interactive:
            self.interact()
        else:
            if self.args.distributed == "worker":
                if ( not self.args.output.endswith(".root") ) or os.path.isdir(self.args.output):
                    raise RuntimeError("Output for worker processes needs to be a ROOT file")
                logger.info("Worker process: calling processTrees for {mod} with ({0}, {1}, treeName={tree}, certifiedLumiFile={certifiedLumiFile}, runRange={runRange}".format(self.args.input, self.args.output, mod=self.args.module, treeName=self.args.treeName, certifiedLumiFile=self.args.certifiedLumiFile, runRange=self.args.runRange))
                self.processTrees(self.args.input, self.args.output, treeName=self.args.treeName, certifiedLumiFile=self.args.certifiedLumiFile, runRange=self.args.runRange)
            elif ( not self.args.distributed ) or self.args.distributed == "driver":
                if len(self.args.input) != 1:
                    raise RuntimeError("Main process (driver or non-distributed) needs exactly one argument (analysis description YAML file)")
                anaCfgName = self.args.input[0]
                workdir = self.args.output
                envConfig = readEnvConfig(self.args.envConfig)
                analysisCfg = parseAnalysisConfig(anaCfgName, redodbqueries=self.args.redodbqueries, overwritesamplefilelists=self.args.overwritesamplefilelists, envConfig=envConfig)
                taskArgs = self.getTasks(analysisCfg, tree=analysisCfg.get("tree", "Events"))
                taskArgs, certifLumiFiles = downloadCertifiedLumiFiles(taskArgs, workdir=workdir)
                resultsdir = os.path.join(workdir, "results")
                if os.path.exists(resultsdir):
                    logger.warning("Output directory {0} exists, previous results may be overwritten".format(resultsdir))
                os.makedirs(resultsdir)
                ##
                if not self.args.distributed: ## sequential mode
                    for (inputs, output), kwargs in taskArgs:
                        output = os.path.join(resultsdir, output)
                        logger.info("Sequential mode: calling processTrees for {mod} with ({0}, {1}, {2}".format(inputs, output, ", ".join("{0}={1}".format(k,v) for k,v in kwargs.items()), mod=self.args.module))
                        self.processTrees(inputs, output, **kwargs)
                else:
                    from .batch import splitTask
                    backend = envConfig["batch"]["backend"]
                    tasks = [ splitTask(["bambooRun", "--module={0}".format(self.args.module), "--distributed=worker", "--output={0}".format(output)]+self.specificArgv+
                                        ["--{0}={1}".format(key, value) for key, value in kwargs], inputs, outdir=resultsdir, config=envConfig.get("splitting"))
                                for (inputs, output), kwargs in taskArgs ]
                    if backend == "slurm":
                        from . import batch_slurm as batchBackend
                        backendOpts = {
                                "sbatch_time"     : "0-00:20",
                                "sbatch_mem"      : "2048",
                                "stageoutFiles"   : ["*.root"],
                                "sbatch_workdir"  : os.getcwd(),
                                "sbatch_additionalOptions" : [ "--export=ALL" ],
                                }
                    elif backend == "htcondor":
                        from . import batch_htcondor as batchBackend
                        backendOpts = {
                                "cmd" : [
                                    "universe     = vanilla",
                                    "+MaxRuntime  = {0:d}".format(20*60), # 20 minutes
                                    "getenv       = True"
                                    ]
                                }
                    else:
                        raise RuntimeError("Unknown backend: {0}".format(backend))
                    batchworkdir = os.path.join(workdir, "batch")
                    if os.path.exists(batchworkdir):
                        raise RuntimeError("Directory '{0}' already exists, previous results would be overwritten".format(resultsdir))
                    os.makedirs(batchworkdir)
                    clusJobs = batchBackend.jobsFromTasks(tasks, workdir=batchworkdir, batchConfig=envConfig.get(backend), configOpts=backendOpts)
                    for j in clusJobs:
                        j.submit()
                    clusMon = batchBackend.makeTasksMonitor(clusJobs, tasks, interval=120)
                    clusMon.collect() ## wait for batch jobs to finish and finalize
                self.postProcess(taskArgs, config=analysisCfg, workdir=workdir, resultsdir=resultsdir)
            else:
                raise RuntimeError("--distributed should be either worker, driver, or be unspecified (for sequential mode)")

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None):
        """ worker method """
        pass
    def getTasks(self, analysisCfg, **extraOpts):
        """ Get tasks from args (for driver or sequential mode) """
        tasks = []
        for sName, sConfig in analysisCfg["samples"].items():
            opts = {
                "certifiedLumiFile" : sConfig.get("certified_lumi_file"),
                "runRange"          : sConfig.get("run_range")
                }
            opts.update(extraOpts)
            tasks.append(((sConfig["files"], "{0}.root".format(sName)), opts))
        return tasks

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Do postprocessing on the results of the tasks, if needed """
        pass
    def interact(self):
        """ Interactive mode (load some things and embed IPython) """
        pass ## define things and embed IPython

class HistogramsModule(AnalysisModule):
    """ Base histogram analysis module """
    def __init__(self, args):
        super(HistogramsModule, self).__init__(args)
        self.systVars = []
        self.plotList = []
    def initialize(self):
        if self.args.distributed == "worker" and len(self.args.input) == 0:
            raise RuntimeError("Worker task needs at least one input file")

    def interact(self):
        tup = self.getATree()
        tree, noSel, backend, (runExpr, lumiBlockExpr) = self.prepareTree(tup)
        import IPython
        IPython.embed()

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None):
        """ Worker sequence: open inputs, define histograms, fill them and save to output file """
        import ROOT
        tup = ROOT.TChain(tree)
        for fName in inputFiles:
            tup.Add(fName)
        tree, noSel, backend, runAndLS = self.prepareTree(tup)
        if certifiedLumiFile:
            noSel = addLumiMask(noSel, certifiedLumiFile, runRange=runRange, runAndLS=runAndLS)

        outF = ROOT.TFile.Open(outputFile, "RECREATE")
        plotList = self.definePlots(tree, noSel, systVar="nominal")
        if not self.plotList:
            self.plotList = list(plotList) ## shallow copy
        for systVar in self.systVars:
            plotList += self.definePlots(tree, noSel, systVar=systVar)

        outF.cd()
        for p in plotList:
            backend.getPlotResult(p).Write()
        outF.Close()
    # processTrees customisation points
    def prepareTree(self, tree):
        """ Create decorated tree, selection root (noSel), backend, and (run,LS) expressions """
        return tree, None, None, None
    def definePlots(self, tree, systVar="nominal"):
        """ Main method: define plots on the trees """
        return None, [] ## backend, and plot list

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        if not self.plotList: ## get plots if not already done so
            tup = self.getATree()
            tree, noSel, backend, runAndLS = self.prepareTree(tup)
            self.plotList = self.definePlots(tree, noSel, systVar="nominal")
        runPlotIt(config, self.plotList, workdir=workdir, resultsdir=resultsdir, plotIt=self.args.plotIt)
