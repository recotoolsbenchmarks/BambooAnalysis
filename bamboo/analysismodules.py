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
import logging
logger = logging.getLogger(__name__)
import os.path
from .analysisutils import addLumiMask, downloadCertifiedLumiFiles, parseAnalysisConfig, readEnvConfig, runPlotIt

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

class AnalysisModule(object):
    """ Base analysis module
    
    Adds common infrastructure for parsing analysis config files
    and running on a batch system, with customization points for
    concrete classes to implement (most importantly
    :py:meth:`~bamboo.analysismodules.AnalysisModule.processTrees`
    and :py:meth:`~bamboo.analysismodules.AnalysisModule.postProcess`)
    """
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
        parser.add_argument("--distributed", type=str, help="Role in distributed mode (sequential mode if not specified)", choices=["worker", "driver"])
        parser.add_argument("input", nargs="*", help="Input: analysis description yml file (driver mode) or files to process (worker mode)")
        parser.add_argument("-o", "--output", type=str, default=".", help="Output directory (driver mode) or file (worker mode) name")
        parser.add_argument("--interactive", "-i", action="store_true", help="Interactive mode (initialize to an IPython shell for exploration)")
        driver = parser.add_argument_group("driver mode only (--distributed=driver or unspecified) optional arguments")
        driver.add_argument("--redodbqueries", action="store_true", help="Redo all DAS/SAMADhi queries even if results can be read from cache files")
        driver.add_argument("--overwritesamplefilelists", action="store_true", help="Write DAS/SAMADhi results to files even if files exist (meaningless without --redodbqueries)")
        driver.add_argument("--envConfig", type=str, help="Config file to read computing environment configuration from (batch system, storage site etc.)")
        driver.add_argument("--onlypost", action="store_true", help="Only run postprocessing step on previous results")
        worker = parser.add_argument_group("worker mode only (--distributed=worker) arguments")
        worker.add_argument("--treeName", type=str, default="Events", help="Tree name (default: Events)")
        worker.add_argument("--runRange", type=parseRunRange, help="Run range (format: 'firstRun,lastRun')")
        worker.add_argument("--certifiedLumiFile", type=str, help="(local) path of a certified lumi JSON file")
        worker.add_argument("--sample", type=str, help="Sample name (as in the samples section of the analysis configuration)")
        worker.add_argument("--era", type=str, help="Era of the sample")
        specific = parser.add_argument_group("module-specific arguments")
        self.addArgs(specific)
        self.args = parser.parse_args(args)
        self.specificArgv = reproduceArgv(self.args, specific)
        self.initialize()
    def addArgs(self, parser):
        """ Hook for adding module-specific argument parsing (receives an argument group), parsed arguments are available in ``self.args`` afterwards """
        pass
    def initialize(self):
        """ Hook for module-specific initialization (called from the constructor after parsing arguments) """
        pass
    def getATree(self):
        """ Retrieve a representative TTree, e.g. for defining the plots or interactive inspection, and a dictionary with metadata """
        if self.args.distributed == "worker":
            from cppyy import gbl
            tup = gbl.TChain(self.args.treeName)
            tup.Add(self.args.input[0])
            return tup, {}
        elif ( not self.args.distributed ) or self.args.distributed == "driver":
            if len(self.args.input) != 1:
                raise RuntimeError("Main process (driver or non-distributed) needs exactly one argument (analysis description YAML file)")
            anaCfgName = self.args.input[0]
            envConfig = readEnvConfig(self.args.envConfig)
            analysisCfg = parseAnalysisConfig(anaCfgName, redodbqueries=self.args.redodbqueries, overwritesamplefilelists=self.args.overwritesamplefilelists, envConfig=envConfig)
            from cppyy import gbl
            tup = gbl.TChain(analysisCfg.get("tree", "Events"))
            smpNm,smpCfg = next(itm for itm in analysisCfg["samples"].items())
            tup.Add(smpCfg["files"][0])
            return tup, {"name": smpNm, "era": smpCfg["era"]}
        else:
            raise RuntimeError("--distributed should be either worker, driver, or be unspecified (for sequential mode)")
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
                logger.info("Worker process: calling processTrees for {mod} with ({0}, {1}, treeName={treeName}, certifiedLumiFile={certifiedLumiFile}, runRange={runRange}, era={era}, sample={sample})".format(self.args.input, self.args.output, mod=self.args.module, treeName=self.args.treeName, certifiedLumiFile=self.args.certifiedLumiFile, runRange=self.args.runRange, era=self.args.era, sample=self.args.sample))
                self.processTrees(self.args.input, self.args.output, tree=self.args.treeName, certifiedLumiFile=self.args.certifiedLumiFile, runRange=self.args.runRange, era=self.args.era, sample=self.args.sample)
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
                if self.args.onlypost:
                    if not os.path.exists(resultsdir):
                        raise RuntimeError("Results directory {0} does not exist".format(resultsdir))
                    ## TODO check for all output files?
                else:
                    if os.path.exists(resultsdir):
                        logger.warning("Output directory {0} exists, previous results may be overwritten".format(resultsdir))
                    else:
                        os.makedirs(resultsdir)
                    ##
                    if not self.args.distributed: ## sequential mode
                        for (inputs, output), kwargs in taskArgs:
                            output = os.path.join(resultsdir, output)
                            logger.info("Sequential mode: calling processTrees for {mod} with ({0}, {1}, {2}".format(inputs, output, ", ".join("{0}={1}".format(k,v) for k,v in kwargs.items()), mod=self.args.module))
                            if "runRange" in kwargs:
                                kwargs["runRange"] = parseRunRange(kwargs["runRange"])
                            self.processTrees(inputs, output, **kwargs)
                    else:
                        from .batch import splitTask
                        backend = envConfig["batch"]["backend"]
                        tasks = [ splitTask(["bambooRun", "--module={0}".format(modAbsPath(self.args.module)), "--distributed=worker", "--output={0}".format(output)]+self.specificArgv+
                                            ["--{0}={1}".format(key, value) for key, value in kwargs.items()], inputs, outdir=resultsdir, config=envConfig.get("splitting"))
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
                        clusJobs = batchBackend.jobsFromTasks(tasks, workdir=os.path.join(workdir, "batch"), batchConfig=envConfig.get(backend), configOpts=backendOpts)
                        for j in clusJobs:
                            j.submit()
                        clusMon = batchBackend.makeTasksMonitor(clusJobs, tasks, interval=120)
                        clusMon.collect() ## wait for batch jobs to finish and finalize
                self.postProcess(taskArgs, config=analysisCfg, workdir=workdir, resultsdir=resultsdir)
            else:
                raise RuntimeError("--distributed should be either worker, driver, or be unspecified (for sequential mode)")

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, era=None, sample=None):
        """ worker method: produce results (e.g. histograms or trees) from the input files

        should be implemented by concrete modules

        :param inputFiles: input file names
        :param outputFile: output file name
        :param tree: key name of the tree inside the files
        :param certifiedLumiFile: lumi mask json file name
        :param runRange: run range to consider (for efficiency of the lumi mask)
        :param era: era to which the sample belongs (allows for era-based customization)
        :param sample: sample name (key in the samples block of the configuration file)
        """
        pass
    def getTasks(self, analysisCfg, **extraOpts):
        """ Get tasks from analysis configs (and args), called in for driver or sequential mode

        Should return a list of ``(inputs, output), kwargs``
        """
        tasks = []
        for sName, sConfig in analysisCfg["samples"].items():
            opts = dict(extraOpts)
            if "certified_lumi_file" in sConfig:
                opts["certifiedLumiFile"] = sConfig.get("certified_lumi_file")
            if "run_range" in sConfig:
                opts["runRange"] = ",".join(str(rn) for rn in sConfig.get("run_range"))
            opts["sample"] = sName
            if "era" in sConfig:
                opts["era"] = sConfig["era"]
            tasks.append(((sConfig["files"], "{0}.root".format(sName)), opts))
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

        Defines ``plotList`` and ``systVars`` member variables. The former will store a list of plots
        (only nominal), the second can be used by concrete classes to pass a list of systematic variations
        (:py:meth:`~bamboo.analysismodules.HistogramsModule.processTrees` will call
        :py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots` with each of them to produce all histograms)
        """
        super(HistogramsModule, self).__init__(args)
        self.systVars = []
        self.plotList = []

    def addArgs(self, parser):
        parser.add_argument("--plotIt", type=str, default="plotIt", help="plotIt executable to use (default is taken from $PATH)")

    def initialize(self):
        """ initialize """
        if self.args.distributed == "worker" and len(self.args.input) == 0:
            raise RuntimeError("Worker task needs at least one input file")

    def interact(self):
        """ Interactively inspect a decorated input tree

        Available variables: ``tree`` (decorated tree), ``tup`` (raw tree),
        ``noSel`` (root selection), ``backend``, ``runExpr`` and ``lumiBlockExpr``
        (the inputs for the lumi mask), and ``op`` (:py:mod:`bamboo.treefunctions`).
        """
        tup, tupInfo = self.getATree()
        tree, noSel, backend, (runExpr, lumiBlockExpr) = self.prepareTree(tup, era=tupInfo.get("era"), sample=tupInfo.get("name"))
        import bamboo.treefunctions as op
        import IPython
        IPython.embed()

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, era=None, sample=None):
        """ Worker sequence: produce histograms from the input files

        More in detail, this will load the inputs, call :py:meth:`~bamboo.analysismodules.HistogramsModule.prepareTree`,
        add a lumi mask if requested, call :py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots`
        (with ``systVar="nominal"`` and each of the variations defined in ``self.systVars``),
        run over all files, and write the produced histograms to the output file.
        """
        from cppyy import gbl
        tup = gbl.TChain(tree)
        for fName in inputFiles:
            tup.Add(fName)
        tree, noSel, backend, runAndLS = self.prepareTree(tup, era=era, sample=sample)
        if certifiedLumiFile:
            noSel = addLumiMask(noSel, certifiedLumiFile, runRange=runRange, runAndLS=runAndLS)

        outF = gbl.TFile.Open(outputFile, "RECREATE")
        plotList = self.definePlots(tree, noSel, systVar="nominal", era=era, sample=sample)
        if not self.plotList:
            self.plotList = list(plotList) ## shallow copy
        for systVar in self.systVars:
            plotList += self.definePlots(tree, noSel, systVar=systVar, era=era, sample=sample)

        outF.cd()
        for p in plotList:
            backend.getPlotResult(p).Write()
        self.mergeCounters(outF, inputFiles)
        outF.Close()
    # processTrees customisation points
    def prepareTree(self, tree, era=None, sample=None):
        """ Create decorated tree, selection root (noSel), backend, and (run,LS) expressions

        should be implemented by concrete modules

        :param tree: decorated tree
        :param era: era name, to allow era-based customization
        :param sample: sample name (as in the samples section of the analysis configuration file)
        """
        return tree, None, None, None
    def definePlots(self, tree, noSel, systVar="nominal", era=None, sample=None):
        """ Main method: define plots on the trees (for a give systematic variation)

        should be implemented by concrete modules, and return a list of
        :py:class:`bamboo.plots.Plot` objects.
        The structure (name, binning) of the histograms should not depend on the sample, era,
        or the systematic variation, and the list should be the same for all values, with
        one exception: for ``systVar`` different from ``"nominal"`` some histograms may
        be left out for some samples (e.g. for data none should be filled, and certain
        systematic effects may only affect some of the MC samples).

        :param tree: decorated tree
        :param noSel: base selection
        :param systVar: ``"nominal"``, or the name of a systematic variation
        :param era: era name, to allow era-based customization
        :param sample: sample name (as in the samples section of the analysis configuration file)
        """
        return [] ## plot list
    def mergeCounters(self, outF, infileNames):
        """ Merge counters

        should be implemented by concrete modules

        :param outF: output file (TFile pointer)
        :param infileNames: input file names
        """
        pass
    def readCounters(self, resultsFile):
        """ Read counters from results file

        should be implemented by concrete modules, and return a dictionary with
        counter names and the corresponding sums

        :param resultsFile: TFile pointer to the results file
        """
        return dict()

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: run plotIt

        The list of plots is created if needed (from a representative file,
        this enables rerunning the postprocessing step on the results files),
        and then plotIt is executed
        """
        if not self.plotList:
            tup, tupInfo = self.getATree()
            tree, noSel, backend, runAndLS = self.prepareTree(tup, era=tupInfo.get("era"), sample=tupInfo.get("name"))
            self.plotList = self.definePlots(tree, noSel, systVar="nominal", era=tupInfo.get("era"), sample=tupInfo.get("name"))
        for eraName, eraCfg in config.get("eras", {}).items():
            if eraName == "combined": ## TODO think harder what needs to go into plotit and what not
                for cmbCfg in eraCfg:
                    runPlotIt(config, self.plotList, workdir=workdir, resultsdir=resultsdir, plotIt=self.args.plotIt, readCounters=self.readCounters, era=cmbCfg)
            else:
                runPlotIt(config, self.plotList, workdir=workdir, resultsdir=resultsdir, plotIt=self.args.plotIt, readCounters=self.readCounters, era=eraName)

class NanoAODModule(AnalysisModule):
    """ A :py:class:`~bamboo.analysismodules.AnalysisModule` extension for NanoAOD, adding decorations and merging of the counters """
    def __init__(self, args):
        super(NanoAODModule, self).__init__(args)
    def isMC(self, sampleName):
        return not any(sampleName.startswith(pd) for pd in ("SingleMuon", "SingleElectron", "DoubleMuon", "DoubleEG", "MuonEG"))
    def prepareTree(self, tree, era=None, sample=None):
        """ Add NanoAOD decorations, and create an RDataFrame backend """
        from bamboo.treedecorators import decorateNanoAOD
        from bamboo.dataframebackend import DataframeBackend
        t = decorateNanoAOD(tree, isMC=self.isMC(sample))
        be, noSel = DataframeBackend.create(t)
        ## force definition of the jet variations calculator such that
        ## it can be configured by calling its member methods through t.Jet.calc
        from bamboo import treefunctions as op
        from cppyy import gbl
        jetcalcName = be.symbol("JMESystematicsCalculator <<name>>{{}}; // for {0}".format(sample), nameHint="bamboo_jmeSystCalc{0}".format("".join(c for c in sample if c.isalnum())))
        t.Jet.initCalc(op.extVar("JMESystematicsCalculator", jetcalcName), calcHandle=getattr(gbl, jetcalcName))
        return t, noSel, be, (t.run, t.luminosityBlock)
    def mergeCounters(self, outF, infileNames):
        """ Merge the ``Runs`` trees """
        from cppyy import gbl
        cruns = gbl.TChain("Runs")
        for fn in infileNames:
            cruns.Add(fn)
        outF.cd()
        runs = cruns.CloneTree()
        runs.Write("Runs")
    def readCounters(self, resultsFile):
        """ Sum over each leaf of the (merged) ``Runs`` tree (except ``run``) """
        runs = resultsFile.Get("Runs")
        from cppyy import gbl
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
                    for i in range(len(vals)):
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
        parser.add_argument("--maxSelected", type=int, default=-1, help="Maximum number of accepted events (default: -1 for all)")

    def initialize(self):
        """ initialize """
        if self.args.distributed == "worker" and len(self.args.input) == 0:
            raise RuntimeError("Worker task needs at least one input file")

    def interact(self):
        """ Inte ## TODO fully genericractively inspect a decorated input tree

        Available variables: ``tree`` (decorated tree), ``tup`` (raw tree),
        ``noSel`` (root selection), ``backend``, ``runExpr`` and ``lumiBlockExpr``
        (the inputs for the lumi mask), and ``op`` (:py:mod:`bamboo.treefunctions`).
        """
        tup, tupInfo = self.getATree()
        tree, noSel, backend, (runExpr, lumiBlockExpr) = self.prepareTree(tup, era=tupInfo.get("era"), sample=tupInfo.get("name"))
        import bamboo.treefunctions as op
        import IPython
        IPython.embed()

    def processTrees(self, inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None, era=None, sample=None):
        """ Worker sequence: produce histograms from the input files

        More in detail, this will load the inputs, call :py:meth:`~bamboo.analysismodules.SkimmerModule.prepareTree`,
        add a lumi mask if requested, call :py:meth:`~bamboo.analysismodules.SkimmerModule.defineSkimSelection`,
        run over all files, and write the skimmed trees to the output file.
        """
        treeName = tree
        from cppyy import gbl
        tup = gbl.TChain(treeName)
        for fName in inputFiles:
            tup.Add(fName)
        tree, noSel, backend, runAndLS = self.prepareTree(tup, era=era, sample=sample)
        if certifiedLumiFile:
            noSel = addLumiMask(noSel, certifiedLumiFile, runRange=runRange, runAndLS=runAndLS)

        finalSel, defBrToKeep = self.defineSkimSelection(tree, noSel, era=era, sample=sample)
        selDF = backend.selDFs[finalSel.name].df

        if self.args.maxSelected and self.args.maxSelected > 0:
            selDF = selDF.Range(self.args.maxSelected)

        # exclude defined columns
        allcolN = selDF.GetColumnNames()
        defcolN = selDF.GetDefinedColumnNames()
        origcolN = type(allcolN)()
        for cn in allcolN:
            if cn not in defcolN or cn in defBrToKeep:
                origcolN.push_back(cn)

        selDF.Snapshot(treeName, outputFile, origcolN)

        outF = gbl.TFile.Open(outputFile, "UPDATE")
        self.mergeCounters(outF, inputFiles)
        outF.Close()

    # processTrees customisation points
    def prepareTree(self, tree, era=None, sample=None):
        """ Create decorated tree, selection root (noSel), backend, and (run,LS) expressions

        should be implemented by concrete modules

        :param tree: decorated tree
        :param era: era name, to allow era-based customization
        :param sample: sample name (as in the samples section of the analysis configuration file)
        """
        return tree, None, None, None
    def defineSkimSelection(self, tree, noSel, era=None, sample=None):
        """ Main method: define a selection for the skim

        should be implemented by concrete modules, and return a
        :py:class:`bamboo.plots.Selection` object

        :param tree: decorated tree
        :param noSel: base selection
        :param era: era name, to allow era-based customization
        :param sample: sample name (as in the samples section of the analysis configuration file)

        :returns: the skim :py:class:`bamboo.plots.Selection`, and a list of defined branch names to save
        """
        return noSel, []

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: write the equivalent analysis.yml file """
        pass ## TODO implement

class NanoAODSkimmerModule(NanoAODModule, SkimmerModule):
    """ A :py:class:`~bamboo.analysismodules.SkimmerModule` implementation for NanoAOD, adding decorations and merging of the counters """
    def __init__(self, args):
        super(NanoAODSkimmerModule, self).__init__(args)
