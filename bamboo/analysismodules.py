"""
Analysis module base classes
"""
import logging
logger = logging.getLogger(__name__)

class AnalysisModule(object):
    """
    Base analysis module
    
    construct from an arguments list and provide a run method (called by bambooRun)
    """
    def __init__(self, args):
        import argparse
        parser = argparse.ArgumentParser(description=self.__doc__, prog="bambooRun --module={0}:{1}".format(__name__, self.__class__.__name__), add_help=False)
        parser.add_argument("--module-help", action="store_true", help="show this help message and exit")
        parser.add_argument("--module", type=str, help="Module specification")
        parser.add_argument("--distributed", type=str, help="Role in distributed mode (worker or driver; sequential if not specified)")
        parser.add_argument("-o", "--output", type=str, default="out.root", help="Output file name")
        parser.add_argument("input", nargs="*")
        parser.add_argument("--tree", type=str, default="Events", help="Tree name")
        parser.add_argument("--interactive", "-i", action="store_true", help="Interactive mode (initialize to an IPython shell for exploration)")
        parser.add_argument("--redodbqueries", action="store_true", help="Redo all DAS/SAMADhi queries even if results can be read from cache files")
        parser.add_argument("--overwritesamplefilelists", action="store_true", help="Write DAS/SAMADhi results to files even if files exist (meaningless without --redodbqueries)")
        self.addArgs(parser)
        self.args = parser.parse_args(args)
        if self.args.module_help:
            parser.print_help()
            import sys; sys.exit(0)
        self.initialize()
    def addArgs(self, parser):
        """ Customize the ArgumentParser (set description and add arguments, if needed) """
        pass
    def initialize(self):
        """ Do more initialization that depends on command-line arguments """
        pass
    def run(self):
        """ Main method """
        if self.args.interactive:
            self.interact()
        else:
            if self.args.distributed == "worker":
                self.processTrees(self.args.input, self.args.output)
            elif ( not self.args.distributed ) or self.args.distributed == "driver":
                if len(self.args.input) != 1:
                    raise RuntimeError("Main process (driver or non-distributed) needs exactly one argument (analysis descriptin YAML file)")
                anaCfgName = self.args.input[0]
                analysisCfg = self.parseAnalysisConfig(anaCfgName)
                import IPython, sys; IPython.embed(), sys.exit(1)
                inOutList = self.getTasks()
                if not self.args.distributed: ## sequential mode
                    for inputs, output in inOutList:
                        self.processTrees(inputs, output)
                else:
                    pass ## TODO figure out how to distribute them (local/slurm/...), splitting etc. and configure that
                self.postProcess(inOutList)
    def processTrees(self, inputFiles, outputFile):
        """ worker method """
        pass
    def parseAnalysisConfig(self, anaCfgName):
        import os.path
        import copy
        cfgDir = os.path.dirname(os.path.abspath(anaCfgName))
        import yaml
        with open(anaCfgName) as anaCfgF:
            analysisCfg = yaml.load(anaCfgF)
        ## finish loading samples (file lists)
        samples = dict()
        for smpName, smpCfg in analysisCfg["samples"].items():
            smp = copy.deepcopy(smpCfg)
            ## read cache, if it's there
            listfile, cachelist = None, []
            if "files" in smpCfg and str(smpCfg["files"]) == smpCfg["files"]:
                listfile = smpCfg["files"] if os.path.isabs(smpCfg["files"]) else os.path.join(cfgDir, smpCfg["files"])
                if os.path.isfile(listfile):
                    with open(listfile) as smpF:
                        cachelist = [ fn for fn in [ ln.strip() for ln in smpF ] if len(fn) > 0 ]

            if "db" in smpCfg and ( "files" not in smpCfg or len(cachelist) == 0 or self.args.redodbqueries ):
                if ":" not in smpCfg["db"]:
                    raise RuntimeError("'db' entry should be of the format 'protocol:location', e.g. 'das:/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD'")
                protocol, dbLoc = smpCfg["db"].split(":")
                files = []
                if protocol == "das":
                    dasQuery = "file dataset={0}".format(dbLoc)
                    files = [ fn for fn in [ ln.strip() for ln in subprocess.check_output(["dasgoclient", "-query", dasQuery]).split() ] if len(fn) > 0 ]
                    if len(files) == 0:
                        raise RuntimeError("No files found with DAS query {0}".format(dasQuery))
                elif protocol == "samadhi":
                    logger.warning("SAMADhi queries are not implemented yet")
                else:
                    raise RuntimeError("Unsupported protocol in '{0}': {1}".format(smpCfg["db"], protocol))
                smp["files"] = files
                if listfile and ( len(cachelist) == 0 or self.args.overwritesamplefilelists ):
                    with open(listfile, "w") as listF:
                        listF.writelines(files)
            elif "files" not in smpCfg:
                raise RuntimeError("Cannot load files for {0}: neither 'db' nor 'files' specified".format(smpName))
            elif listfile:
                if len(cachelist) == 0:
                    raise RuntimeError("No file names read from {0}".format())
                smp["files"] = cachelist
            else: ## list in yml
                smp["files"] = smpCfg["files"]
            samples[smpName] = smp
        analysisCfg["samples"] = samples
        return analysisCfg
    def getTasks(self):
        """ Get tasks from args (for driver or sequential mode) """
        return []
    def postProcess(self, taskList):
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

    def initialize(self):
        if self.args.distributed == "worker" and len(self.args.input) == 0:
            raise RuntimeError("Worker task needs at least one input file")

    def interact(self):
        if self.args.distributed == "worker":
            import ROOT
            tup = ROOT.TChain(self.args.tree)
            tup.Add(self.args.input[0])
            tree, noSel, backend = self.prepareTree(tup)
        import IPython
        IPython.embed()

    def processTrees(self, inputFiles, outputFile):
        """ Worker sequence: open inputs, define histograms, fill them and save to output file """
        import ROOT
        tup = ROOT.TChain(self.args.tree)
        for fName in inputFiles:
            tup.Add(fName)
        tree, noSel, backend = self.prepareTree(tup)

        outF = ROOT.TFile.Open(outputFile, "RECREATE")
        plotList = self.definePlots(tree, noSel, systVar="nominal")
        for systVar in self.systVars:
            plotList += self.definePlots(tree, noSel, systVar=systVar)

        outF.cd()
        for p in plotList:
            backend.getPlotResult(p).Write()
        outF.Close()
    # processTrees customisation points
    def prepareTree(self, tree):
        """ Create decorated tree, selection root (noSel) and backend"""
        return tree, None, None
    def definePlots(self, tree, systVar="nominal"):
        """ Main method: define plots on the trees """
        return None, [] ## backend, and plot list

    def postprocess(self, taskList):
        pass ## TODO write plotit and run it, and configure "postOnly"
