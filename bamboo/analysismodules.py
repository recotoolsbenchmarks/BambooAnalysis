"""
Analysis module base classes
"""
import argparse

class AnalysisModule(object):
    """
    Base analysis module
    
    construct from an arguments list and provide a run method (called by bambooRun)
    """
    def __init__(self, args):
        parser = argparse.ArgumentParser(description=self.__doc__, prog="bambooRun --module={0}:{1}".format(__name__, self.__class__.__name__), add_help=False)
        parser.add_argument("--module-help", action="store_true", help="show this help message and exit")
        parser.add_argument("--module", type=str, help="Module specification")
        parser.add_argument("--worker", action="store_true", dest="isWorker", help="Worker task for distributed running")
        parser.add_argument("-o", "--output", type=str, default="out.root", help="Output file name")
        parser.add_argument("input", nargs="*")
        parser.add_argument("--tree", type=str, default="Events", help="Tree name")
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
        if self.args.isWorker:
            self.runWorker()
        else:
            self.runDriver()
    def runWorker(self):
        """ Worker method - to be implemented by concrete class (process input files and write output) """
        pass
    def runDriver(self):
        """ Driver method - to be implemented by concrete class (e.g. interpret input file, submit slurm jobs and combine) """
        pass

class HistogramsModule(AnalysisModule):
    """ Base histogram analysis module """
    def __init__(self, args):
        super(HistogramsModule, self).__init__(args)
    def initialize(self):
        if self.args.isWorker and len(self.args.input) == 0:
            raise RuntimeError("Worker task needs at least one input file")
    def runWorker(self):
        """ Worker sequence: open inputs, define histograms, fill them and save to output file """
        import ROOT
        tup = ROOT.TChain(self.args.tree)
        for fName in self.args.input:
            tup.Add(fName)
        outF = ROOT.TFile.Open(self.args.output, "RECREATE")
        be, plotList = self.definePlots(tup)
        outF.cd()
        for p in plotList:
            be.getPlotResult(p).Write()
        outF.Close()
    def runDriver(self):
        ## TODO implement (json(s), plotIt from definePlots with a placeholder tree)
        pass
    def definePlots(self, tree):
        """ Main method: define plots on the trees """
        return None, [] ## backend, and plot list
