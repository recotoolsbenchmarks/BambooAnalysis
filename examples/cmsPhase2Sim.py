""" Example CMS Phase2 simulation analysis module """

#############################################
##                                         ##
## methods and classes for a common module ##
##                                         ##
#############################################

from bamboo.analysismodules import AnalysisModule, HistogramsModule

class CMSPhase2SimModule(AnalysisModule):
    """ Base module for processing Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimModule, self).__init__(args)
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import decorateCMSPhase2SimTree
        from bamboo.dataframebackend import DataframeBackend
        t = decorateCMSPhase2SimTree(tree, isMC=True)
        be, noSel = DataframeBackend.create(t)
        return t, noSel, be, tuple()
    def mergeCounters(self, outF, infileNames, sample=None):
        from bamboo.root import gbl
        outF.cd()
        hNEvts = gbl.TH1F("nEvents", "Number of events", 1, 0., 1.)
        for fn in infileNames:
            f = gbl.TFile.Open(fn)
            tup = f.Get("myana/mytree") ## FIXME
            if not tup:
                raise RuntimeError("File {0} does not have a tree {1}".format(fn, self.args.treeName))
            hNEvts.Fill(0, tup.GetEntries())
        outF.cd()
        hNEvts.Write("nEvents")
    def readCounters(self, resultsFile):
        hNEvts = resultsFile.Get("nEvents")
        return {"nEvents": hNEvts.GetBinContent(1)}

class CMSPhase2SimHistoModule(CMSPhase2SimModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimHistoModule, self).__init__(args)


################################
##                            ##
## The actual analysis module ##
##                            ##
################################

class CMSPhase2SimTest(CMSPhase2SimHistoModule):
    """ Example plotter module for Phase2 flat trees """
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        plots = []

        electrons = op.select(t.elec, lambda el : op.AND(el.pt > 20., op.abs(el.eta) < 2.5, el.idpass, el.isopass))
        muons = op.select(t.muon, lambda mu : op.AND(mu.pt > 20., op.abs(mu.eta) < 2.4, mu.idpass, mu.isopass))
        cleanedJets = op.select(t.jet, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda el : op.deltaR(el.p4, j.p4) < 0.4)), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(mu.p4, j.p4) < 0.4))))

        plots.append(Plot.make1D("nElectron", op.rng_len(electrons), noSel, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuon", op.rng_len(muons), noSel, EqB(10, 0., 10.), title="Number of muons"))
        plots.append(Plot.make1D("nJet", op.rng_len(cleanedJets), noSel, EqB(10, 0., 10.), title="Number of jets"))

        hasOneJet = noSel.refine("hasOneJet", cut=(op.rng_len(cleanedJets) > 0))
        plots.append(Plot.make1D("jet1PT", cleanedJets[0].pt, hasOneJet, EqB(50, 0., 250.), title="Leading jet PT"))

        plots.append(Plot.make1D("gen_nDaug", op.map(t.genpart, lambda gp : op.rng_len(gp.children)), noSel, EqB(100, 0., 100.), title="Number of children"))

        g_tquarks = op.select(t.genpart, lambda gp : op.abs(gp.pid) == 6)
        plots.append(Plot.make1D("gen_nTop", op.rng_len(g_tquarks), noSel, EqB(100, 0., 100.), title="Number of top quark gen-particles"))
        g_dzero = op.select(t.genpart, lambda gp : op.abs(gp.pid) == 421)
        g_kshort = op.select(t.genpart, lambda gp : op.abs(gp.pid) == 310)
        plots += [
            Plot.make1D("gen_nD0S", op.rng_len(g_dzero), noSel, EqB(20, 0., 20.), title="Number of D0 gen-particles"),
            Plot.make1D("gen_D0ndaugs", op.map(g_dzero, lambda gp : op.rng_len(gp.children)), noSel, EqB(20, 0., 20.), title="D0 number of daughters"),
            Plot.make1D("gen_D0ndesc", op.map(g_dzero, lambda gp : op.rng_count(gp.descendants)), noSel, EqB(20, 0., 20.), title="D0 number of descendants"),
            Plot.make2D("gen_D0ndaugsndesc", (op.map(g_dzero, lambda gp : op.rng_len(gp.children)), op.map(g_dzero, lambda gp : op.rng_count(gp.descendants))), noSel, (EqB(20, 0., 20.), EqB(20, 0., 20.)), title="D0 number of daughters vs descendants"),
            Plot.make1D("gen_K0sndesc", op.map(g_kshort, lambda gp : op.rng_count(gp.descendants)), noSel, EqB(20, 0., 20.), title="K0_s number of descendants"),
            ]

        return plots
