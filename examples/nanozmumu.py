"""
Example analysis module: make a dimuon mass plot from a NanoAOD
"""
from bamboo.analysismodules import NanoAODHistoModule

class NanoZMuMu(NanoAODHistoModule):
    """ Example module: Z->MuMu histograms from NanoAOD """
    def __init__(self, args):
        super(NanoZMuMu, self).__init__(args)

    def prepareTree(self, tree, era=None, sample=None):
        ## initializes tree.Jet.calc so should be called first (better: use super() instead)
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, era=era, sample=sample)
        from bamboo.analysisutils import configureJets
        if era == "2016":
            if self.isMC(sample):
                configureJets(tree.Jet.calc, "AK4PFchs",
                    jec="Summer16_07Aug2017_V20_MC",
                    smear="Summer16_25nsV1_MC",
                    jesUncertaintySources=["Total"])
            else:
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    configureJets(tree.Jet.calc, "AK4PFchs",
                        jec="Summer16_07Aug2017BCD_V11_DATA")
                elif "2016E" in sample or "2016F" in sample:
                    configureJets(tree.Jet.calc, "AK4PFchs",
                        jec="Summer16_07Aug2017EF_V11_DATA")
                elif "2016G" in sample or "2016H" in sample:
                    configureJets(tree.Jet.calc, "AK4PFchs",
                        jec="Summer16_07Aug2017GH_V11_DATA")

        return tree,noSel,be,lumiArgs

    def definePlots(self, t, noSel, systVar="nominal", era=None, sample=None):
        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight)

        from bamboo.plots import Plot, EquidistantBinning
        from bamboo import treefunctions as op

        plots = []

        muons = op.select(t.Muon, lambda mu : mu.p4.Pt() > 20.)

        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])
        plots.append(Plot.make1D("dimu_M", op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel,
                EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

        ## evaluate jets for all events passing twoMuSel
        ## more optimization will be needed with systematics etc.
        from bamboo.analysisutils import forceDefine
        forceDefine(t.Jet.calcProd, twoMuSel)

        jets = op.select(t.Jet["nominal"], lambda j : j.p4.Pt() > 20.)
        plots.append(Plot.make1D("nJets", op.rng_len(jets), twoMuSel,
                EquidistantBinning(10, 0., 10.), title="Number of jets"))

        twoMuTwoJetSel = twoMuSel.refine("twoMuonsTwoJets", cut=[ op.rng_len(jets) > 1 ])

        plots.append(Plot.make1D("leadJetPT", jets[0].p4.Pt(), twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="Leading jet PT"))
        plots.append(Plot.make1D("subleadJetPT", jets[1].p4.Pt(), twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="Subleading jet PT"))

        return plots
