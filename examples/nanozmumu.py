"""
Example analysis module: make a dimuon mass plot from a NanoAOD
"""
from bamboo.analysismodules import NanoAODHistoModule, NanoAODSkimmerModule
import os.path

class NanoZMuMu(NanoAODHistoModule):
    """ Example module: Z->MuMu histograms from NanoAOD """
    def __init__(self, args):
        super(NanoZMuMu, self).__init__(args)
        self.calcToAdd += ["nJet", "nMuon"] ## will do Jet variations and Rochester correction

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        ## initializes tree._Jet.calc so should be called first (better: use super() instead)
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg)
        from bamboo.analysisutils import makePileupWeight, configureJets, configureRochesterCorrection
        isNotWorker = (self.args.distributed != "worker")
        era = sampleCfg["era"]
        if self.isMC(sample):
            jecTag = None
            smearTag = None
            puWeightsFile = None
            if era == "2016":
                jecTag = "Summer16_07Aug2017_V20_MC"
                smearTag = "Summer16_25nsV1_MC"
                puWeightsFile = os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "data", "puweights.json")

            configureJets(tree, "Jet", "AK4PFchs", jec=jecTag, smear=smearTag, jesUncertaintySources=["Total"], mayWriteCache=isNotWorker)

            mcWgts = [ tree.genWeight ]
            if puWeightsFile:
                mcWgts.append(makePileupWeight(puWeightsFile, tree.Pileup_nTrueInt, variation="Nominal",
                    nameHint="bamboo_puWeight{0}".format("".join(c for c in sample if c.isalnum()))))

            noSel = noSel.refine("mcWeight", weight=mcWgts)

        else: ## DATA
            if era == "2016":
                jecTag = None
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    jecTag = "Summer16_07Aug2017BCD_V11_DATA"
                elif "2016E" in sample or "2016F" in sample:
                    jecTag = "Summer16_07Aug2017EF_V11_DATA"
                elif "2016G" in sample or "2016H" in sample:
                    jecTag = "Summer16_07Aug2017GH_V11_DATA"

            configureJets(tree, "Jet", "AK4PFchs", jec=jecTag, mayWriteCache=isNotWorker)

        if era == "2016":
            configureRochesterCorrection(tree._Muon.calc, os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "data", "RoccoR2016.txt"))

        return tree,noSel,be,lumiArgs

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, EquidistantBinning
        from bamboo import treefunctions as op
        from bamboo.analysisutils import forceDefine

        plots = []

        ## calculate (corrected) muon 4-momenta before accessing them
        forceDefine(t._Muon.calcProd, noSel)

        muons = op.select(t.Muon, lambda mu : mu.p4.Pt() > 20.)

        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])
        plots.append(Plot.make1D("dimu_M", op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel,
                EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

        ## evaluate jets for all events passing twoMuSel
        ## more optimization will be needed with systematics etc.
        forceDefine(t._Jet.calcProd, twoMuSel)

        jets = op.select(t.Jet, lambda j : j.p4.Pt() > 20.)
        plots.append(Plot.make1D("nJets", op.rng_len(jets), twoMuSel,
                EquidistantBinning(10, 0., 10.), title="Number of jets"))

        twoMuTwoJetSel = twoMuSel.refine("twoMuonsTwoJets", cut=[ op.rng_len(jets) > 1 ])

        plots.append(Plot.make1D("leadJetPT", jets[0].p4.Pt(), twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="Leading jet PT"))
        plots.append(Plot.make1D("subleadJetPT", jets[1].p4.Pt(), twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="Subleading jet PT"))

        return plots

class SkimNanoZMuMu(NanoAODSkimmerModule):
    def __init__(self, args):
        super(SkimNanoZMuMu, self).__init__(args)
    def defineSkimSelection(self, tree, noSel, sample=None, sampleCfg=None):
        from bamboo import treefunctions as op
        muons = op.select(tree.Muon, lambda mu : op.AND(mu.p4.Pt() > 20., op.abs(mu.p4.Eta()) < 2.4))
        hasTwoMu = noSel.refine("hasTwoMu", cut=(op.rng_len(muons) >= 2))
        varsToKeep = {"nMuon": None, "Muon_eta": None, "Muon_pt": None} ## from input file
        varsToKeep["nSelMuons"] = op.static_cast("UInt_t", op.rng_len(muons)) ## TBranch doesn't accept size_t
        varsToKeep["selMu_miniPFRelIsoNeu"] = op.map(muons, lambda mu : mu.miniPFRelIso_all - mu.miniPFRelIso_chg)
        return hasTwoMu, varsToKeep
