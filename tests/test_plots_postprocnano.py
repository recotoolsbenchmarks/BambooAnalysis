import pytest

import logging
logger = logging.getLogger(__name__)

import os.path
testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def decoNanoMC16():
    from bamboo import treefunctions as op
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape_puWeight.root"))
    tree = f.Get("Events")
    from bamboo.treedecorators import decorateNanoAOD, NanoAODDescription, nanoPUWeightVar, ReadJetMETVar
    from bamboo.dataframebackend import DataframeBackend
    nanoReadJetMETVar_MC = ReadJetMETVar("Jet", "MET", jetsExclVars=["raw", "jesHF", "jesHF_2018", "jesHF_2017", "jesHF_2016"], metNomName="jer", metExclVars=["raw", "nom"], bTaggers=["csvv2", "deepcsv", "deepjet", "cmva"], bTagWPs=["L", "M", "T", "shape"])
    tup = decorateNanoAOD(tree, NanoAODDescription.get("v5", year="2016", isMC=True, systVariations=[nanoPUWeightVar, nanoReadJetMETVar_MC]))
    be, noSel = DataframeBackend.create(tup)
    yield tup, noSel, be

@pytest.fixture(scope="module")
def decoNanoData16():
    from bamboo import treefunctions as op
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DoubleMuon_2016Gpostproc_JMEKin.root"))
    tree = f.Get("Events")
    from bamboo.treedecorators import decorateNanoAOD, NanoAODDescription, ReadJetMETVar
    from bamboo.dataframebackend import DataframeBackend
    tup = decorateNanoAOD(tree, NanoAODDescription.get("v5", year="2016", isMC=False, systVariations=[
        ReadJetMETVar("Jet", "MET", jetsNomName="nom", metNomName="nom")]))
    be, noSel = DataframeBackend.create(tup)
    yield tup, noSel, be

def definePlots(tup, noSel, isMC=False):
    # a somewhat realistic (but not very sensible) combination of selections and plots
    from bamboo import treefunctions as op
    from bamboo.plots import Plot, Selection, EquidistantBinning
    plots = []
    electrons = op.select(tup.Electron, lambda ele : op.AND(ele.cutBased_Sum16 >= 3, ele.pt > 15., op.abs(ele.eta) < 2.4))
    muons = op.select(tup.Muon, lambda mu : op.AND(mu.pt > 10., op.abs(mu.eta) < 2.4, mu.mediumId, mu.pfRelIso04_all < 0.15))
    plots.append(Plot.make1D("nElectrons", op.rng_len(electrons), noSel, EquidistantBinning(10, 0., 10.), title="Number of electrons", xTitle="N_{e}"))
    hasMuon = noSel.refine("hasMuon", cut=(op.rng_len(muons) > 0))
    plots.append(Plot.make1D("hasMuon_leadMuPT", muons[0].pt, hasMuon, EquidistantBinning(50, 0., 100.), title="Leading muon PT", xTitle="p_{T}(mu_{1})"))
    jets = op.select(tup.Jet, lambda j : op.AND(j.pt > 20., op.abs(j.eta) < 2.4, j.jetId & 2))
    cleanedJets = op.select(jets, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
    plots.append(Plot.make1D("nCleanedJets", op.rng_len(cleanedJets), noSel, EquidistantBinning(10, 0., 10.), title="Number of cleaned jets", xTitle="N_{j}"))
    if isMC:
        plots.append(Plot.make1D("btagSF", op.rng_product(cleanedJets, lambda j : j.btagSF), noSel, EquidistantBinning(100, 0., 1.2), title="Product of cleaned jets b-tag SF", xTitle="prod(btagSF)"))
        plots.append(Plot.make1D("btagSF_shape", op.rng_product(cleanedJets, lambda j : j.btagSF_shape), noSel, EquidistantBinning(100, 0., 1.2), title="Product of cleaned jets b-tag shape SF", xTitle="prod(btagSF_shape)"))
    cleanedJetsByDeepFlav = op.sort(cleanedJets, lambda jet: jet.btagDeepFlavB)
    hasMuJ = hasMuon.refine("hasMuonJ", cut=(op.rng_len(cleanedJets) > 0), weight=op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB))
    plots.append(Plot.make1D("hasMuonJ_prodBTags", op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB), hasMuJ, EquidistantBinning(1, 0., 1.), title="Product of jet b-tags", xTitle="X"))
    plots.append(Plot.make1D("MET", tup.MET.pt, noSel, EquidistantBinning(50, 0., 100.), title="MET pt", xTitle="MET"))
    plots.append(Plot.make1D("hasMuonJ_MET", tup.MET.pt, hasMuJ, EquidistantBinning(50, 0., 100.), title="MET pt", xTitle="MET"))
    if isMC:
        addsyst5p = op.systematic(op.c_float(1.), name="addSyst", up=op.c_float(1.05), down=op.c_float(0.95))
        hasTwoJets = noSel.refine("hasTwoJets", cut=(op.rng_len(cleanedJets) > 1), weight=addsyst5p)
        plots.append(Plot.make1D("2J_cleanedProdBRegCorr", op.rng_product(cleanedJetsByDeepFlav[:2], lambda j : j.bRegCorr), hasTwoJets, EquidistantBinning(50, 0.5, 1.5)))
    return plots

def testMC16(decoNanoMC16):
    tup, noSel, be = decoNanoMC16
    noSel = noSel.refine("mcWeight", weight=[ tup.genWeight, tup.puWeight ])
    plots = definePlots(tup, noSel, isMC=True)
    histos_per_plot = { p.name : list(be.getResults(p)) for p in plots }
    logger.debug("Plots for MC16")
    for pN, histos in histos_per_plot.items():
        logger.debug("Plot {0} -> {1}".format(pN, ", ".join(h.GetName() for h in histos)))
    assert all(h for p,histos in histos_per_plot.items() for h in histos)

from bamboo.root import gbl
@pytest.mark.skipif(int(gbl.gROOT.GetVersion().split("/")[0].split(".")[1]) < 18, reason="Test not supported for ROOT older than 6.18")
def testData16(decoNanoData16):
    tup, noSel, be = decoNanoData16
    plots = definePlots(tup, noSel)
    histos_per_plot = { p.name : list(be.getResults(p)) for p in plots }
    logger.debug("Plots for data16")
    for pN, histos in histos_per_plot.items():
        logger.debug("Plot {0} -> {1}".format(pN, ", ".join(h.GetName() for h in histos)))
    assert all(len(histos) == 1 for histos in histos_per_plot.values())
    assert all(h for p,histos in histos_per_plot.items() for h in histos)
