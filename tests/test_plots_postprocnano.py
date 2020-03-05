import pytest

import os.path
testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def decoNano():
    from bamboo import treefunctions as op
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape_puWeight.root"))
    tree = f.Get("Events")
    from bamboo.treedecorators import decorateNanoAOD, nanoPUWeightVar, nanoReadJetMETVar
    from bamboo.dataframebackend import DataframeBackend
    tup = decorateNanoAOD(tree, isMC=True, systVariations=[nanoPUWeightVar, nanoReadJetMETVar])
    be, noSel = DataframeBackend.create(tup)
    yield tup, noSel, be

def test(decoNano):
    # a somewhat realistic (but not very sensible) combination of selections and plots
    tup, noSel, be = decoNano
    noSel = noSel.refine("mcWeight", weight=[ tup.genWeight, tup.puWeight ])
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
    plots.append(Plot.make1D("btagSF", op.rng_product(cleanedJets, lambda j : j.btagSF), noSel, EquidistantBinning(100, 0., 1.2), title="Product of cleaned jets b-tag SF", xTitle="prod(btagSF)"))
    plots.append(Plot.make1D("btagSF_shape", op.rng_product(cleanedJets, lambda j : j.btagSF_shape), noSel, EquidistantBinning(100, 0., 1.2), title="Product of cleaned jets b-tag shape SF", xTitle="prod(btagSF_shape)"))
    cleanedJetsByDeepFlav = op.sort(cleanedJets, lambda jet: jet.btagDeepFlavB)
    hasMuJ = hasMuon.refine("hasMuonJ", cut=(op.rng_len(cleanedJets) > 0), weight=op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB))
    plots.append(Plot.make1D("hasMuonJ_prodBTags", op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB), hasMuJ, EquidistantBinning(1, 0., 1.), title="Product of jet b-tags", xTitle="X"))
    plots.append(Plot.make1D("MET", tup.MET.pt, noSel, EquidistantBinning(50, 0., 100.), title="MET pt", xTitle="MET"))
    plots.append(Plot.make1D("hasMuonJ_MET", tup.MET.pt, hasMuJ, EquidistantBinning(50, 0., 100.), title="MET pt", xTitle="MET"))
    addsyst5p = op.systematic(op.c_float(1.), name="addSyst", up=op.c_float(1.05), down=op.c_float(0.95))
    hasTwoJets = noSel.refine("hasTwoJets", cut=(op.rng_len(cleanedJets) > 1), weight=addsyst5p)
    plots.append(Plot.make1D("2J_cleanedProdBRegCorr", op.rng_product(cleanedJetsByDeepFlav[:2], lambda j : j.bRegCorr), hasTwoJets, EquidistantBinning(50, 0.5, 1.5)))

    assert all(h for p in plots for h in be.getPlotResults(p))
