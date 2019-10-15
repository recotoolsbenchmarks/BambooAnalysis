import pytest

import os.path
testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def decoNano():
    from bamboo import treefunctions as op
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016.root"))
    tree = f.Get("Events")
    from bamboo.treedecorators import decorateNanoAOD
    from bamboo.dataframebackend import DataframeBackend
    tup = decorateNanoAOD(tree, isMC=True, addCalculators=("nJet", "nMuon"))
    be, noSel = DataframeBackend.create(tup)
    from bamboo.analysisutils import addJMESystematicsCalculator, addRochesterCorrectionCalculator
    addJMESystematicsCalculator(be, tup._Jet, isMC=True)
    addRochesterCorrectionCalculator(be, tup._Muon, isMC=True)
    yield tup, noSel, be

def test(decoNano):
    # a somewhat realistic (but not very sensible) combination of selections and plots
    tup, noSel, be = decoNano
    from bamboo import treefunctions as op
    from bamboo.plots import Plot, Selection, EquidistantBinning
    from bamboo.analysisutils import forceDefine
    forceDefine(tup._Muon.calcProd, noSel)
    plots = []
    electrons = op.select(tup.Electron, lambda ele : op.AND(ele.cutBased_Sum16 >= 3, ele.pt > 15., op.abs(ele.eta) < 2.4))
    muons = op.select(tup.Muon, lambda mu : op.AND(mu.pt > 10., op.abs(mu.eta) < 2.4, mu.mediumId, mu.pfRelIso04_all < 0.15))
    plots.append(Plot.make1D("nElectrons", op.rng_len(electrons), noSel, EquidistantBinning(10, 0., 10.), title="Number of electrons", xTitle="N_{e}"))
    hasMuon = noSel.refine("hasMuon", cut=(op.rng_len(muons) > 0))
    plots.append(Plot.make1D("hasMuon_leadMuPT", muons[0].pt, hasMuon, EquidistantBinning(50, 0., 100.), title="Leading muon PT", xTitle="p_{T}(mu_{1})"))
    forceDefine(tup._Jet.calcProd, noSel)
    jets = op.select(tup.Jet, lambda j : op.AND(j.pt > 20., op.abs(j.eta) < 2.4, j.jetId & 2))
    cleanedJets = op.select(jets, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
    plots.append(Plot.make1D("nCleanedJets", op.rng_len(cleanedJets), noSel, EquidistantBinning(10, 0., 10.), title="Number of cleaned jets", xTitle="N_{j}"))
    cleanedJetsByDeepFlav = op.sort(cleanedJets, lambda jet: jet.btagDeepFlavB)
    hasMuJ = hasMuon.refine("hasMuonJ", cut=(op.rng_len(cleanedJets) > 0), weight=op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB))
    plots.append(Plot.make1D("hasMuonJ_prodBTags", op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB), hasMuJ, EquidistantBinning(1, 0., 1.), title="Product of jet b-tags", xTitle="X"))
    plots.append(Plot.make1D("cleanedjet_pt", op.map(cleanedJets, lambda j : j.pt), noSel, EquidistantBinning(30, 30., 730.), title="Jet p_{T} (GeV)"))
    assert all(h for p in plots for h in be.getPlotResults(p))
