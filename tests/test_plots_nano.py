import pytest
import logging
logging.basicConfig(level=logging.DEBUG)

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
    tup = decorateNanoAOD(tree, isMC=True, addCalculators=("nJet", "MET", "nMuon"))
    be, noSel = DataframeBackend.create(tup)
    from bamboo.analysisutils import configureJets, configureType1MET, configureRochesterCorrection
    configureJets(tup._Jet, "AK4PFchs", isMC=True, backend=be, uName="test_plots_nano")
    configureType1MET(tup._MET, isMC=True, backend=be, uName="test_plots_nano")
    configureRochesterCorrection(tup._Muon, os.path.join(testData, "RoccoR2016.txt"), isMC=True, backend=be, uName="test_plots_nano")
    yield tup, noSel, be

def test(decoNano):
    # a somewhat realistic (but not very sensible) combination of selections and plots
    tup, noSel, be = decoNano
    noSel = noSel.refine("mcWeight", weight=tup.genWeight)
    from bamboo import treefunctions as op
    from bamboo.plots import Plot, Selection
    from bamboo.plots import EquidistantBinning as EqBin
    from bamboo.analysisutils import forceDefine
    from functools import partial
    forceDefine(tup._Muon.calcProd, noSel)
    plots = []
    electrons = op.select(tup.Electron, lambda ele : op.AND(ele.cutBased_Sum16 >= 3, ele.pt > 15., op.abs(ele.eta) < 2.4))
    muons = op.select(tup.Muon, lambda mu : op.AND(mu.pt > 10., op.abs(mu.eta) < 2.4, mu.mediumId, mu.pfRelIso04_all < 0.15))
    plots.append(Plot.make1D("nElectrons", op.rng_len(electrons), noSel, EqBin(10, 0., 10.), title="Number of electrons", xTitle="N_{e}"))
    hasMuon = noSel.refine("hasMuon", cut=(op.rng_len(muons) > 0))
    plots.append(Plot.make1D("hasMuon_leadMuPT", muons[0].pt, hasMuon, EqBin(50, 0., 100.), title="Leading muon PT", xTitle="p_{T}(mu_{1})"))
    forceDefine(tup._Jet.calcProd, noSel)
    jets = op.select(tup.Jet, lambda j : op.AND(j.pt > 20., op.abs(j.eta) < 2.4, j.jetId & 2))
    cleanedJets = op.select(jets, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
    plots.append(Plot.make1D("nCleanedJets", op.rng_len(cleanedJets), noSel, EqBin(10, 0., 10.), title="Number of cleaned jets", xTitle="N_{j}"))
    cleanedJetsByDeepFlav = op.sort(cleanedJets, lambda jet: jet.btagDeepFlavB)
    hasMuJ = hasMuon.refine("hasMuonJ", cut=(op.rng_len(cleanedJets) > 0), weight=op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB))
    plots.append(Plot.make1D("hasMuonJ_prodBTags", op.rng_product(cleanedJetsByDeepFlav, lambda jet: jet.btagDeepB), hasMuJ, EqBin(1, 0., 1.), title="Product of jet b-tags", xTitle="X"))
    plots.append(Plot.make1D("cleanedjet_pt", op.map(cleanedJets, lambda j : j.pt), noSel, EqBin(30, 30., 730.), title="Jet p_{T} (GeV)"))
    muRecoilJets = op.select(cleanedJets, partial((lambda l,j : op.deltaR(l.p4, j.p4) > 0.7), muons[0]))
    hasMuRecJ = hasMuon.refine("hasMuonRecJ", cut=(op.rng_len(muRecoilJets) > 0))
    dijets = op.combine(cleanedJets, N=2, pred=lambda j1,j2 : op.deltaR(j1.p4, j2.p4) > 0.6, samePred=lambda j1,j2 : j1.pt > j2.pt)
    plots.append(Plot.make1D("nCleanediJets", op.rng_len(dijets), noSel, EqBin(50, 0., 50.), title="Number of cleaned dijets", xTitle="N_{jj}"))
    alljetpairs = op.combine(cleanedJets, N=2, samePred=lambda j1,j2 : j1.pt > j2.pt)
    minJetDR = op.rng_min(alljetpairs, lambda pair: op.deltaR(pair[0].p4, pair[1].p4))
    plots.append(Plot.make1D("minJetDR", minJetDR, noSel, EqBin(40, 0.4, 2.)))
    assert all(h for p in plots for h in be.getPlotResults(p))
