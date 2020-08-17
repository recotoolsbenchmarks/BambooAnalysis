import pytest
import os.path
import ROOT ## 6.22: enable getattr(TTree, branchname)

testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def nanojetargs():
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016.root"))
    tup = f.Get("Events")
    tup.GetEntry(0)
    i = 0
    while tup.nMuon < 2:
        i += 1
        tup.GetEntry(i)
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    RVec_int   = getattr(gbl, "ROOT::VecOps::RVec<Int_t>")
    muon_pt   = RVec_float(tup.Muon_pt, tup.nMuon)
    muon_eta  = RVec_float(tup.Muon_eta, tup.nMuon)
    muon_phi  = RVec_float(tup.Muon_phi, tup.nMuon)
    muon_mass = RVec_float(tup.Muon_mass, tup.nMuon)
    muon_charge = RVec_int(tup.Muon_charge, tup.nMuon)
    muon_nLayers = RVec_int(tup.Muon_nTrackerLayers, tup.nMuon)
    muon_genIdx = RVec_int(tup.Muon_genPartIdx, tup.nMuon)
    gen_pt = RVec_float(tup.GenPart_pt, tup.nGenPart)
    seed = 5489
    yield ((muon_pt, muon_eta, muon_phi, muon_mass, muon_charge, muon_nLayers, muon_genIdx, gen_pt, seed),
           (muon_pt, muon_eta, muon_phi, muon_mass, muon_charge, muon_nLayers, RVec_int(), RVec_float(), seed))

@pytest.fixture(scope="module")
def roccorcalc_empty():
    from bamboo.root import gbl, loadRochesterCorrectionCalculator
    import bamboo.treefunctions
    loadRochesterCorrectionCalculator()
    calc = gbl.RochesterCorrectionCalculator()
    yield calc

@pytest.fixture(scope="module")
def roccorcalc_2016():
    from bamboo.root import gbl, loadRochesterCorrectionCalculator
    import bamboo.treefunctions
    loadRochesterCorrectionCalculator()
    calc = gbl.RochesterCorrectionCalculator()
    calc.setRochesterCorrection(os.path.join(testData, "RoccoR2016.txt"))
    yield calc

def test_rochester_empty(roccorcalc_empty):
    assert roccorcalc_empty

def test_rochester_2016(roccorcalc_empty):
    assert roccorcalc_2016

def test_rochester_nano_off_mc(roccorcalc_empty, nanojetargs):
    nanojet_mc, nanojet_data = nanojetargs
    res_mc = roccorcalc_empty.produce(*nanojet_mc)
    assert res_mc

def test_rochester_nano_off_data(roccorcalc_empty, nanojetargs):
    nanojet_mc, nanojet_data = nanojetargs
    res_data = roccorcalc_empty.produce(*nanojet_data)
    assert res_data

def test_rochester_nano_off_mc(roccorcalc_2016, nanojetargs):
    nanojet_mc, nanojet_data = nanojetargs
    res_mc = roccorcalc_2016.produce(*nanojet_mc)
    assert res_mc

def test_rochester_nano_off_data(roccorcalc_2016, nanojetargs):
    nanojet_mc, nanojet_data = nanojetargs
    res_data = roccorcalc_2016.produce(*nanojet_data)
    assert res_data
