import pytest

@pytest.fixture(scope="module")
def nanojetargs():
    from cppyy import gbl
    res_t = getattr(gbl, "JMESystematicsCalculator::result_t") ## trigger dictionary generation
    f = gbl.TFile.Open("/storage/data/cms/store/mc/RunIISummer16NanoAODv4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PUMoriond17_Nano14Dec2018_102X_mcRun2_asymptotic_v6_ext2-v1/40000/8B742F43-C711-5140-9999-44FA54B3A21F.root")
    tup = f.Get("Events")
    tup.GetEntry(0)
    i = 0
    while tup.nJet < 5:
        i += 1
        tup.GetEntry(i)
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    jet_pt   = RVec_float(tup.Jet_pt, tup.nJet)
    jet_eta  = RVec_float(tup.Jet_eta, tup.nJet)
    jet_phi  = RVec_float(tup.Jet_phi, tup.nJet)
    jet_mass = RVec_float(tup.Jet_mass, tup.nJet)
    jet_rawFactor = RVec_float(tup.Jet_rawFactor, tup.nJet)
    jet_area = RVec_float(tup.Jet_area, tup.nJet)
    genjet_pt   = RVec_float(tup.GenJet_pt, tup.nJet)
    genjet_eta  = RVec_float(tup.GenJet_eta, tup.nJet)
    genjet_phi  = RVec_float(tup.GenJet_phi, tup.nJet)
    genjet_mass = RVec_float(tup.GenJet_mass, tup.nJet)
    yield (jet_pt, jet_eta, jet_phi, jet_mass,
           jet_rawFactor, jet_area, tup.fixedGridRhoFastjetAll,
           tup.MET_phi, tup.MET_pt, tup.MET_sumEt,
           genjet_pt, genjet_eta, genjet_phi, genjet_mass)

@pytest.fixture(scope="module")
def jmesystcalc_empty():
    from cppyy import gbl
    import bamboo.treefunctions ## loads and includes
    calc = gbl.JMESystematicsCalculator()
    yield calc

@pytest.fixture(scope="module")
def jmesystcalc_smear():
    from cppyy import gbl
    import bamboo.treefunctions
    calc = gbl.JMESystematicsCalculator()
    calc.setSmearing("tests/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt", "tests/Summer16_25nsV1_MC_SF_AK4PFchs.txt", True, 0.2, 3.)
    yield calc

@pytest.fixture(scope="module")
def jmesystcalc_jec():
    from cppyy import gbl
    import bamboo.treefunctions
    calc = gbl.JMESystematicsCalculator()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters("tests/Summer16_07Aug2017_V20_MC_L1FastJet_AK4PFchs.txt")
    l2Param = gbl.JetCorrectorParameters("tests/Summer16_07Aug2017_V20_MC_L2Relative_AK4PFchs.txt")
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setJEC(jecParams)
    yield calc

@pytest.fixture(scope="module")
def jmesystcalc_jesunc():
    from cppyy import gbl
    import bamboo.treefunctions
    calc = gbl.JMESystematicsCalculator()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters("tests/Summer16_07Aug2017_V20_MC_L1FastJet_AK4PFchs.txt")
    l2Param = gbl.JetCorrectorParameters("tests/Summer16_07Aug2017_V20_MC_L2Relative_AK4PFchs.txt")
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setJEC(jecParams)
    ## uncertaintysources?
    for jus in ["AbsoluteStat", "AbsoluteScale"]:
        param = gbl.JetCorrectorParameters("tests/Summer16_07Aug2017_V20_MC_UncertaintySources_AK4PFchs.txt", jus)
        calc.addJESUncertainty(jus, param)
    yield calc

def test_jmesystcalc_empty(jmesystcalc_empty):
    assert jmesystcalc_empty

def test_jmesystcalc_smear(jmesystcalc_smear):
    assert jmesystcalc_smear

def test_jmesystcalc_nano_smear(jmesystcalc_smear, nanojetargs):
    res = jmesystcalc_smear.produceModifiedCollections(*nanojetargs)
    assert res
    # assert None ## to see the printout

def test_jmesystcalc_nano_jec(jmesystcalc_jec, nanojetargs):
    res = jmesystcalc_jec.produceModifiedCollections(*nanojetargs)
    assert res
    # assert None ## to see the printout

def test_jmesystcalc_nano_jesunc(jmesystcalc_jesunc, nanojetargs):
    res = jmesystcalc_jesunc.produceModifiedCollections(*nanojetargs)
    assert res
    # assert None ## to see the printout
