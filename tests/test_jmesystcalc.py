import pytest
import os.path

testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def nanojetargsMC16():
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016.root"))
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
    seed = (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0)
    genjet_pt   = RVec_float(tup.GenJet_pt, tup.nJet)
    genjet_eta  = RVec_float(tup.GenJet_eta, tup.nJet)
    genjet_phi  = RVec_float(tup.GenJet_phi, tup.nJet)
    genjet_mass = RVec_float(tup.GenJet_mass, tup.nJet)
    yield (jet_pt, jet_eta, jet_phi, jet_mass,
           jet_rawFactor, jet_area, tup.fixedGridRhoFastjetAll,
           seed, genjet_pt, genjet_eta, genjet_phi, genjet_mass)

@pytest.fixture(scope="module")
def nanojetargsMC16_postvalues():
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape.root"))
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
    seed = (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0)
    genjet_pt   = RVec_float(tup.GenJet_pt, tup.nJet)
    genjet_eta  = RVec_float(tup.GenJet_eta, tup.nJet)
    genjet_phi  = RVec_float(tup.GenJet_phi, tup.nJet)
    genjet_mass = RVec_float(tup.GenJet_mass, tup.nJet)
    ## results to compare to
    jet_vars = {
        "nominal" : (RVec_float(tup.Jet_pt_nom    , tup.nJet), RVec_float(tup.Jet_mass_nom    , tup.nJet)),
        "jerup"   : (RVec_float(tup.Jet_pt_jerUp  , tup.nJet), RVec_float(tup.Jet_mass_jerUp  , tup.nJet)),
        "jerdown" : (RVec_float(tup.Jet_pt_jerDown, tup.nJet), RVec_float(tup.Jet_mass_jerDown, tup.nJet))
        }
    from itertools import chain
    jet_vars.update(dict(chain.from_iterable(
        { "jes{0}up".format(src) : tuple(RVec_float(getattr(tup, "Jet_{0}_jes{1}Up".format(ivar, src)), tup.nJet) for ivar in ("pt", "mass")),
          "jes{0}down".format(src) : tuple(RVec_float(getattr(tup, "Jet_{0}_jes{1}Down".format(ivar, src)), tup.nJet) for ivar in ("pt", "mass")),
        }.items() for src in ("AbsoluteStat", "AbsoluteScale")
        )))
    yield ((jet_pt, jet_eta, jet_phi, jet_mass,
            jet_rawFactor, jet_area, tup.fixedGridRhoFastjetAll,
            seed, genjet_pt, genjet_eta, genjet_phi, genjet_mass),
           jet_vars
          )

@pytest.fixture(scope="module")
def nanoMETargsMC16_postvalues():
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape.root"))
    tup = f.Get("Events")
    tup.GetEntry(0)
    i = 0
    while tup.nJet < 5:
        i += 1
        tup.GetEntry(i)
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    args = ([ RVec_float(getattr(tup, "Jet_{0}".format(vNm)), tup.nJet)
            for vNm in ("pt", "eta", "phi", "mass", "rawFactor", "area", "muonSubtrFactor", "neEmEF", "chEmEF") ]
        + [ tup.fixedGridRhoFastjetAll, (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0) ]
        + [ RVec_float(getattr(tup, "GenJet_{0}".format(vNm)), tup.nGenJet) for vNm in ("pt", "eta", "phi", "mass") ]
        + [ tup.RawMET_phi, tup.RawMET_pt, tup.MET_MetUnclustEnUpDeltaX, tup.MET_MetUnclustEnUpDeltaY ]
        + [ RVec_float(getattr(tup, "CorrT1METJet_{0}".format(vNm)), tup.nJet)
            for vNm in ("rawPt", "eta", "phi", "area", "muonSubtrFactor") ] + [ RVec_float(), RVec_float() ]
        )
    met_vars = {
        "nominal" : (tup.MET_pt_nom, tup.MET_phi_nom),
        "jer"     : (tup.MET_pt_jer, tup.MET_phi_jer)
        }
    from itertools import chain
    met_vars.update(dict(
        ("{0}{1}".format(nm, var.lower()), (getattr(tup, "MET_pt_{0}{1}".format(nm, var)), getattr(tup, "MET_phi_{0}{1}".format(nm, var))))
        for var in ("Up", "Down") for nm in ["jer", "unclustEn"]+[ "jes{0}".format(jsnm) for jsnm in ("AbsoluteScale", "AbsoluteStat") ]
        ))
    yield tuple(args), met_vars

@pytest.fixture(scope="module")
def jmesystcalc_empty():
    from bamboo.root import gbl, loadJMESystematicsCalculators
    import bamboo.treefunctions
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    yield calc

@pytest.fixture(scope="module")
def jmesystcalcMC16_smear():
    from bamboo.root import gbl, loadJMESystematicsCalculators
    import bamboo.treefunctions
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    calc.setSmearing(os.path.join(testData, "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"), os.path.join(testData, "Summer16_25nsV1_MC_SF_AK4PFchs.txt"), True, 0.2, 3.)
    yield calc

@pytest.fixture(scope="module")
def jmesystcalcMC16_jec():
    from bamboo.root import gbl, loadJMESystematicsCalculators
    import bamboo.treefunctions
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt"))
    l2Param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt"))
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setJEC(jecParams)
    calc.setSmearing(os.path.join(testData, "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"), os.path.join(testData, "Summer16_25nsV1_MC_SF_AK4PFchs.txt"), True, 0.2, 3.)
    yield calc

@pytest.fixture(scope="module")
def jmesystcalcMC16_jesunc():
    from bamboo.root import gbl, loadJMESystematicsCalculators
    import bamboo.treefunctions
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt"))
    l2Param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt"))
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setJEC(jecParams)
    calc.setSmearing(os.path.join(testData, "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"), os.path.join(testData, "Summer16_25nsV1_MC_SF_AK4PFchs.txt"), True, 0.2, 3.)
    ## uncertaintysources?
    for jus in ["AbsoluteStat", "AbsoluteScale"]:
        param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt"), jus)
        calc.addJESUncertainty(jus, param)
    yield calc

@pytest.fixture(scope="module")
def metvarcalcMC16_jesunc():
    from bamboo.root import gbl, loadJMESystematicsCalculators
    import bamboo.treefunctions
    loadJMESystematicsCalculators()
    calc = gbl.Type1METVariationsCalculator()
    calc.setUnclusteredEnergyTreshold(15.)
    jecParams_L1 = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt"))
    l2Param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt"))
    jecParams_L1.push_back(l1Param)
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setL1JEC(jecParams_L1)
    calc.setJEC(jecParams)
    calc.setSmearing(os.path.join(testData, "Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"), os.path.join(testData, "Summer16_25nsV1_MC_SF_AK4PFchs.txt"), True, 0.2, 3.)
    ## uncertaintysources?
    for jus in ["AbsoluteStat", "AbsoluteScale"]:
        param = gbl.JetCorrectorParameters(os.path.join(testData, "Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt"), jus)
        calc.addJESUncertainty(jus, param)
    yield calc

def test_jmesystcalc_empty(jmesystcalc_empty):
    assert jmesystcalc_empty

def test_jmesystcalcMC16_smear(jmesystcalcMC16_smear):
    assert jmesystcalcMC16_smear

def test_jmesystcalcMC16_nano_smear(jmesystcalcMC16_smear, nanojetargsMC16):
    res = jmesystcalcMC16_smear.produceModifiedCollections(*nanojetargsMC16)
    assert res

def test_jmesystcalcMC16_nano_jec(jmesystcalcMC16_jec, nanojetargsMC16):
    res = jmesystcalcMC16_jec.produceModifiedCollections(*nanojetargsMC16)
    assert res

def test_jmesystcalcMC16_nano_jesunc(jmesystcalcMC16_jesunc, nanojetargsMC16):
    res = jmesystcalcMC16_jesunc.produceModifiedCollections(*nanojetargsMC16)
    assert res

import math
def isclose_float(a, b, tol=1.):
    from bamboo.root import gbl
    return math.isclose(a, b, rel_tol=tol*getattr(gbl, "std::numeric_limits<float>").epsilon())

def test_jmesystcalc_nanopost_jesunc(jmesystcalcMC16_jesunc, nanojetargsMC16_postvalues):
    nanojetargsMC16, postValues = nanojetargsMC16_postvalues
    res = jmesystcalcMC16_jesunc.produceModifiedCollections(*nanojetargsMC16)
    for ky,(post_pt, post_mass) in postValues.items():
        print(ky, res.at(ky).pt(), post_pt)
        print(ky, res.at(ky).mass(), post_mass)
        assert all(isclose_float(a,b, tol=2.) for a,b in zip(post_pt, res.at(ky).pt())) and all(isclose_float(a,b, tol=2.) for a,b in zip(post_mass, res.at(ky).mass()))

def test_metvarcalc_nanopost_jesunc(metvarcalcMC16_jesunc, nanoMETargsMC16_postvalues):
    nanoMETargsMC16, postValues = nanoMETargsMC16_postvalues
    res = metvarcalcMC16_jesunc.produceModifiedCollections(*nanoMETargsMC16)
    names = list(metvarcalcMC16_jesunc.availableProducts())
    for ky,(post_pt, post_phi) in postValues.items():
        idx = names.index(ky)
        print(ky, res.pt(idx), post_pt)
        print(ky, res.phi(idx), post_phi)
        assert isclose_float(res.pt(idx), post_pt, tol=3.) and isclose_float(res.phi(idx), post_phi, tol=3.)
