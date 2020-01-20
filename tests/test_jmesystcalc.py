import pytest
import os.path

testData = os.path.join(os.path.dirname(__file__), "data")

def toRVecFloat(values):
    from bamboo.root import gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    res = RVec_float(len(values), 0.)
    for i,val in enumerate(values):
        res[i] = val
    return res

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
    jet_pt   = toRVecFloat(tup.Jet_pt)
    jet_eta  = toRVecFloat(tup.Jet_eta)
    jet_phi  = toRVecFloat(tup.Jet_phi)
    jet_mass = toRVecFloat(tup.Jet_mass)
    jet_rawFactor = toRVecFloat(tup.Jet_rawFactor)
    jet_area = toRVecFloat(tup.Jet_area)
    seed = (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0)
    genjet_pt   = toRVecFloat(tup.GenJet_pt)
    genjet_eta  = toRVecFloat(tup.GenJet_eta)
    genjet_phi  = toRVecFloat(tup.GenJet_phi)
    genjet_mass = toRVecFloat(tup.GenJet_mass)
    yield (jet_pt, jet_eta, jet_phi, jet_mass,
           jet_rawFactor, jet_area, tup.fixedGridRhoFastjetAll,
           seed, genjet_pt, genjet_eta, genjet_phi, genjet_mass)

@pytest.fixture(scope="module")
def nanojetargsMC16_postvalues():
    from bamboo.root import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        jet_pt   = toRVecFloat(tup.Jet_pt)
        jet_eta  = toRVecFloat(tup.Jet_eta)
        jet_phi  = toRVecFloat(tup.Jet_phi)
        jet_mass = toRVecFloat(tup.Jet_mass)
        jet_rawFactor = toRVecFloat(tup.Jet_rawFactor)
        jet_area = toRVecFloat(tup.Jet_area)
        seed = (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0)
        genjet_pt   = toRVecFloat(tup.GenJet_pt)
        genjet_eta  = toRVecFloat(tup.GenJet_eta)
        genjet_phi  = toRVecFloat(tup.GenJet_phi)
        genjet_mass = toRVecFloat(tup.GenJet_mass)
        ## results to compare to
        jet_vars = {
            "nominal" : (toRVecFloat(tup.Jet_pt_nom    ), toRVecFloat(tup.Jet_mass_nom    )),
            "jerup"   : (toRVecFloat(tup.Jet_pt_jerUp  ), toRVecFloat(tup.Jet_mass_jerUp  )),
            "jerdown" : (toRVecFloat(tup.Jet_pt_jerDown), toRVecFloat(tup.Jet_mass_jerDown))
            }
        from itertools import chain
        jet_vars.update(dict(chain.from_iterable(
            { "jes{0}up".format(src) : tuple(toRVecFloat(getattr(tup, "Jet_{0}_jes{1}Up".format(ivar, src))) for ivar in ("pt", "mass")),
              "jes{0}down".format(src) : tuple(toRVecFloat(getattr(tup, "Jet_{0}_jes{1}Down".format(ivar, src))) for ivar in ("pt", "mass")),
            }.items() for src in ("AbsoluteStat", "AbsoluteScale")
            )))
        res.append(((jet_pt, jet_eta, jet_phi, jet_mass,
                     jet_rawFactor, jet_area, tup.fixedGridRhoFastjetAll,
                     seed, genjet_pt, genjet_eta, genjet_phi, genjet_mass),
                     jet_vars))
    yield res

@pytest.fixture(scope="module")
def nanoMETargsMC16_postvalues():
    from bamboo.root import gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        args = ([ toRVecFloat(getattr(tup, "Jet_{0}".format(vNm)))
                for vNm in ("pt", "eta", "phi", "mass", "rawFactor", "area", "muonSubtrFactor", "neEmEF", "chEmEF") ]
            + [ tup.fixedGridRhoFastjetAll, (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0) ]
            + [ toRVecFloat(getattr(tup, "GenJet_{0}".format(vNm))) for vNm in ("pt", "eta", "phi", "mass") ]
            + [ tup.RawMET_phi, tup.RawMET_pt, tup.MET_MetUnclustEnUpDeltaX, tup.MET_MetUnclustEnUpDeltaY ]
            + [ toRVecFloat(getattr(tup, "CorrT1METJet_{0}".format(vNm)))
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
        res.append((tuple(args), met_vars))
    yield res

@pytest.fixture(scope="module")
def nanoMETFixEE2017argsMC17_postvalues():
    from bamboo.root import gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2017postproc_JMEKin_METFixEE2017.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        args = ([ toRVecFloat(getattr(tup, "Jet_{0}".format(vNm)))
                for vNm in ("pt", "eta", "phi", "mass", "rawFactor", "area", "muonSubtrFactor", "neEmEF", "chEmEF") ]
            + [ tup.fixedGridRhoFastjetAll, (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0) ]
            + [ toRVecFloat(getattr(tup, "GenJet_{0}".format(vNm))) for vNm in ("pt", "eta", "phi", "mass") ]
            + [ tup.RawMET_phi, tup.RawMET_pt, tup.METFixEE2017_MetUnclustEnUpDeltaX, tup.METFixEE2017_MetUnclustEnUpDeltaY ]
            + [ toRVecFloat(getattr(tup, "CorrT1METJet_{0}".format(vNm)))
                for vNm in ("rawPt", "eta", "phi", "area", "muonSubtrFactor") ] + [ RVec_float(), RVec_float() ]
            + [ tup.MET_phi, tup.MET_pt, tup.METFixEE2017_phi, tup.METFixEE2017_pt ]
            )
        met_vars = {
            "nominal" : (tup.METFixEE2017_pt_nom, tup.METFixEE2017_phi_nom),
            "jer"     : (tup.METFixEE2017_pt_jer, tup.METFixEE2017_phi_jer)
            }
        from itertools import chain
        met_vars.update(dict(
            ("{0}{1}".format(nm, var.lower()), (getattr(tup, "METFixEE2017_pt_{0}{1}".format(nm, var)), getattr(tup, "METFixEE2017_phi_{0}{1}".format(nm, var))))
            for var in ("Up", "Down") for nm in ["jer", "unclustEn"]+[ "jes{0}".format(jsnm) for jsnm in ("AbsoluteScale", "AbsoluteStat") ]
            ))
        res.append((tuple(args), met_vars))
    yield res

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
    jecTag = "Summer16_07Aug2017_V11_MC"
    jerTag = "Summer16_25nsV1_MC"
    jetType = "AK4PFchs"
    jecParams_L1 = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters(os.path.join(testData, "{0}_L1FastJet_{1}.txt".format(jecTag, jetType)))
    l2Param = gbl.JetCorrectorParameters(os.path.join(testData, "{0}_L2Relative_{1}.txt".format(jecTag, jetType)))
    jecParams_L1.push_back(l1Param)
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setL1JEC(jecParams_L1)
    calc.setJEC(jecParams)
    calc.setSmearing(os.path.join(testData, "{0}_PtResolution_{1}.txt".format(jerTag, jetType)), os.path.join(testData, "{0}_SF_{1}.txt".format(jerTag, jetType)), True, 0.2, 3.)
    ## uncertaintysources?
    for jus in ["AbsoluteStat", "AbsoluteScale"]:
        param = gbl.JetCorrectorParameters(os.path.join(testData, "{0}_UncertaintySources_{1}.txt".format(jecTag, jetType)), jus)
        calc.addJESUncertainty(jus, param)
    yield calc

@pytest.fixture(scope="module")
def metvarcalcMC17_FixEE():
    ## a better test even would be to have a different JEC than production
    from bamboo.root import gbl, loadJMESystematicsCalculators
    import bamboo.treefunctions
    loadJMESystematicsCalculators()
    calc = gbl.FixEE2017Type1METVariationsCalculator()
    calc.setUnclusteredEnergyTreshold(15.)
    jecTag = "Fall17_17Nov2017_V32_MC"
    jerTag = "Fall17_V3_MC"
    jetType = "AK4PFchs"
    jecParams_L1 = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
    l1Param = gbl.JetCorrectorParameters(os.path.join(testData, "{0}_L1FastJet_{1}.txt".format(jecTag, jetType)))
    l2Param = gbl.JetCorrectorParameters(os.path.join(testData, "{0}_L2Relative_{1}.txt".format(jecTag, jetType)))
    jecParams_L1.push_back(l1Param)
    jecParams.push_back(l1Param)
    jecParams.push_back(l2Param)
    calc.setL1JEC(jecParams_L1)
    calc.setJEC(jecParams)
    calc.setSmearing(os.path.join(testData, "{0}_PtResolution_{1}.txt".format(jerTag, jetType)), os.path.join(testData, "{0}_SF_{1}.txt".format(jerTag, jetType)), True, 0.2, 3.)
    ## uncertaintysources?
    for jus in ["AbsoluteStat", "AbsoluteScale"]:
        param = gbl.JetCorrectorParameters(os.path.join(testData, "{0}_UncertaintySources_{1}.txt".format(jecTag, jetType)), jus)
        calc.addJESUncertainty(jus, param)
    yield calc

def test_jmesystcalc_empty(jmesystcalc_empty):
    assert jmesystcalc_empty

def test_jmesystcalcMC16_smear(jmesystcalcMC16_smear):
    assert jmesystcalcMC16_smear

def test_jmesystcalcMC16_nano_smear(jmesystcalcMC16_smear, nanojetargsMC16):
    res = jmesystcalcMC16_smear.produce(*nanojetargsMC16)
    assert res

def test_jmesystcalcMC16_nano_jec(jmesystcalcMC16_jec, nanojetargsMC16):
    res = jmesystcalcMC16_jec.produce(*nanojetargsMC16)
    assert res

def test_jmesystcalcMC16_nano_jesunc(jmesystcalcMC16_jesunc, nanojetargsMC16):
    res = jmesystcalcMC16_jesunc.produce(*nanojetargsMC16)
    assert res

import math
def isclose_float(a, b, tol=1.):
    from bamboo.root import gbl
    return math.isclose(a, b, rel_tol=tol*getattr(gbl, "std::numeric_limits<float>").epsilon())

def test_jmesystcalc_nanopost_jesunc(jmesystcalcMC16_jesunc, nanojetargsMC16_postvalues):
    for nanojetargsMC16, postValues in nanojetargsMC16_postvalues:
        res = jmesystcalcMC16_jesunc.produce(*nanojetargsMC16)
        names = list(jmesystcalcMC16_jesunc.available())
        for ky,(post_pt, post_mass) in postValues.items():
            idx = names.index(ky)
            print(ky, res.pt(idx), post_pt)
            print(ky, res.mass(idx), post_mass)
            assert all(isclose_float(a,b) for a,b in zip(post_pt, res.pt(idx))) and all(isclose_float(a,b) for a,b in zip(post_mass, res.mass(idx)))

def test_metvarcalc_nanopost_jesunc(metvarcalcMC16_jesunc, nanoMETargsMC16_postvalues):
    for nanoMETargsMC16, postValues in nanoMETargsMC16_postvalues:
        res = metvarcalcMC16_jesunc.produce(*nanoMETargsMC16)
        names = list(metvarcalcMC16_jesunc.available())
        for ky,(post_pt, post_phi) in postValues.items():
            idx = names.index(ky)
            print(ky, res.pt(idx), post_pt)
            print(ky, res.phi(idx), post_phi)
            assert math.isclose(res.pt(idx), post_pt, rel_tol=1.e-6) and math.isclose(res.phi(idx), post_phi, rel_tol=1.e-6, abs_tol=1.e-6)

def test_metvarcalc_nanopost_jesunc_MCFixEE2017(metvarcalcMC17_FixEE, nanoMETFixEE2017argsMC17_postvalues):
    for nanoMETargsMC17FixEE, postValues in nanoMETFixEE2017argsMC17_postvalues:
        res = metvarcalcMC17_FixEE.produce(*nanoMETargsMC17FixEE)
        names = list(metvarcalcMC17_FixEE.available())
        for ky,(post_pt, post_phi) in postValues.items():
            idx = names.index(ky)
            print(ky, res.pt(idx), post_pt)
            print(ky, res.phi(idx), post_phi)
            assert math.isclose(res.pt(idx), post_pt, rel_tol=2.e-6) and math.isclose(res.phi(idx), post_phi, rel_tol=2.e-6, abs_tol=2.e-6)
