import pytest
import os.path
import math

def isclose_float(a, b):
    from cppyy import gbl
    return math.isclose(a, b, rel_tol=getattr(gbl, "std::numeric_limits<float>").epsilon())

def makeParameters(**kwargs):
    """ Construct Parameters argument for BinnedValues """ 
    import bamboo.treefunctions
    from cppyy import gbl
    gbl.BinningVariable ## somehow loads Pt, Eta etc.
    params = gbl.Parameters()
    for k,v in kwargs.items():
        params.set(getattr(gbl, k), v)
    return params

@pytest.fixture(scope="module")
def sf_leptonSingle():
    from cppyy import gbl
    import bamboo.treefunctions
    gbl.SystVariation ## somehow loads Nominal etc.
    elSFJSON = os.path.join(os.path.dirname(__file__), "data", "Electron_EGamma_SF2D_loose_moriond17.json")
    return gbl.ScaleFactor(elSFJSON)

@pytest.fixture(scope="module")
def puWeight():
    from cppyy import gbl
    import bamboo.treefunctions
    puWeightJSON = os.path.join(os.path.dirname(__file__), "data", "puweights.json")
    return gbl.ScaleFactor(puWeightJSON)

def test_lepSingle_constructEval(sf_leptonSingle):
    from cppyy import gbl
    assert isclose_float(sf_leptonSingle.get(makeParameters(Pt=20., Eta=1.5), gbl.Nominal), 0.9901639223098755)
    assert isclose_float(sf_leptonSingle.get(makeParameters(Pt=20., Eta=1.5), gbl.Up), 0.9901639223098755+0.19609383660010074)
    assert isclose_float(sf_leptonSingle.get(makeParameters(Pt=20., Eta=1.5), gbl.Down), 0.9901639223098755-0.19609383660010074)

def test_puWeight_constructEvalInRange(puWeight):
    from cppyy import gbl
    assert isclose_float(puWeight.get(makeParameters(NumTrueInteractions=20.5), gbl.Nominal), 1.0656023337493363)

def test_puWeight_evalBinEdge(puWeight):
    from cppyy import gbl
    assert isclose_float(puWeight.get(makeParameters(NumTrueInteractions=20.), gbl.Nominal), 1.0656023337493363)

def test_puWeight_evalOutOfRangeBelow(puWeight):
    from cppyy import gbl
    assert isclose_float(puWeight.get(makeParameters(NumTrueInteractions=-.5), gbl.Nominal), 0.36607730074755906)

def test_puWeight_evalOutOfRangeAbove(puWeight):
    from cppyy import gbl
    assert isclose_float(puWeight.get(makeParameters(NumTrueInteractions=100.), gbl.Nominal), 0.001723281482149061)
