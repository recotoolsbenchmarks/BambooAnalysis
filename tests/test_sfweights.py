import pytest
import os.path

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
    assert sf_leptonSingle.get(makeParameters(Pt=20., Eta=1.5), gbl.Nominal) != 1.

def test_puWeight_constructEval(puWeight):
    from cppyy import gbl
    assert puWeight.get(makeParameters(NumTrueInteractions=20.), gbl.Nominal) != 1.
