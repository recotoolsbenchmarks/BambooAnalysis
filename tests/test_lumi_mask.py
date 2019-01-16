import pytest

@pytest.fixture(scope="module")
def json_lumimask():
    import os.path
    jsonName = os.path.join(os.path.dirname(__file__), "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt")
    runRange = (276831, 277420) ## 2016E
    import bamboo.treefunctions ## loads and includes, minimal alternative below
    from cppyy import gbl
    #ROOT.gSystem.Load("libBambooLumiMask.so")
    #ROOT.gROOT.ProcessLine('#include "../ext/include/LumiMask.h"')
    yield gbl.LumiMask.fromJSON(jsonName, *runRange)

def test_lumimask_excluded(json_lumimask):
    assert not json_lumimask.accept(283681, 10 ) ## excluded run
def test_lumimask_included(json_lumimask):
    assert     json_lumimask.accept(276948, 10 ) ## included run
def test_lumimask_first(json_lumimask):
    assert     json_lumimask.accept(276948, 1  ) ## first for a run
def test_lumimask_last(json_lumimask):
    assert     json_lumimask.accept(276834, 720) ## last for that run
def test_lumimask_later_excluded(json_lumimask):
    assert not json_lumimask.accept(277070, 310) ## excluded part
def test_lumimask_later_included(json_lumimask):
    assert     json_lumimask.accept(277070, 311) ## later block
