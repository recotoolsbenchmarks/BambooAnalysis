import pytest

import os.path
testData = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="module")
def decoNano():
    from bamboo import treefunctions as op
    from cppyy import gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016.root"))
    tree = f.Get("Events")
    from bamboo.treedecorators import decorateNanoAOD
    from bamboo.dataframebackend import DataframeBackend
    tup = decorateNanoAOD(tree, isMC=True)
    be, noSel = DataframeBackend.create(tup)
    jetcalcName = be.symbol("JMESystematicsCalculator <<name>>{};", nameHint="bamboo_jmeSystCalc")
    tup._Jet.initCalc(op.extVar("JMESystematicsCalculator", jetcalcName), calcHandle=getattr(gbl, jetcalcName))
    roccorName = be.symbol("RochesterCorrectionCalculator <<name>>{};", nameHint="bamboo_roccorCalc")
    tup._Muon.initCalc(op.extVar("RochesterCorrectionCalculator", roccorName), calcHandle=getattr(gbl, roccorName))
    yield tup

def test_getSimpleObjects(decoNano):
    assert decoNano.Electron[0]
    assert decoNano.Muon[0]
    assert decoNano.Photon[0]
    assert decoNano.SubJet[0]
    assert decoNano.FatJet[0]

def test_getJet(decoNano):
    assert decoNano.Jet[0]

def test_simpleRef(decoNano):
    assert decoNano.Electron[0].photon
    assert decoNano.FatJet[0].subJet1

def test_toJetRef(decoNano):
    assert decoNano.Electron[0].jet
    assert decoNano.Muon[0].jet

def test_fromJetRef(decoNano):
    jet = decoNano.Jet[0]
    assert jet.muon1
    assert jet.muon2
    assert jet.electron1
    assert jet.electron2
