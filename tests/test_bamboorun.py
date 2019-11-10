import pytest
import os.path
import subprocess
import shutil

testData = os.path.join(os.path.dirname(__file__), "data")
examples = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")

def test_bambooRun_help():
    assert subprocess.run(["bambooRun", "--help"]).returncode == 0

def test_bambooRun_help_nanozmm():
    assert subprocess.run(["bambooRun", "--module={0}:NanoZMuMu".format(os.path.join(examples, "nanozmumu.py")), "--help"]).returncode == 0

@pytest.mark.skipif(shutil.which("plotIt") is None, reason="plotIt was not found")
def test_bambooRun_nanozmm_1():
    assert subprocess.run(["bambooRun", "--module={0}:NanoZMuMu".format(os.path.join(examples, "nanozmumu.py")), os.path.join(examples, "test1.yml"), "--output=nanozmm_test1"]).returncode == 0

def test_bambooRun_nanozmm_worker():
    assert subprocess.run(["bambooRun", "--module={0}:NanoZMuMu".format(os.path.join(examples, "nanozmumu.py")), "--distributed=worker", "--sample=DY_2016", os.path.join(testData, "DY_M50_2016.root"), "--output=test_worker_1.root"]).returncode == 0
