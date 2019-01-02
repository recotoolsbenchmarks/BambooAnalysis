import os.path
from urllib.parse import urlparse
import subprocess
jsonPath = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
runRange = (276831, 277420) ## 2016E 
jsonName = urlparse(jsonPath).path.split("/")[-1]
if not os.path.exists(jsonName):
    print("Downloading JSON file from {0}".format(jsonPath))
    subprocess.check_call(["wget", jsonPath])
import bamboo.treefunctions ## loads and includes
import ROOT
#ROOT.gSystem.Load("libBambooLumiMask.so")
#ROOT.gROOT.ProcessLine('#include "../ext/include/LumiMask.h"')
lm = ROOT.LumiMask.fromJSON(jsonName, *runRange)
assert not lm.accept(283681, 10 ) ## excluded run
assert     lm.accept(276948, 10 ) ## included run
assert     lm.accept(276948, 1  ) ## first for a run
assert     lm.accept(276834, 720) ## last for that run
assert not lm.accept(277070, 310) ## excluded part
assert     lm.accept(277070, 311) ## later block
print("Tests pass!")
