#!/bin/bash
CMSSW_REL_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0_pre2"
JMOPKG="${CMSSW_REL_BASE}/src/CondFormats/JetMETObjects"
cp "${JMOPKG}/interface/Utilities.h" include/
jmoClassNames=(JetCorrectionUncertainty JetCorrectorParameters JetCorrectorParametersHelper JetResolutionObject SimpleJetCorrectionUncertainty SimpleJetCorrector FactorizedJetCorrector FactorizedJetCorrectorCalculator)
for className in "${jmoClassNames[@]}"
do
  cp "${JMOPKG}/interface/${className}.h" include/
  cp "${JMOPKG}/src/${className}.cc" src/
done
JMMPKG="${CMSSW_REL_BASE}/src/JetMETCorrections/Modules"
cp "${JMMPKG}/interface/JetResolution.h" include/
cp "${JMMPKG}/src/JetResolution.cc" src/
CMUPKG="${CMSSW_REL_BASE}/src/CommonTools/Utils"
cp "${CMUPKG}/interface/FormulaEvaluator.h" include/
cp "${CMUPKG}/src/FormulaEvaluator.cc" src/
for ffn in $(find "${CMUPKG}/src" -type f -name "formula*")
do
  cp "${ffn}" src/
done
patch -p2 -i jetclasses.patch
