#!/bin/bash
CMSSW_REL_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0_pre2"
JMOPKG="${CMSSW_REL_BASE}/src/CondFormats/JetMETObjects"
DEST_INC="CMSJet/include"
DEST_SRC="CMSJet/src"
mkdir -p "${DEST_INC}"
mkdir -p "${DEST_SRC}"

cp "${JMOPKG}/interface/Utilities.h" "${DEST_INC}"
jmoClassNames=(JetCorrectionUncertainty JetCorrectorParameters JetCorrectorParametersHelper JetResolutionObject SimpleJetCorrectionUncertainty SimpleJetCorrector FactorizedJetCorrector FactorizedJetCorrectorCalculator)
for className in "${jmoClassNames[@]}"
do
  cp "${JMOPKG}/interface/${className}.h" "${DEST_INC}"
  cp "${JMOPKG}/src/${className}.cc" "${DEST_SRC}"
done
JMMPKG="${CMSSW_REL_BASE}/src/JetMETCorrections/Modules"
cp "${JMMPKG}/interface/JetResolution.h" "${DEST_INC}"
cp "${JMMPKG}/src/JetResolution.cc" "${DEST_SRC}"
CMUPKG="${CMSSW_REL_BASE}/src/CommonTools/Utils"
cp "${CMUPKG}/interface/FormulaEvaluator.h" "${DEST_INC}"
cp "${CMUPKG}/src/FormulaEvaluator.cc" "${DEST_SRC}"
for ffn in $(find "${CMUPKG}/src" -type f -name "formula*")
do
  cp "${ffn}" "${DEST_SRC}"
done
pushd CMSJet
patch -p2 -i "$(dirname $0)/jetclasses.patch"
popd
