#!/bin/bash
CMSSW_REL_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_5_0_pre2"

DEST_INC="CMSJet/include"
DEST_SRC="CMSJet/src"
mkdir -p "${DEST_INC}"
mkdir -p "${DEST_SRC}"
if [ -d "/cvmfs/cms.cern.ch" ]
then
  function get_cms {
    cp "${CMSSW_REL_BASE}/src/${1}" "${2}"
  }
else
  echo "Should download from github now..."
  exit 0
fi
anyWasModified=""

function getheader {
  local pkg="${1}"
  local cls="${2}"
  if [ ! -e "${DEST_INC}/${cls}.h" ]; then
    get_cms "${pkg}/interface/${cls}.h" "${DEST_INC}"
    anyWasModified="yes"
  fi
}
function getheadernp {
  local pkg="${1}"
  local cls="${2}"
  if [ ! -e "${DEST_INC}/${cls}.h" ]; then
    get_cms "${pkg}/src/${cls}.h" "${DEST_INC}"
    anyWasModified="yes"
  fi
}
function getimplem {
  local pkg="${1}"
  local cls="${2}"
  if [ ! -e "${DEST_SRC}/${cls}.cc" ]; then
    get_cms "${pkg}/src/${cls}.cc" "${DEST_SRC}"
    anyWasModified="yes"
  fi
}
function getclass {
  local pkg="${1}"
  local cls="${2}"
  getheader "${pkg}" "${cls}"
  getimplem "${pkg}" "${cls}"
}
function getclassnp {
  local pkg="${1}"
  local cls="${2}"
  getheadernp "${pkg}" "${cls}"
  getimplem "${pkg}" "${cls}"
}

## JEC/JER-related classes
getheader "CondFormats/JetMETObjects" "Utilities"
jmoClassNames=(JetCorrectionUncertainty JetCorrectorParameters JetCorrectorParametersHelper JetResolutionObject SimpleJetCorrectionUncertainty SimpleJetCorrector FactorizedJetCorrector FactorizedJetCorrectorCalculator)
for className in "${jmoClassNames[@]}"
do
  getclass "CondFormats/JetMETObjects" "${className}"
done
getclass "JetMETCorrections/Modules" "JetResolution"
## formula evaluator
getclass "CommonTools/Utils" "FormulaEvaluator"
getheadernp "CommonTools/Utils" "formulaBinaryOperatorEvaluator"
getheadernp "CommonTools/Utils" "formulaFunctionOneArgEvaluator"
getheadernp "CommonTools/Utils" "formulaFunctionTwoArgsEvaluator"
getheadernp "CommonTools/Utils" "formulaUnaryMinusEvaluator"
getclassnp "CommonTools/Utils" "formulaConstantEvaluator"
getclassnp "CommonTools/Utils" "formulaEvaluatorBase"
getclassnp "CommonTools/Utils" "formulaParameterEvaluator"
getclassnp "CommonTools/Utils" "formulaVariableEvaluator"

## now patch them
if [ "${anyWasModified}" != "" ]; then
  pushd CMSJet > /dev/null
  patch -p2 -i "$(dirname $0)/jetclasses.patch" > /dev/null
  popd > /dev/null
fi
