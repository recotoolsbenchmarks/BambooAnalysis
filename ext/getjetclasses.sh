#!/bin/bash
CMSSW_VERSION="CMSSW_10_5_0_pre2"
CMSSW_CVMFS_BASE="/cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/${CMSSW_VERSION}"
CMSSW_GITHUB_BASE="https://raw.githubusercontent.com/cms-sw/cmssw/${CMSSW_VERSION}"

DEST_INC="CMSJet/include"
DEST_SRC="CMSJet/src"
mkdir -p "${DEST_INC}"
mkdir -p "${DEST_SRC}"
if [ -d "/cvmfs/cms.cern.ch/slc7_amd64_gcc700" ]
then
  function get_cms {
    cp "${CMSSW_CVMFS_BASE}/src/${1}" "${2}"
  }
else
  function get_cms {
    wget -P "${2}" "${CMSSW_GITHUB_BASE}/${1}" 2> /dev/null || exit 1
  }
fi

anyWasModified=""
function get_file {
  local fname="${1}"
  local dest="${2}"
  if [ ! -e "${dest}/$(basename ${fname})" ]; then
    get_cms "${fname}" "${dest}"
    anyWasModified="yes"
  fi
}

here=$(realpath $(dirname $0))

## copy headers and implementation files
while read filename
do
  if [[ "${filename}" == *"/interface/"* ]]; then
    get_file "${filename}" "${DEST_INC}"
  else
    get_file "${filename}" "${DEST_SRC}"
  fi
done < "${here}/jetclasses_filenames.txt"

## now patch them
if [ "${anyWasModified}" != "" ]; then
  pushd CMSJet > /dev/null
  patch -p2 -i "${here}/jetclasses.patch" > /dev/null
  popd > /dev/null
fi
