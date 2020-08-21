#!/usr/bin/env bash
if [[ -z "${VIRTUAL_ENV}" ]]; then
  echo "This script is only for installing tensorflow-c binaries in a virtualenv, please activate or set $VIRTUAL_ENV first"
  return 1
fi
tfcversion="1.14.0"
tfdllink="https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-${tfcversion}.tar.gz"
workdir=$(mktemp -d)
echo "---> Downloading in ${workdir}"
pushd "${workdir}" &>/dev/null
wget "${tfdllink}"
popd &>/dev/null
instdir="${VIRTUAL_ENV}"
echo "---> Installing in ${instdir}"
if [ ! -d "${instdir}" ]; then
  mkdir -p "${instdir}"
fi
pushd "${instdir}" &>/dev/null
tar xzf "${workdir}/$(basename "${tfdllink}")"
popd &>/dev/null
if [ -n "${workdir}" ]; then
  rm -rf "${workdir}"
fi
