#!/usr/bin/env bash
if [[ -z "${VIRTUAL_ENV}" ]]; then
  echo "This script is only for installing lwtnn in a virtualenv, please activate or set $VIRTUAL_ENV first"
  return 1
fi
lwtnnversion="2.10"
lwtnndllink="https://github.com/lwtnn/lwtnn/archive/v${lwtnnversion}.tar.gz"
workdir=$(mktemp -d)
echo "---> Downloading in ${workdir}"
pushd "${workdir}" &>/dev/null
wget "${lwtnndllink}"
instdir="${VIRTUAL_ENV}"
echo "---> Installing in ${instdir}"
if [ ! -d "${instdir}" ]; then
  mkdir -p "${instdir}"
fi
extract "${workdir}/v${lwtnnversion}.tar.gz"
mkdir build
pushd build &>/dev/null
cmake -DCMAKE_INSTALL_PREFIX="${VIRTUAL_ENV}" -DCMAKE_BUILD_TYPE=Release -DBUILTIN_BOOST=true -DBUILTIN_EIGEN=ON "../lwtnn-${lwtnnversion}" || return 1
make || return 1
make install || return 1
popd &>/dev/null # out of build
popd &>/dev/null # back to start
if [ -n "${workdir}" ]; then
  rm -rf "${workdir}"
fi
