#!/bin/bash

sudo apt install -y git cmake libboost-thread-dev libmpfr-dev libmpfrc++-dev libcgal-dev libeigen3-dev protobuf-compiler mercurial

git submodule update --recursive --init

#https://stackoverflow.com/questions/59895/get-the-source-directory-of-a-bash-script-from-within-the-script-itself
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

MANDOLINE_DIR="${SCRIPT_DIR}/../"
IMGUI_DIR="${MANDOLINE_DIR}/extern/core/extern/imgui"

CUR_DIR="$(pwd)"


mkdir mosra; pushd mosra;
for repo in corrade magnum magnum-integration;
do git clone https://github.com/mosra/$repo; done

pushd magnum-integration;
pwd
sed "s|IMDIR|${IMGUI_DIR}|g" ${SCRIPT_DIR}/magnum_integration.patch > magnum_integration_imgui.patch
patch -b ./package/debian/rules  magnum_integration_imgui.patch
popd #magnum integration

for repo in corrade magnum magnum-integration;
do pushd $repo;
    ln -s package/debian .
    # ignoring the dsc errors taht happen from magnum
    dpkg-buildpackage || true
    popd
done

sudo dpkg -i *.deb

popd #mosra

pushd extern
hg clone http://bitbucket.org/eigen/eigen
popd

mkdir build
pushd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DMTAO_CUSTOM_EIGEN_PATH=$(pwd)/../extern/eigen
make -j$( nproc ) make_cutmesh_gui exploded_mesh 
popd

