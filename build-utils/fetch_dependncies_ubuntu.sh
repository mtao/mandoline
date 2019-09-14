#!/bin/bash

sudo apt install -y git cmake libboost-thread-dev libmpfr-dev libmpfrc++-dev libcgal-dev libeigen3-dev protobuf-compiler mercurial

git submodule update --recursive --init


mkdir mosra; pushd mosra;
for repo in corrade magnum magnum-integration;
do git clone https://github.com/mosra/$repo; done

pushd magnum-integration;
pwd
patch -b ./package/debian/rules  ../../build-utils/magnum_integration.patch 
popd #magnum integration

for repo in corrade magnum magnum-integration;
do pushd $repo;
    ln -s package/debian .
    dpkg-buildpackage
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

