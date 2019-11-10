#!/bin/bash

IMGUI_DIR="$1"
SCRIPT_DIR="$2"
TMPDIR_DIR="$3"


mkdir -p ${TMPDIR}/mosra; pushd ${TMPDIR}/mosra;
for repo in corrade magnum magnum-integration;
do git clone https://github.com/mosra/$repo; done

pushd magnum-integration;
pwd
sed "s|IMDIR|${IMGUI_DIR}|g" ${SCRIPT_DIR}/magnum_integration.patch > magnum_integration_imgui.patch
patch -b ./package/debian/rules  magnum_integration_imgui.patch
popd #magnum integration
popd #mosra
