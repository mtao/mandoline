#!/bin/bash
realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

TMPDIR="$1"
IMGUI_DIR=$( realpath "$2" )
SCRIPT_DIR=$( realpath "$3" )
pushd "${TMPDIR}/mosra"
pushd magnum-integration;

sed "s|IMDIR|${IMGUI_DIR}|g" ${SCRIPT_DIR}/magnum_integration_deb.patch > magnum_integration_imgui.patch
patch -b ./package/debian/rules  magnum_integration_imgui.patch
popd #magnum integration

for repo in corrade magnum magnum-integration;
do pushd $repo;
    ln -s package/debian .
    # ignoring the dsc errors taht happen from magnum
    dpkg-buildpackage || true
    popd
    dpkg -i ${repo}*.deb
done
popd
