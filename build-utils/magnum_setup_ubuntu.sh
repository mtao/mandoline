#!/bin/bash

pushd "$1"

for repo in corrade magnum magnum-integration;
do pushd $repo;
    ln -s package/debian .
    # ignoring the dsc errors taht happen from magnum
    dpkg-buildpackage || true
    popd
    dpkg -i ${repo}*.deb
done
popd
