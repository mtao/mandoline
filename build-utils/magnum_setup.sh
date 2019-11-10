#!/bin/bash

TMPDIR="$1"


mkdir -p ${TMPDIR}/mosra; pushd ${TMPDIR}/mosra;
for repo in corrade magnum magnum-integration;
do git clone https://github.com/mosra/$repo;
    pushd $repo
    git checkout tags/v2019.01
    popd
done

popd #mosra
