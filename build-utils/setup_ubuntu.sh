#/bin/bash

git submodule update --recursive --init

#https://stackoverflow.com/questions/59895/get-the-source-directory-of-a-bash-script-from-within-the-script-itself
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

MANDOLINE_DIR="${SCRIPT_DIR}/../"
IMGUI_DIR="${MANDOLINE_DIR}/extern/core/extern/imgui"

CUR_DIR="$(pwd)"
BUILD_DIR="${MANDOLINE_DIR}/build"

mkdir -p "${BUILD_DIR}"

bash ${SCRIPT_DIR}/fetch_dependencies_ubuntu.sh

bash ${SCRIPT_DIR}/magnum_setup.sh "${SCRIPT_DIR}" "${SCRIPT_DIR}"
sudo bash ${SCRIPT_DIR}/magnum_setup_ubuntu.sh "${SCRIPT_DIR}" "${IMGUI_DIR}" "${SCRIPT_DIR}" 

pushd ${MANDOLINE_DIR}/extern
hg clone http://bitbucket.org/eigen/eigen
popd

pushd "${BUILD_DIR}"
cmake .. -DCMAKE_BUILD_TYPE=Release -DMTAO_CUSTOM_EIGEN_PATH=${MANDOLINE_DIR}/extern/eigen
make -j$( nproc ) make_cutmesh_gui exploded_mesh 
popd
