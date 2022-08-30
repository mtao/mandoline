#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {

    balsa::eigen::ColVecs3d cell_centroids(const CutCellMesh<3>& ccm);
    balsa::eigen::ColVecs3d face_centroids(const CutCellMesh<3>& ccm);

}// namespace mandoline::operators
