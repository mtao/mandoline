#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {

    mtao::ColVecs3d cell_centroids(const CutCellMesh<3>& ccm);
    mtao::ColVecs3d face_centroids(const CutCellMesh<3>& ccm);

}// namespace mandoline::operators
