#pragma once
#include "mandoline/mesh3.hpp"
namespace mandoline::operators {
    mtao::VecXd mesh_face_mask(const CutCellMesh<3>& ccm);//for removing mesh faces
    mtao::VecXd grid_boundary_mask(const CutCellMesh<3>& ccm);//for removing mesh faces
}
