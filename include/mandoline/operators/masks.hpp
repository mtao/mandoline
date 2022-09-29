#pragma once
#include "mandoline/mesh3.hpp"
namespace mandoline::operators {
balsa::eigen::VecXd mesh_face_mask(const CutCellMesh<3> &ccm);//for removing mesh faces
balsa::eigen::VecXd grid_boundary_mask(const CutCellMesh<3> &ccm);//for removing mesh faces
}// namespace mandoline::operators
