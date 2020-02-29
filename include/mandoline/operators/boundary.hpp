#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {

Eigen::SparseMatrix<double> boundary(const CutCellMesh<3> &ccm);
std::set<int> grid_boundary_faces(const CutCellMesh<3> &ccm);//for removing mesh faces
}// namespace mandoline::operators
