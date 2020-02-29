#pragma once
#include "mandoline/mesh2.hpp"


namespace mandoline::operators {

Eigen::SparseMatrix<double> boundary(const CutCellMesh<2> &ccm, bool include_domain_boundary = false);
std::set<int> grid_boundary_faces(const CutCellMesh<2> &ccm);//for removing mesh faces
}
