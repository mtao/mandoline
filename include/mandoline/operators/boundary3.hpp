#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::operators {

std::set<int> grid_boundary_faces(const AdaptiveGrid& db, int offset = 0);
    // standard finite volume boundary operator
Eigen::SparseMatrix<double> boundary(const CutCellMesh<3> &ccm, bool include_domain_boundary = false);
std::vector<Eigen::Triplet<double>> boundary_triplets(const AdaptiveGrid& ag, int offset = 0, bool domain_boundary = false);
// select the cut-faces that happen to lie along axis-aligned grid cells
std::set<int> grid_boundary_faces(const CutCellMesh<3> &ccm);//for removing mesh faces
}// namespace mandoline::operators
