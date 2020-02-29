#pragma once
#include "mandoline/domain_boundary.hpp"
#include "mandoline/adaptive_grid.hpp"


namespace mandoline::operators {
    // outputs the indices of faces in the domain boundary which are on the boundary of the domain
    // as this is typically used to augment a collection of cut-cells, we allow for an optional offset
std::set<int> grid_boundary_faces(const DomainBoundary& db, int offset = 0);
std::vector<Eigen::Triplet<double>> boundary_triplets(const DomainBoundary& db, bool include_domain_boundary = false);
}

