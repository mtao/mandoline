#include "mandoline/operators/boundary.hpp"
#include <mtao/iterator/enumerate.hpp>

namespace mandoline::operators {
std::set<int> grid_boundary_faces(const DomainBoundary& db, int offset) {
    std::set<int> ret;
    for(auto&& [idx,bfp]: mtao::iterator::enumerate(db.boundary_facet_pairs())) {
        if(db.is_boundary_facet(bfp)) {
            ret.emplace(idx + offset);
        }
    }
    return ret;
}
std::vector<Eigen::Triplet<double>> boundary_triplets(const DomainBoundary& db, bool include_domain_boundary) {
    std::vector<Eigen::Triplet<double>> ret;
    for(auto&& [idx,bfp]: mtao::iterator::enumerate(db.boundary_facet_pairs())) {
        if(include_domain_boundary || !db.is_boundary_facet(bfp)) {
            if(int a = bfp[0]; a >= 0) { ret.emplace_back(idx,a, -1); }
            if(int a = bfp[1]; a >= 0) { ret.emplace_back(idx,a, 1); }
        }
    }
    return ret;
}

}
