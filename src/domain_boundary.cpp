#include "mandoline/domain_boundary.hpp"
#include <mtao/iterator/enumerate.hpp>
#include "mandoline/operators/boundary.hpp"

namespace mandoline {


bool DomainBoundary::is_boundary_facet(int idx) const {
    return DomainBoundary::is_boundary_facet(m_boundary_facet_pairs.at(idx));
}

std::vector<Eigen::Triplet<double>> DomainBoundary::boundary_triplets(bool use_domain_boundaries) const {
    return operators::boundary_triplets(*this,use_domain_boundaries);
}
balsa::eigen::VecXd DomainBoundary::boundary_face_mask() const {
    balsa::eigen::VecXd R(boundary_facet_size());
    R.setZero();
    for (auto &&ind : boundary_facets()) {
        R(ind) = 1;
    }
    return R;
}
balsa::eigen::VecXd DomainBoundary::interior_face_mask() const {
    balsa::eigen::VecXd R(boundary_facet_size());
    R.setOnes();
    // we'll use boundary facets because there should be less of them
    for (auto &&ind : boundary_facets()) {
        R(ind) = 0;
    }
    return R;
}
std::set<int> DomainBoundary::boundary_facets() const {
    std::set<int> ret;
    for (auto &&[idx, pr] : boundary_facet_pairs()) {
        if (is_boundary_facet(pr)) {
            ret.insert(idx);
        }
    }
    return ret;
}
std::set<int> DomainBoundary::interior_facets() const {
    std::set<int> ret;
    for (auto &&[idx, pr] : boundary_facet_pairs()) {
        if (!is_boundary_facet(pr)) {
            ret.insert(idx);
        }
    }
    return ret;
}
const std::vector<std::array<int, 2>> &DomainBoundary::boundary_facet_pairs() const { return m_boundary_facet_pairs; }
size_t DomainBoundary::boundary_facet_size() const {
    return m_boundary_facet_pairs.size();
}
const std::array<int, 2> &DomainBoundary::boundary_facet_pair(size_t idx) const {
    return m_boundary_facet_pairs.at(idx);
}
}// namespace mandoline
