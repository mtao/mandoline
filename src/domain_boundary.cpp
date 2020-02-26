#include "mandoline/domain_boundary.hpp"
#include <mtao/iterator/enumerate.hpp>

namespace mandoline {


bool DomainBoundary::is_boundary_facet(int idx) const {
    return DomainBoundary::is_boundary_facet(m_boundary_facet_pairs.at(idx));
}

std::vector<Eigen::Triplet<double>> DomainBoundary::boundary_triplets(bool use_domain_boundaries) const {
    std::vector<Eigen::Triplet<double>> trips;
    auto include = [use_domain_boundaries](int idx) -> bool {

        if (use_domain_boundaries) {
            return idx >= 0 || idx == -2;
        } else {
            return (idx >= 0);
        }
    };
    for (auto &&[idx, pr] : mtao::iterator::enumerate(boundary_facet_pairs())) {
        auto [p, n] = pr;
        assert(n != -1 && p != -1);
        if (include(p) && n >= 0) {
            assert(n >= 0);
            trips.emplace_back(idx,n, 1);
        }
        if (include(n) && p >= 0) {
            assert(p >= 0);
            trips.emplace_back(idx,p, -1);
        }
    }
    return trips;
}
mtao::VecXd DomainBoundary::boundary_face_mask() const {
    mtao::VecXd R(boundary_facet_size());
    R.setZero();
    for (auto &&ind : boundary_facets()) {
        R(ind) = 1;
    }
    return R;
}
mtao::VecXd DomainBoundary::interior_face_mask() const {
    mtao::VecXd R(boundary_facet_size());
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
