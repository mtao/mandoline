#pragma once
#include <mtao/types.hpp>
#include <Eigen/Sparse>
#include <set>
#include <vector>
#include <array>
namespace mandoline {

class DomainBoundary {
  public:
    // -1 = masked off, -2 == boundary of the domain itself

    // boundary facet array should be [positive,negative] to make 1D look like
    // [- cell A +] |FACE| [- cell B +]

    // ignores masked off boundaries and optionally domain boundaries
    std::vector<Eigen::Triplet<double>> boundary_triplets(bool use_domain_boundaries = false) const;
    // computes the boundary of the adaptive grid, including the interior boundaries
    mtao::VecXd boundary_face_mask() const;
    // computes facets on the interior
    mtao::VecXd interior_face_mask() const;
    // includes the boundary of hte adaptive grid, ignoring hte interior boundaries
    std::set<int> boundary_facets() const;// grid boundary facets using grid indices
    std::set<int> interior_facets() const;// grid boundary facets using grid indices
    const std::vector<std::array<int, 2>> &boundary_facet_pairs() const;
    const std::array<int, 2> &boundary_facet_pair(size_t idx) const;

    size_t boundary_facet_size() const;
    bool is_boundary_facet(int idx) const;
    static bool is_boundary_facet(const std::array<int, 2> &pr) { return pr[0] == -2 || pr[1] == -1; }

  protected:
    std::vector<std::array<int, 2>> m_boundary_facet_pairs;
};
}// namespace mandoline
