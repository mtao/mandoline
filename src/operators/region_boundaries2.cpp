#include "mandoline/operators/region_boundaries2.hpp"

#include "mandoline/operators/boundary2.hpp"
namespace mandoline::operators {
std::set<int> region_boundaries(const CutCellMesh<2>& ccm) {
    std::set<int> ret;

    for (auto&& [eidx, cutedge] : mtao::iterator::enumerate(ccm.cut_edges())) {
        if (cutedge.is_mesh_edge() ||
            (cutedge.external_boundary &&
             std::get<0>(*cutedge.external_boundary) < 0)) {
            ret.emplace(eidx);
        }
    }
    int fidx_offset = ccm.cut_edge_size();
    for (auto&& [bidx, pr] :
         mtao::iterator::enumerate(ccm.exterior_grid.boundary_facet_pairs())) {
        if (mandoline::DomainBoundary::is_boundary_facet(pr)) {
            ret.emplace(fidx_offset + bidx);
        }
    }
    return ret;
}
std::map<int, bool> boundary_from_cells(const CutCellMesh<2>& ccm,
                                        const std::set<int>& cells) {
    auto B = boundary(ccm, true);
    balsa::eigen::VecXd C(ccm.num_cells());
    C.setZero();
    for (auto&& c : cells) {
        C(c) = 1;
    }
    balsa::eigen::VecXd I = B * C;
    std::map<int, bool> boundary;
    for (auto&& [idx, v] : mtao::iterator::enumerate(I)) {
        if (v != 0) {
            if (v == 1) {
                boundary[idx] = false;
            } else {
                boundary[idx] = true;
            }
        }
    }
    return boundary;
}

std::map<int, bool> region_boundaries(const CutCellMesh<2>& ccm,
                                      const std::set<int>& regions) {
    std::set<int> cells;
    for (auto&& [idx, region] : mtao::iterator::enumerate(ccm.regions())) {
        if (regions.contains(region)) {
            cells.emplace(idx);
        }
    }
    return boundary_from_cells(ccm, cells);
}
std::map<int, bool> region_boundaries(const CutCellMesh<2>& ccm, int region) {
    return region_boundaries(ccm, std::set<int>({region}));
}
}  // namespace mandoline::operators
