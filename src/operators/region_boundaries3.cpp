#include "mandoline/operators/region_boundaries3.hpp"

#include "mandoline/operators/boundary3.hpp"
namespace mandoline::operators {
std::set<int> region_boundaries(const CutCellMesh<3>& ccm) {
    std::set<int> ret;

    for (auto&& [eidx, cutface] : mtao::iterator::enumerate(ccm.cut_faces())) {
        if (cutface.is_mesh_face() ||
            (cutface.external_boundary &&
             std::get<0>(*cutface.external_boundary) < 0)) {
            ret.emplace(eidx);
        }
    }
    int fidx_offset = ccm.cut_face_size();
    for (auto&& [i, face] :
         mtao::iterator::enumerate(ccm.exterior_grid().faces())) {
        const auto& e = face.dual_edge;
        // if we are not the domain boundary
        if (ccm.exterior_grid().is_boundary_face(e)) {
            ret.emplace(i + fidx_offset);
        }
    }
    return ret;
}
std::map<int, bool> boundary_from_cells(const CutCellMesh<3>& ccm,
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

std::map<int, bool> region_boundaries(const CutCellMesh<3>& ccm,
                                      const std::set<int>& regions) {
    std::set<int> cells;
    for (auto&& [idx, region] : mtao::iterator::enumerate(ccm.regions())) {
        if (regions.contains(region)) {
            cells.emplace(idx);
        }
    }
    return boundary_from_cells(ccm, cells);
}
std::map<int, bool> region_boundaries(const CutCellMesh<3>& ccm, int region) {
    return region_boundaries(ccm, std::set<int>({region}));
}
}  // namespace mandoline::operators
