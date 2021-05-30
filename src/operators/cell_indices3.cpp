#include "mandoline/operators/cell_indices.hpp"
namespace mandoline::operators {

/// int nearest_vertex(const CutCellMesh<3>& ccm, Eigen::Ref<const mtao::Vec3d>

mtao::VecXi get_cell_indices(const CutCellMesh<3>& ccm,
                             Eigen::Ref<const mtao::ColVecs3d> p) {
    return get_cell_indices(ccm, ccm.exterior_grid().cell_ownership_grid(), p);
}
mtao::VecXi get_cell_indices(
    const CutCellMesh<3>& ccm,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const mtao::ColVecs3d> P) {
    mtao::VecXi I(P.cols());

#if defined(MTAO_TBB_ENABLED)
    tbb::parallel_for<int>(0, I.size(), [&](const int j) {
#else
    for (int j = 0; j < I.size(); ++j) {
#endif
        I(j) = get_cell_index(ccm, cell_ownership_grid, P);
#if defined(MTAO_TBB_ENABLED)
    });
#else
    }
#endif
    return I;
}

mtao::VecXi get_cell_indices(
    const CutCellMesh<3>::ExteriorGridType& exterior_grid,
    Eigen::Ref<const mtao::ColVecs3d> P) {
    return get_cell_indices(exterior_grid, exterior_grid.cell_ownership_grid(),
                            P);
}

mtao::VecXi get_cell_indices(
    const AdaptiveGrid& exterior_grid,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const mtao::ColVecs3d> P) {
    mtao::VecXi I(P.cols());

#if defined(MTAO_TBB_ENABLED)
    tbb::parallel_for<int>(0, I.size(), [&](const int j) {
#else
    for (int j = 0; j < I.size(); ++j) {
#endif
        I(j) = get_cell_index(exterior_grid, cell_ownership_grid, P);
#if defined(MTAO_TBB_ENABLED)
    });
#else
    }
#endif
    return I;
}

int get_cell_index(
    const CutCellMesh<3>& ccm,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const mtao::Vec3d> p) {
    auto [c, q] = ccm.vertex_grid().coord(p);
    // check if its  in an adaptive grid cell, tehn we can just use that cell
    // if the grid returns -2 we pass that through
    if (int ret = get_cell_index(ccm.exterior_grid(), cell_ownership_grid, c);
        ret != -1) {
        if (ret == -2) {
            mtao::logging::warn() << "Point lies outside the grid";
        }
        return ret;
    } else {
        auto cell_indices = ccm.cut_cells_in_grid_cell(c);
        mtao::ColVecs3d V;
        const auto& CV = ccm.cached_vertices();
        const mtao::ColVecs3d* VP = CV ? &*CV : &V;
        if (!CV) {
            V = ccm.vertices();
        }
        for (auto&& ci : cell_indices) {
            auto&& cell = ccm.cells().at(ci);
            if (cell.contains(*VP, ccm.faces(), p, ccm.folded_faces())) {
                return ci;
            }
        }
        // if I haven't returned yet then either this ccm is bad
        // or i'm too close to an edge. lets not assume that for now
        // if (!quiet_failures) {
        {
            // spdlog::warn(
            //    "There are {} cells in grid cell ({}) but Mandoline "
            //    "failed to find point ({}) in it. Falling back to max "
            //    "solid angles",
            //    cell_indices.size(), fmt::join(c, ","), fmt::join(p, ","));
            int max_cell = -1;
            double max_solid_angle = std::numeric_limits<double>::min();
            for (auto&& ci : cell_indices) {
                const auto& cell = ccm.cells().at(ci);
                double solid_angle =
                    cell.solid_angle(*VP, ccm.faces(), p, ccm.folded_faces());
                bool contains =
                    cell.contains(*VP, ccm.faces(), p, ccm.folded_faces());
                if (solid_angle > max_solid_angle) {
                    max_cell = ci;
                    max_solid_angle = solid_angle;
                }
                // spdlog::info(
                //    "Cell not found debug: got a solid angle of {} "
                //    "from "
                //    "cell "
                //    "{} with contains = {}",
                //    solid_angle, ci, contains);
            }
            return max_cell;
        }
    }
    return -1;
}
int get_cell_index(
    const AdaptiveGrid& exterior_grid,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const mtao::Vec3d> p) {
    const auto& vg = exterior_grid.Base::vertex_grid();
    auto [coord, quot] = vg.coord(p);
    return get_cell_index(exterior_grid, cell_ownership_grid, coord);
}

int get_cell_index(
    const AdaptiveGrid& exterior_grid,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    const std::array<int, 3>& coord) {
    auto ec = mtao::eigen::stl2eigen(coord);
    //if (exterior_grid.Base::cell_grid().valid_index(coord)) {
    if (ec.minCoeff() < 0 || (ec.array() >= mtao::eigen::stl2eigen(exterior_grid.Base::cell_grid().shape()).array()).any()) {
        return cell_ownership_grid(coord);
    } else {
        return -2;
    }
}
}  // namespace mandoline::operators
