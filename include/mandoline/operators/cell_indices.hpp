#pragma once
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/geometry/grid/staggered_grid.hpp>

#include "mandoline/adaptive_grid.hpp"
#include "mandoline/mesh2.hpp"
#include "mandoline/mesh3.hpp"

namespace mandoline::operators {

/// int nearest_vertex(const CutCellMesh<3>& ccm, Eigen::Ref<const balsa::eigen::Vec3d>

balsa::eigen::VecXi get_cell_indices(const CutCellMesh<3>& ccm,
                             Eigen::Ref<const balsa::eigen::ColVecs3d> p);
balsa::eigen::VecXi get_cell_indices(
    const CutCellMesh<3>& ccm,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const balsa::eigen::ColVecs3d> p);

balsa::eigen::VecXi get_cell_indices(
    const CutCellMesh<3>::ExteriorGridType& exterior_grid,
    Eigen::Ref<const balsa::eigen::ColVecs3d> p);

balsa::eigen::VecXi get_cell_indices(
    const AdaptiveGrid& exterior_grid,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const balsa::eigen::ColVecs3d> p);

int get_cell_index(
    const CutCellMesh<3>& ccm,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const balsa::eigen::Vec3d> p);
int get_cell_index(
    const AdaptiveGrid& exterior_grid,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    Eigen::Ref<const balsa::eigen::Vec3d> p);

int get_cell_index(
    const AdaptiveGrid& exterior_grid,
    const mtao::geometry::grid::GridDataD<int, 3>& cell_ownership_grid,
    const std::array<int, 3>& coord);
}  // namespace mandoline::operators
