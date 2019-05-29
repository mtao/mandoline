#pragma once
#include "mandoline/mesh3.hpp"


namespace mandoline::construction {
    CutCellMesh<3> from_grid_unnormalized(const mtao::ColVecs3d& V, const std::array<int,3>& cell_shape, int adaptive_level = 0, std::optional<double> threshold = {});
    CutCellMesh<3> from_bbox(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const Eigen::AlignedBox<double,3>& bbox, const std::array<int,3>& cell_shape, int adaptive_level = 0, std::optional<double> threshold = {});
    CutCellMesh<3> from_grid(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F, const mtao::geometry::grid::StaggeredGrid3d& grid, int adaptive_level = 0, std::optional<double> threshold = {});
}
