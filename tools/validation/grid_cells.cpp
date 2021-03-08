#define _gmp_const const
#include <fmt/format.h>
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#include <spdlog/spdlog.h>

#include "cutmesh_validation.hpp"
using namespace mandoline;
std::array<int, 2> grid_cells_fully_utilized_count(const CutCellMesh<3>& ccm) {
    auto V = ccm.vertices();

    Eigen::MatrixXd IV = V.transpose();
    // int is_utilized = 0;
    // int is_not_utilized = 0;
    int unused_grid_cells = 0;
    int invalid_parent_cells = 0;
    mtao::VecXd vols = ccm.cell_volumes();
    if (vols.minCoeff() < 0) {
        spdlog::error("Negative volumes found");
    }

    spdlog::info("Creating sets");
    std::map<std::array<int, 3>, std::set<int>> grid_cell_children;
    // how the hell did i mess up implementing this in a way that breaks gdb?
    // auto grid_cell_children =
    //    mtao::geometry::grid::GridDataD<std::set<int>, 3>(ccm.cell_shape());

    auto grid_cell_volumes =
        mtao::geometry::grid::GridData3d::Constant(0.0, ccm.cell_shape());

    // for each cell add its volume to the grid
    for (auto&& [cid, cell] : mtao::iterator::enumerate(ccm.cells())) {
        double v = vols(cid);
        if (grid_cell_volumes.valid_index(cell.grid_cell)) {
            grid_cell_volumes(cell.grid_cell) += v;
            grid_cell_children[cell.grid_cell].emplace(cid);
        } else {
            invalid_parent_cells++;
        }
    }

    double cell_vol = ccm.dx().prod();
    // for (auto&& [c, v] : grid_cell_volumes) {
    //    spdlog::info("{} {}", fmt::join(c, ","), v);
    //}
    for (auto&& [cell_index, cell] : ccm.exterior_grid().cells()) {
        auto corner = cell.corner();
        int width = cell.width();
        std::array<int, 3> ijk;
        auto& [i, j, k] = ijk;

        double per_cell_vol = vols(cell_index) / (width * width * width);
        for (i = corner[0]; i < corner[0] + width; ++i) {
            for (j = corner[1]; j < corner[1] + width; ++j) {
                for (k = corner[2]; k < corner[2] + width; ++k) {
                    if (grid_cell_volumes.valid_index(ijk)) {
                        if (grid_cell_children.contains(ijk)) {
                            spdlog::warn(
                                "Grid cell {} shares a exterior cell {} with "
                                "normal cut-cells {}",
                                fmt::join(ijk, ","), cell_index,
                                fmt::join(grid_cell_children[ijk], ","));
                        }
                    }
                }
            }
        }
    }
    for (auto&& [cell_index, cell] : ccm.exterior_grid().cells()) {
        auto corner = cell.corner();
        int width = cell.width();
        std::array<int, 3> ijk;
        auto& [i, j, k] = ijk;

        double per_cell_vol = vols(cell_index) / (width * width * width);
        for (i = corner[0]; i < corner[0] + width; ++i) {
            for (j = corner[1]; j < corner[1] + width; ++j) {
                for (k = corner[2]; k < corner[2] + width; ++k) {
                    if (grid_cell_volumes.valid_index(ijk)) {
                        grid_cell_volumes(ijk) += per_cell_vol;
                        grid_cell_children[ijk].emplace(cell_index);
                    } else {
                        invalid_parent_cells++;
                    }
                }
            }
        }
    }
    // auto& active_cells = ccm.active_grid_cell_mask();
    // int active_count =
    //    std::count(active_cells.begin(), active_cells.end(), true);
    // int unused_cells = active_count - grid_cell_volumes.size();
    // if (unused_cells > 0) {
    //    std::cout << "Missing " << unused_cells << " cells" << std::endl;
    //}

    int bad_vols = 0;
    grid_cell_volumes.loop([&](const std::array<int, 3>& c, const double& v) {
        if (std::abs(1 - v / cell_vol) > 1e-5) {
            if (v == 0) {
                unused_grid_cells++;
            } else {
                const auto& children = grid_cell_children[c];
                spdlog::warn(
                    "Grid cell {} had an invalid volume inside: {} / {} "
                    "(children were {})",
                    fmt::join(c, ","), v, cell_vol, fmt::join(children, ","));
                for (auto&& c : children) {
                    std::cout << fmt::format("cell[{},vol={}] ", c, vols(c));
                }
                std::cout << std::endl;
                bad_vols += 1;
            }
        }
    });
    if (invalid_parent_cells > 0) {
        spdlog::error("Invalid parent cells: {}", invalid_parent_cells);
    }
    return std::array<int, 2>{{bad_vols, unused_grid_cells}};
}
bool grid_cells_fully_utilized(const mandoline::CutCellMesh<3>& ccm) {
    auto [a, b] = grid_cells_fully_utilized_count(ccm);
    return a == 0 && b == 0;
}
