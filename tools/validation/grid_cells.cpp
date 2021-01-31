#define __gmp_const const
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#include <spdlog/spdlog.h>

#include "cutmesh_validation.hpp"
using namespace mandoline;
std::array<int, 2> grid_cells_fully_utilized_count(const CutCellMesh<3>& ccm) {
    auto V = ccm.vertices();

    Eigen::MatrixXd IV = V.transpose();
    int is_utilized = 0;
    int is_not_utilized = 0;
    mtao::VecXd vols = ccm.cell_volumes();
    std::map<std::array<int, 3>, double> grid_cell_volumes;
    for (auto&& [cid, cell] : mtao::iterator::enumerate(ccm.cells())) {
        double v = vols(cid);
        auto [it, happened] = grid_cell_volumes.try_emplace(cell.grid_cell, v);
        if (!happened) {
            it->second += v;
        }
    }

    //for (auto&& [c, v] : grid_cell_volumes) {
    //    spdlog::info("{} {}", fmt::join(c, ","), v);
    //}
    auto& active_cells = ccm.active_grid_cell_mask();
    int active_count =
        std::count(active_cells.begin(), active_cells.end(), true);
    int unused_cells = active_count - grid_cell_volumes.size();
    if (unused_cells > 0) {
        std::cout << "Missing " << unused_cells << " cells" << std::endl;
    }

    int bad_vols = 0;
    double cell_vol = ccm.dx().prod();
    for (auto&& [c, v] : grid_cell_volumes) {
        if (std::abs(1 - v / cell_vol) > 1e-5) {
            bad_vols += 1;
        } else if (!active_cells.valid_index(c) || !active_cells(c)) {
            //spdlog::info("{} : {} {}", fmt::join(c, ","),
            //             active_cells.valid_index(c), active_cells(c));
            unused_cells++;
        }
    }

    return {{bad_vols, unused_cells}};
}
bool grid_cells_fully_utilized(const mandoline::CutCellMesh<3>& ccm) {
    auto [a, b] = grid_cells_fully_utilized_count(ccm);
    return a == 0 && b == 0;
}
