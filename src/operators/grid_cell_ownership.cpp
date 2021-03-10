#include "mandoline/operators/grid_cell_ownership.hpp"

#include <array>
#include <set>
#include <vector>

namespace mandoline::operators {

std::vector<std::set<int>> grid_cell_ownership(const CutCellMesh<3>& ccm) {
    auto& cell_grid = ccm.StaggeredGrid::cell_grid();
    std::vector<std::set<int>> ret(cell_grid.size());

    for (auto&& [cid, cell] : mtao::iterator::enumerate(ccm.cells())) {
        ret[cell_grid.index(cell.grid_cell)].emplace(cid);
    }
    for (auto&& [cell_index, cell] : ccm.exterior_grid().cells()) {
        auto corner = cell.corner();
        int width = cell.width();
        std::array<int, 3> ijk;
        auto& [i, j, k] = ijk;

        for (i = corner[0]; i < corner[0] + width; ++i) {
            for (j = corner[1]; j < corner[1] + width; ++j) {
                for (k = corner[2]; k < corner[2] + width; ++k) {
                    ret[cell_grid.index(ijk)].emplace(cell_index);
                }
            }
        }
    }
    return ret;
}
}  // namespace mandoline::operators
