#include "cutmesh_validation.hpp"
using namespace mandoline;

using CoordType = std::array<int,3>;
auto possible_cells(const CutCellMesh<3>& ccm, const std::vector<int>& face) -> std::set<CoordType> {

    if(face.empty()) { return {}; }
    auto possible = [&](int idx) -> std::set<std::array<int,3>> {
        if(ccm.is_grid_vertex(idx)) {
            return Vertex<3>(ccm.vertex_grid().unindex(idx)).possible_cells();
        } else {
            return {};
        }
    };

    std::set<CoordType> possibles = possible(face[0]);

    for(auto&& f: face) {
        auto s = possible(f);
        std::set<CoordType> i;
        std::set_intersection(possibles.begin(),possibles.end(),s.begin(),s.end(),std::inserter(i,i.end()));
        possibles = std::move(i);
        if(possibles.empty()) {
            return {};
        }
    }

    return possibles;
}
