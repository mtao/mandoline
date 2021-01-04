#include "mandoline/mesh3.hpp"
#include <fmt/format.h>

std::map<int,int> exterior_grid_valences(const mandoline::CutCellMesh<3>& ccm) {

    std::map<int,int> valences;
    auto cell_grid = ccm.exterior_grid().cell_ownership_grid();
    for(auto&& f: ccm.cut_faces()) {
        if(f.external_boundary) {
            auto [exterior_cell, sgn] = *f.external_boundary;
            if(exterior_cell < 0) continue;
            int agindex = cell_grid.get(exterior_cell);
            //fmt::print("Stencil boundary ({1}({0}) {2})\n", exterior_cell, agindex, sgn);
            if(valences.find(agindex) != valences.end()) {
                valences[agindex]++;
            } else {
                valences[agindex] = 1;
            }
        }
    }
    std::cout << std::endl;
    for(auto&& face: ccm.exterior_grid().faces()) {
        auto&& [n,p] = face.dual_edge;
        //fmt::print("dual edge ({} {})\n", n, p);
        if(n >= 0) {
            
            if(valences.find(n) != valences.end()) {
                //valences[cell_grid.get(n)]++;
                valences[n]++;
            } else {
                //valences[cell_grid.get(n)] = 1;
                valences[n] = 1;
            }
        }
        if(p >= 0) {
            
            if(valences.find(p) != valences.end()) {
                valences[p]++;
            } else {
                valences[p] = 1;
            }
        }
    }
    return valences;

    

}
// exterior cells should have valence 6
bool exterior_cell_valence_counts(const mandoline::CutCellMesh<3>& ccm) {
    auto valences = exterior_grid_valences(ccm);
    bool ret = true;
    for(auto&& [v,c]: valences) {
        if(c != 6) {
            ret = false;
            fmt::print("Bad exterior cell valence ({} had valence {})\n", v,c);
        }
    }
    return ret;
}

