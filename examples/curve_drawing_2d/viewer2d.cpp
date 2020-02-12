#include "viewer2d.h"
#include "imgui.h"
#include <GL/gl.h>
#include <iterator>
#include <mtao/iterator/enumerate.hpp>

#include "mandoline/diffgeo_utils.hpp"

void CutcellViewer::initialize() {
}
void CutcellViewer::set_data(const CutCellMesh<2>& ccm) {
    has_data = true;

    /*
    m_ccm = ccm;
    auto&& V = ccm.vertices();
    auto&& hem = ccm.hem;
    //auto&& bE = ccm.boundary_edges;
    auto&& mask = ccm.active_cell_mask;



    vertices = V;
    dual_vertices = ccm.dual_vertices();

    cutC.resize(ccm.cell_size());
    edges.clear();
    cutregion_boundary.clear();
    dual_edges.clear();
    for(int i = 0; i < hem.size(); ++i) {
        auto e = hem.edge(i);
        edges.push_back(std::array<int,2>{{e.vertex(),e.get_dual().vertex()}});
        edgesCell.push_back(e.cell());
        edgesDual.push_back(e.dual_index());
        if(e.cell() >= 0 && e.cell() < cutC.size()) {
        cutC[e.cell()].push_back(i);
        }
        if(e.cell() < ccm.grid_size()) {
            edgesboundary.push_back(i);
        }
    }
    for(auto&& [idx,c]: mtao::iterator::enumerate(cutC)) {
        if(c.size() > 20) {
            draw_cell_idx = idx;
            std::cout << idx << ":";
            std::copy(c.begin(),c.end(),
                    std::ostream_iterator<int>(std::cout,","));
            std::cout << std::endl;
        }
    }
    //std::copy(bE.begin(),bE.end(),std::back_inserter(edges));
    //std::copy(bE.begin(),bE.end(),std::back_inserter(edgesRigidBoundary));

    {
        auto V = ccm.dual_vertices();

        cell_centers.clear();
        std::transform(V.begin(),V.end(),std::back_inserter(cell_centers),[](const Vec2d& v) -> Vec2d{ return Vec2d(v); });
    }
    boundary_op = ccm.boundary();
    pressure_solve();
    */
}

void CutcellViewer::draw() const {

}


void CutcellViewer::gui() {
}

