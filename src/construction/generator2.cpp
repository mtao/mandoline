#include "mandoline/construction/generator2.hpp"
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/colvector_loop.hpp>

namespace mandoline::construction {
template<>
CutCellMesh<2> CutCellEdgeGenerator<2>::generate() const {
    return generate_edges();
}


CutCellMesh<2> CutCellGenerator<2>::generate() const {
    auto ccm = generate_faces();
    return ccm;
}


void CutCellGenerator<2>::add_boundary_elements(const BoundaryElements &E) {
    //BoundaryElements
    add_edges(E);
}


void CutCellGenerator<2>::bake_faces() {

    std::tie(hem, std::ignore) = compute_planar_hem(all_GV(), mtao::eigen::stl2eigen(edges()), m_active_grid_cell_mask, StaggeredGrid::cell_size());

    std::vector<bool> active(data().nE());// edges that sit on an axis
    auto &&ti = data().edge_intersections();
    std::transform(ti.begin(), ti.end(), active.begin(), [](auto &&ti) {
        return ti.active();
    });
    cut_cell_to_primal_map.clear();


    mtao::logging::debug() << "Number of faces: " << hem.cell_halfedges().size();

    for (auto &&[i, f] : mtao::iterator::enumerate(m_cut_faces)) {
        cut_cell_to_primal_map[i] = f.parent_fid;
        if (active[f.parent_fid]) {
            axial_primal_faces[i].insert(smallest_ordered_edge(f.indices));
        }
    }
}
template<>
CutCellMesh<2> CutCellEdgeGenerator<2>::generate_faces() const {
    // generate the edge structure
    auto ret = generate_edges();


    // make a halfedge mesh and store it
    std::tie(ret.hem, std::ignore) = compute_planar_hem(all_V(), mtao::eigen::stl2eigen(edges()), m_active_grid_cell_mask, 0);


    {
        //auto t = mtao::logging::timer("Gluing halfedges to edges");
        mtao::map<std::array<int, 2>, int> emap;
        for (int i = 0; i < ret.m_cut_edges.size(); ++i) {
            std::array<int, 2> v = ret.m_cut_edges[i].indices;
            emap[v] = i;
            std::swap(v[0], v[1]);
            emap[v] = i;
        }
        ret.halfedges_per_edge.resize(2, ret.m_cut_edges.size());


        for (int i = 0; i < ret.hem.size(); ++i) {
            std::array<int, 2> v;
            auto e = ret.hem.edge(i);
            v[0] = e.vertex();
            e.dual();
            int dhe = int(e);
            v[1] = e.vertex();
            int idx = ret.halfedge_to_edge_index[i] = emap[v];
            auto hpee = ret.halfedges_per_edge.col(idx);
            auto &&ee = ret.m_cut_edges[idx].indices;
            if (v[0] == ee[0]) {
                hpee(0) = i;
                hpee(1) = dhe;
            } else {
                hpee(1) = i;
                hpee(0) = dhe;
            }
        }
    }

    ret.m_active_grid_cell_mask = m_active_grid_cell_mask;
    ret.exterior_grid = ExteriorGrid<2>(m_active_grid_cell_mask);
    ret.active_cell_mask.resize(ret.cell_size());
    ret.active_cell_mask.setConstant(1);

    ret.active_cell_mask.topRows(m_active_grid_cell_mask.size()) = m_active_grid_cell_mask.as_eigen_vector().cast<double>();


    {

        auto VV = all_GV();
        auto vols = ret.hem.signed_areas(VV);
        auto mesh_faces = ret.hem.cells_multi_component_map();
        for (auto &&[i, v] : mesh_faces) {
            CutFace<2> F;
            F.indices = v;
            if (!F.indices.empty()) {
                double vol = 0;
                for (auto &&v : F.indices) {
                    vol += mtao::geometry::curve_volume(VV, v);
                }
                if (vol > 0) {
                    ret.m_faces.emplace_back(std::move(F));
                }
            }
        }
    }
    mtao::logging::debug() << "Making cut-faces";
    auto AV = all_vertices();

    for (auto &&[idx, f] : mtao::iterator::enumerate(ret.m_faces)) {
        auto c = f.possible_cells(AV);
        assert(c.size() == 1);
        int grid_cell = StaggeredGrid::cell_index(*c.begin());
        ret.cell_grid_ownership[grid_cell].insert(idx);
    }

    auto make_boundary_pair = [&](auto &&boundary_facet) -> std::optional<std::tuple<int, bool>> {
        auto pc = boundary_facet.possible_cells(AV);
        if (pc.size() != 2) {
            return {};
        }
        std::array<CoordType, 2> pca{ { { {} }, { {} } } };
        std::copy(pc.begin(), pc.end(), pca.begin());
        //find the boundary cells axis
        int idx;
        for (idx = 0; idx < 2; ++idx) {
            if (pca[0][idx] != pca[1][idx]) {
                break;
            }
        }
        if (pca[0][idx] + 1 != pca[1][idx]) {
            mtao::logging::error() << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?";
        }
        if (pca[0][idx] < 0) {
            return {{ -1, 1 }};
        } else if (pca[1][idx] >= cell_shape()[idx]) {
            return {{ -1, 0 }};
        } else {
            int pi = cell_index(pca[0]);
            int ni = cell_index(pca[1]);
            bool pa = m_active_grid_cell_mask.get(pi);
            bool na = m_active_grid_cell_mask.get(ni);
            bool sign = (idx%2==1);
            if (na ^ pa) {//active inactive boundary
                if (pa) {
                    return {{ pi, sign }};
                } else {
                    return {{ ni, !sign }};
                }
            }
        }
        return {};
    };
    mtao::logging::debug() << "Making cut-edge boundaries";
    std::map<Edge,std::pair<int,bool>> diredge_to_edge_orientation;
    for (auto &&[idx, e] : mtao::iterator::enumerate(ret.m_cut_edges)) {
        e.external_boundary = make_boundary_pair(e);

        auto c = e.indices;
        diredge_to_edge_orientation[c] = {idx,false};
        std::swap(c[0],c[1]);
        diredge_to_edge_orientation[c] = {idx,true};
    }
    for (auto &&[idx, f] : mtao::iterator::enumerate(ret.m_faces)) {
        for(auto&& c: f.indices) {
            auto& fm = ret.m_face_boundary_map[idx];
            for(int j = 0; j < c.size(); ++j) {
                int k = (j+1)%c.size();
                Edge me{{c[j],c[k]}};
                auto pr = diredge_to_edge_orientation.at(me);
                auto [eidx,sgn] = pr;
                auto e = ret.m_cut_edges[eidx].indices;
                fm.emplace(pr);
            }
        }
    }

    mtao::logging::debug() << "Original vertices: " << origV().size();
    //extra_metadata(ret);
    ret.m_origV.resize(2, origV().size());
    for (int i = 0; i < origV().size(); ++i) {
        ret.m_origV.col(i) = origV()[i];
    }
    ret.m_origE = data().E();
    return ret;
}


void CutCellGenerator<2>::extra_metadata(CutCellMesh<2> &mesh) const {
}


}// namespace mandoline::construction
