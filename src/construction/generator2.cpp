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
    extra_metadata(ccm);
    return ccm;
}


void CutCellGenerator<2>::add_boundary_elements(const BoundaryElements &E) {
    //BoundaryElements
    add_edges(E);
}


void CutCellGenerator<2>::bake_faces() {

    std::tie(hem, std::ignore) = compute_planar_hem(all_GV(), mtao::eigen::stl2eigen(edges()), m_active_grid_cell_mask);

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


    auto eE = mtao::eigen::stl2eigen(edges());
    if(eE.size() == 0) {
        return ret;
    }
#if defined(_DEBUG)
    assert(eE.minCoeff() >= 0);
    assert(eE.maxCoeff() < total_vertex_size());
#endif

    // make a halfedge mesh and store it
    std::tie(ret.hem, std::ignore) = compute_planar_hem(all_V(), eE , m_active_grid_cell_mask);


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
    ret.exterior_grid = ExteriorGrid<2>(*this,m_active_grid_cell_mask);
    ret.active_cell_mask.resize(ret.cell_size());
    ret.active_cell_mask.setConstant(1);

    ret.active_cell_mask.topRows(m_active_grid_cell_mask.size()) = m_active_grid_cell_mask.as_eigen_vector().cast<double>();


    mtao::logging::debug() << "Making cut-faces";
    auto AV = all_vertices();
    {
        auto is_boundary_facet = [&](const CutFace<2>& f) {
            
            auto c = f.possible_cells([&](int idx) { return AV.at(idx).coord; } );
            if(c.empty() ) {
                return true;
            } else if(c.size() == 1) {
                coord_type a = *c.begin();
                assert(m_active_grid_cell_mask.valid_index(a));

                return m_active_grid_cell_mask(a);
            }
            return false;
        };

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
                if (vol > 0 && !is_boundary_facet(F)) {
                    ret.m_faces.emplace_back(std::move(F));
                }
            }
        }
    }

    mtao::logging::debug() << "Filling cell-grid ownership";
    for (auto &&[idx, f] : mtao::iterator::enumerate(ret.m_faces)) {
        auto c = f.possible_cells([&](int idx) { return AV.at(idx).coord; } );
        assert(c.size() == 1);
        int grid_cell = StaggeredGrid::cell_index(*c.begin());
        ret.cell_grid_ownership[grid_cell].insert(idx);
    }

    auto make_boundary_pair = [&](auto &&boundary_facet) -> std::optional<std::tuple<int, bool>> {
        //auto pc = boundary_facet.possible_cells(AV);
        auto pc = boundary_facet.possible_cells([&](int idx) { return AV.at(idx).coord; } );
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
            return {{ -2, 1 }};
        } else if (pca[1][idx] >= cell_shape()[idx]) {
            return {{ -2, 0 }};
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
    ret.m_origV.resize(2, origV().size());
    for (int i = 0; i < origV().size(); ++i) {
        ret.m_origV.col(i) = origV()[i];
    }
    ret.m_origE = data().E();
    return ret;
}


void CutCellGenerator<2>::extra_metadata(CutCellMesh<2> &mesh) const {
    mtao::data_structures::DisjointSet<int> cell_ds;
    int num_cells = mesh.num_cells();
    int num_cutfaces = mesh.num_cutfaces();
    {
        auto t = mtao::logging::profiler("region disjoint set construction", false, "profiler");


        // exterior grid pairs
        for (int i = 0; i < num_cells + mesh.num_cutedges(); ++i) {
            cell_ds.add_node(i);
        }
        for (auto &&[a, b] : mesh.exterior_grid.boundary_facet_pairs()) {
            if (a >= 0 && b >= 0) {
                cell_ds.join(a +num_cutfaces ,b +num_cutfaces);
            }
        }

        for (auto [cid, edge_pairs] : mesh.face_boundary_map()) {
            for (auto &&[eid, s] : edge_pairs) {
                if (eid >= 0) {
                    auto &&e = mesh.cut_edges()[eid];
                    if (!e.is_mesh_edge()) {
                        cell_ds.join(cid, eid + num_cells);
                    }
                }

            }
        }
        // join boundary faces with teh exterior faces
        for(auto&& [eidx,e]: mtao::iterator::enumerate(mesh.cut_edges())) {

            if(e.external_boundary) {
                auto [of,sgn] = *e.external_boundary;
                if(of >= 0) {
                    cell_ds.join(eidx+num_cells,mesh.exterior_grid.cell_indices().get(of)+num_cutfaces);
                }
            }
        }

        cell_ds.reduce_all();
    }

    int outside_root = -1;
    {// find a grid-boundary cell. this is a substantial amount of computation to make the "outside" the "first region"
        int boundary_idx = -1;// a boundary face
        for (auto &&[a, b] : mesh.exterior_grid.boundary_facet_pairs()) {
            if(a == -2) {
                boundary_idx = b + num_cutfaces;
            } else if(b == -2) {
                boundary_idx = a + num_cutfaces;
            }
        }
        if(boundary_idx == -1) {// try to find a boundary edge, and then a face from it
            // search through cut-edges for something on the boundary
            for(auto&& [eidx,e]: mtao::iterator::enumerate(mesh.cut_edges())) {

                if(e.external_boundary) {
                    auto [of,sgn] = *e.external_boundary;
                    if(of == -2) {
                        boundary_idx = eidx;
                    }
                }
            }
            for (auto [cid, edge_pairs] : mtao::iterator::enumerate(mesh.exterior_grid.boundary_facet_pairs())) {
                auto &&[a,b] = edge_pairs;
                if(a == -2 && b >= 0) {
                    outside_root = b;
                    break;
                }
                if(b == -2 && a >= 0) {
                    outside_root = a;
                    break;
                }
            }
        } else {
            outside_root = cell_ds.get_root(boundary_idx).data;
        }
    }
    std::map<int, int> reindexer;
    if(outside_root >= 0) {
        reindexer[outside_root] = 0;
    }
    for (int i = 0; i < num_cells; ++i) {
        int root = cell_ds.get_root(i).data;
        if (root != outside_root && reindexer.find(root) == reindexer.end()) {
            reindexer[root] = reindexer.size();
        }
    }


    auto& eg_regions = mesh.exterior_grid.m_regions;
    eg_regions.resize(mesh.exterior_grid.num_cells());
    for (int i = 0; i < mesh.exterior_grid.num_cells(); ++i) {
        eg_regions[i] = reindexer[cell_ds.get_root(mesh.num_cutfaces() + i).data];
    }

    for (int i = 0; i < num_cutfaces; ++i) {
        mesh.m_faces[i].region = reindexer[cell_ds.get_root(i).data];
    }
}


}// namespace mandoline::construction
