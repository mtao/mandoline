#include "mandoline/construction/generator2.hpp"

#include <tbb/parallel_for.h>

#include <mtao/colvector_loop.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/geometry/mesh/edge_tangents.hpp>

namespace mandoline::construction {
template <>
CutCellMesh<2> CutCellEdgeGenerator<2>::generate() const {
    return generate_edges();
}

CutCellMesh<2> CutCellGenerator<2>::generate() const {
    auto ccm = generate_faces();
    extra_metadata(ccm);
    return ccm;
}

void CutCellGenerator<2>::add_boundary_elements(const BoundaryElements &E) {
    // BoundaryElements
    add_edges(E);
}

void CutCellGenerator<2>::bake_faces() {
    spdlog::warn("Making edges");
    std::vector<Edge> edges;
    edges.reserve(m_cut_edges.size() + m_grid_edges.size());
    std::transform(m_cut_edges.begin(), m_cut_edges.end(),
                   std::back_inserter(edges),
                   [](auto &&e) { return e.indices; });
    std::transform(m_grid_edges.begin(), m_grid_edges.end(),
                   std::back_inserter(edges),
                   [](auto &&e) { return e.indices; });

    spdlog::warn("Making tangents");
    int nE = data().nE();
    mtao::ColVecs2d T(2, nE + 2);
    T.topLeftCorner(2, nE) = mtao::geometry::mesh::edge_tangents(
        mtao::eigen::stl2eigen(origV()), data().E());
    // std::vector<bool> origESigns(nE);
    auto get_t = [&](size_t parent_idx, size_t idx) {
        if (idx < grid_vertex_size()) {
            return data().get_edge_coord(parent_idx, grid_vertex(idx));
        } else {
            return data().get_edge_coord(parent_idx, crossing(idx));
        }
    };
    /*
    tbb::parallel_for(tbb::blocked_range<size_t>(0,nE),[&](const
    tbb::blocked_range<size_t>& range) { for(size_t idx = range.begin(); idx !=
    range.end(); ++idx) { auto e = data().E().col(idx); origESigns[idx] = e(1) >
    e(0);

            }
    });
    */
    T.col(nE) << 1, 0;
    T.col(nE + 1) << 0, 1;
    spdlog::warn("Making tangent edge map (cutedges)");
    std::map<std::array<int, 2>, std::tuple<int, bool>> edge_map;
    std::transform(
        m_cut_edges.begin(), m_cut_edges.end(),
        std::inserter(edge_map, edge_map.end()),
        [&](auto &&e) -> std::pair<std::array<int, 2>, std::tuple<int, bool>> {
            std::array<int, 2> ei = e.indices;
            bool flip_sign = ei[0] > ei[1];
            if (flip_sign) {
                std::swap(ei[0], ei[1]);
            }
            double a = get_t(e.parent_eid, ei[0]);
            double b = get_t(e.parent_eid, ei[1]);
            bool rel_flip = b < a;
            return {ei, std::make_tuple(e.parent_eid, flip_sign ^ rel_flip)};
            ;
        });
    spdlog::warn("Making tangent edge map (grid)");
    std::transform(
        m_grid_edges.begin(), m_grid_edges.end(),
        std::inserter(edge_map, edge_map.end()),
        [&](auto &&e) -> std::pair<std::array<int, 2>, std::tuple<int, bool>> {
            std::array<int, 2> ei = e.indices;
            bool flip_sign = ei[0] > ei[1];
            if (flip_sign) {
                std::swap(ei[0], ei[1]);
            }
            int ubx = e.unbound_axis();
            auto a = GV(ei[0]);
            auto b = GV(ei[1]);
            bool sgn =
                flip_sign ^ (a.coord[ubx] > b.coord[ubx]) ||
                (a.coord[ubx] == b.coord[ubx] && a.quot(ubx) > b.quot(ubx));
            std::tuple<int, bool> tp{nE + ubx, sgn};
            return {ei, tp};
        });

#if defined(USE_TANGENT_PLANAR_HEM_2)
    spdlog::warn("Tangent-based planar hem compute");
    std::tie(m_hem, std::ignore) = compute_planar_hem(
        edge_map, T, mtao::eigen::stl2eigen(edges), m_active_grid_cell_mask);
#else
    // Eigen::MatrixXi hem_edges = m_hem.edges();
    // spdlog::warn("Done");
    std::tie(m_hem, std::ignore) = compute_planar_hem(
        all_GV(), mtao::eigen::stl2eigen(edges), m_active_grid_cell_mask);
    // if(hem_edges != m_hem.edges()) {
    //    spdlog::warn("FAILED to make equal HEM!");
    //}
#endif
    {
        auto VV = all_GV();

        FaceCollapser fc(mtao::eigen::stl2eigen(edges));
        fc.bake(VV, false, edge_map, T);
        auto f = fc.faces();
        std::cout << "facecollapser size: " << f.size() << std::endl;
    }

    std::vector<bool> active(data().nE());  // edges that sit on an axis
    auto &&ti = data().edge_intersections();
    std::transform(ti.begin(), ti.end(), active.begin(),
                   [](auto &&ti) { return ti.active(); });
    cut_cell_to_primal_map.clear();

    mtao::logging::debug() << "Number of faces: "
                           << m_hem.cell_halfedges().size();

    for (auto &&[i, f] : mtao::iterator::enumerate(m_cut_faces)) {
        cut_cell_to_primal_map[i] = f.parent_fid;
        if (active[f.parent_fid]) {
            axial_primal_faces[i].insert(smallest_ordered_edge(f.indices));
        }
    }
}
template <>
CutCellMesh<2> CutCellEdgeGenerator<2>::generate_faces() const {
    // generate the edge structure
    auto ret = generate_edges();

    /*

        auto eE = mtao::eigen::stl2eigen(edges());
        if (eE.size() == 0) {
            return ret;
        }
#if defined(_DEBUG)
        assert(eE.minCoeff() >= 0);
        assert(eE.maxCoeff() < total_vertex_size());
#endif

        // make a halfedge mesh and store it
        std::tie(ret.hem, std::ignore) = compute_planar_hem(all_V(), eE,
m_active_grid_cell_mask);


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
        */
    mtao::logging::debug() << "Making exterior grid";

    ret.m_active_grid_cell_mask = m_active_grid_cell_mask;
    ret.exterior_grid = ExteriorGrid<2>(*this, m_active_grid_cell_mask);

    mtao::logging::debug() << "Making cut-faces";
    auto AV = all_vertices();
    {
        auto is_boundary_facet = [&](const CutFace<2> &f) {
            auto c =
                f.possible_cells([&](int idx) { return AV.at(idx).coord; });
            // each cell should only lie in one cell. if this doesn't belong in
            // a cell we're on the boundary
            if (c.empty()) {
                return true;
            } else if (c.size() == 1) {
                coord_type a = *c.begin();
                // if the cell is a valid cell
                if (m_active_grid_cell_mask.valid_index(a)) {
                    return m_active_grid_cell_mask(a);
                } else {
                    // toss the external grid cells
                    return true;
                }
            }
            return true;
        };

        auto VV = all_GV();
        auto vols = m_hem.signed_areas(VV);
        auto mesh_faces = m_hem.cells_multi_component_map();
        for (auto &&[i, v] : mesh_faces) {
            CutFace<2> F;
            F.indices = v;
            if (!F.indices.empty()) {
                double vol = 0;
                for (auto &&v : F.indices) {
                    vol += mtao::geometry::curve_volume(VV, v);
                }
                // only count things that are _fairly_ negative as negative
                // (this isn't very robust)
                const bool has_neg_vol = vol < -.5;
                const bool is_boundary = is_boundary_facet(F);
                if (!has_neg_vol && !is_boundary) {
                    ret.m_faces.emplace_back(std::move(F));
                } else {
                    int gvs = grid_vertex_size();
                    if (has_neg_vol) {
                        std::cout << "Neg area: " << std::string(F) << ": "
                                  << vol << std::endl;
                        for (auto &&v : F.indices) {
                            for (auto &&i : v) {
                                if (i < gvs) {
                                    std::cout << "GV(";
                                    auto v = vertex_unindex(i);
                                    std::copy(v.begin(), v.end(),
                                              std::ostream_iterator<int>(
                                                  std::cout, ","));
                                    std::cout << ")";
                                    std::cout << " => ";
                                } else {
                                    std::cout << std::string(crossing(i))
                                              << " => ";
                                    // std::cout << VV.col(i).transpose() << "
                                    // => ";
                                }
                            }
                            std::cout << std::endl;
                        }
                    }
                }
            }
        }
    }

    mtao::logging::debug() << "Filling cell-grid ownership";
    for (auto &&[idx, f] : mtao::iterator::enumerate(ret.m_faces)) {
        auto c = f.possible_cells([&](int idx) { return AV.at(idx).coord; });
        if (c.size() != 1) {
            std::cout << std::string(f) << std::endl;
        }
        assert(c.size() == 1);
        int grid_cell = StaggeredGrid::cell_index(*c.begin());
        ret.cell_grid_ownership[grid_cell].insert(idx);
    }
    for (auto &&[local_idx, cell_coord] :
         mtao::iterator::enumerate(ret.exterior_grid.cell_coords())) {
        ret.cell_grid_ownership[cell_grid().index(cell_coord)].insert(
            local_idx + ret.num_cutcells());
    }

    auto make_boundary_pair =
        [&](auto &&boundary_facet) -> std::optional<std::tuple<int, bool>> {
        // auto pc = boundary_facet.possible_cells(AV);
        auto pc = boundary_facet.possible_cells(
            [&](int idx) { return AV.at(idx).coord; });
        if (pc.size() != 2) {
            return {};
        }
        std::array<coord_type, 2> pca{{{{}}, {{}}}};
        std::copy(pc.begin(), pc.end(), pca.begin());
        // find the boundary cells axis
        int idx;
        for (idx = 0; idx < 2; ++idx) {
            if (pca[0][idx] != pca[1][idx]) {
                break;
            }
        }
        if (pca[0][idx] + 1 != pca[1][idx]) {
            mtao::logging::error()
                << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?";
        }
        if (pca[0][idx] < 0) {
            return {{-2, 1}};
        } else if (pca[1][idx] >= cell_shape()[idx]) {
            return {{-2, 0}};
        } else {
            int pi = cell_index(pca[0]);
            int ni = cell_index(pca[1]);
            bool pa = m_active_grid_cell_mask.get(pi);
            bool na = m_active_grid_cell_mask.get(ni);
            bool sign = (idx % 2 == 1);
            if (na ^ pa) {  // active inactive boundary
                if (pa) {
                    return {{pi, sign}};
                } else {
                    return {{ni, !sign}};
                }
            }
        }
        return {};
    };
    mtao::logging::debug() << "Making cut-edge boundaries";
    std::map<Edge, std::pair<int, bool>> diredge_to_edge_orientation;
    for (auto &&[idx, e] : mtao::iterator::enumerate(ret.m_cut_edges)) {
        e.external_boundary = make_boundary_pair(e);

        auto c = e.indices;
        diredge_to_edge_orientation[c] = {idx, false};
        std::swap(c[0], c[1]);
        diredge_to_edge_orientation[c] = {idx, true};
    }
    for (auto &&[idx, f] : mtao::iterator::enumerate(ret.m_faces)) {
        for (auto &&c : f.indices) {
            auto &fm = ret.m_face_boundary_map[idx];
            for (int j = 0; j < c.size(); ++j) {
                int k = (j + 1) % c.size();
                Edge me{{c[j], c[k]}};
                auto pr = diredge_to_edge_orientation.at(me);
                auto [eidx, sgn] = pr;
                auto e = ret.m_cut_edges[eidx].indices;
                fm.emplace(pr);
            }
        }
    }

    mtao::logging::debug() << "Original vertices: " << origV().size();
    ret.m_origV.resize(2, origV().size());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, origV().size()),
                      [&](const tbb::blocked_range<size_t> &range) {
                          for (size_t idx = range.begin(); idx != range.end();
                               ++idx) {
                              ret.m_origV.col(idx) = origV()[idx];
                          }
                      });
    ret.m_origE = data().E();
    return ret;
}

void CutCellGenerator<2>::extra_metadata(CutCellMesh<2> &mesh) const {
    mtao::data_structures::DisjointSet<int> cell_ds;
    int num_cells = mesh.num_cells();
    int num_cutfaces = mesh.num_cutfaces();
    {
        auto t = mtao::logging::profiler("region disjoint set construction",
                                         false, "profiler");

        // exterior grid pairs
        for (int i = 0; i < num_cells + mesh.num_cutedges(); ++i) {
            cell_ds.add_node(i);
        }
        for (auto &&[a, b] : mesh.exterior_grid.boundary_facet_pairs()) {
            if (a >= 0 && b >= 0) {
                cell_ds.join(a + num_cutfaces, b + num_cutfaces);
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
        for (auto &&[eidx, e] : mtao::iterator::enumerate(mesh.cut_edges())) {
            if (e.external_boundary) {
                auto [of, sgn] = *e.external_boundary;
                if (of >= 0) {
                    cell_ds.join(eidx + num_cells,
                                 mesh.exterior_grid.cell_indices().get(of) +
                                     num_cutfaces);
                }
            }
        }

        cell_ds.reduce_all();
    }

    int outside_root = -1;
    {  // find a grid-boundary cell. this is a substantial amount of computation
       // to make the "outside" the "first region"
        int boundary_idx = -1;  // a boundary face
        for (auto &&[a, b] : mesh.exterior_grid.boundary_facet_pairs()) {
            if (a == -2) {
                boundary_idx = b + num_cutfaces;
            } else if (b == -2) {
                boundary_idx = a + num_cutfaces;
            }
        }
        if (boundary_idx ==
            -1) {  // try to find a boundary edge, and then a face from it
            // search through cut-edges for something on the boundary
            for (auto &&[eidx, e] :
                 mtao::iterator::enumerate(mesh.cut_edges())) {
                if (e.external_boundary) {
                    auto [of, sgn] = *e.external_boundary;
                    if (of == -2) {
                        boundary_idx = eidx;
                    }
                }
            }
            for (auto [cid, edge_pairs] : mtao::iterator::enumerate(
                     mesh.exterior_grid.boundary_facet_pairs())) {
                auto &&[a, b] = edge_pairs;
                if (a == -2 && b >= 0) {
                    outside_root = b;
                    break;
                }
                if (b == -2 && a >= 0) {
                    outside_root = a;
                    break;
                }
            }
        } else {
            outside_root = cell_ds.get_root(boundary_idx).data;
        }
    }
    std::map<int, int> reindexer;
    if (outside_root >= 0) {
        reindexer[outside_root] = 0;
    }
    for (int i = 0; i < num_cells; ++i) {
        int root = cell_ds.get_root(i).data;
        if (root != outside_root && reindexer.find(root) == reindexer.end()) {
            reindexer[root] = reindexer.size();
        }
    }

    auto &eg_regions = mesh.exterior_grid.m_regions;
    for (int i = 0; i < mesh.exterior_grid.num_cells(); ++i) {
        eg_regions[i] =
            reindexer[cell_ds.get_root(mesh.num_cutfaces() + i).data];
    }

    for (int i = 0; i < num_cutfaces; ++i) {
        mesh.m_faces[i].region = reindexer[cell_ds.get_root(i).data];
    }
}

}  // namespace mandoline::construction
