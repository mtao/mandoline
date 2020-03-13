#include "mandoline/construction/generator2.hpp"
#include <spdlog/spdlog.h>

template <typename GridB>
void print_gridb(const GridB& g) {

    if constexpr(GridB::D == 3) 
    {
        for(int i = 0; i < g.shape(0); ++i) {
            for(int j = 0; j < g.shape(1); ++j) {
                for(int k = 0; k < g.shape(2); ++k) {
                    if(g(i,j,k)) {
                        std::cout << "o";
                    } else {
                        std::cout << ".";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

}
template <typename GridB>
void print_grid(const GridB& g) {

    if constexpr(GridB::D == 3) 
    {
        for(int i = 0; i < g.shape(0); ++i) {
            for(int j = 0; j < g.shape(1); ++j) {
                for(int k = 0; k < g.shape(2); ++k) {
                    std::cout << g(i,j,k) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

}
namespace mandoline::construction {

std::array<int, 2> smallest_ordered_edge(const std::vector<int> &v) {
    assert(v.size() >= 2);
    std::array<int, 2> min{ { v[0], v[1] } };
    for (int i = 0; i < v.size(); ++i) {
        int j = (i + 1) % v.size();
        std::array<int, 2> e{ { v[i], v[j] } };
        min = std::min(e, min);
    }
    return min;
}
template<>
auto CutCellEdgeGenerator<2>::compute_planar_hem(const std::vector<VType> &GV, const Edges &E, const GridDatab &interior_cell_mask) const -> std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> {
    ColVecs V(2, GV.size());
    for (int i = 0; i < GV.size(); ++i) {
        auto v = V.col(i);
        v = GV[i].p();
    }
    return compute_planar_hem(GV, V, E, interior_cell_mask);
}
template<>
auto CutCellEdgeGenerator<2>::compute_planar_hem(const std::vector<VType> &GV, const ColVecs &V, const Edges &E, const GridDatab &interior_cell_mask) const -> std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> {
    spdlog::warn("Doing planar hem GV,V,E,icm");
    {
        auto s = interior_cell_mask.shape();
        spdlog::warn("Interior cell mask shape: {},{}", s[0],s[1]);
    }
    print_gridb(interior_cell_mask);
    bool adaptive = interior_cell_mask.empty();
    auto ret = compute_planar_hem(V, E, interior_cell_mask);
    //the cells that vertices belong to
    std::vector<std::set<coord_type>> vertex_cells(GV.size());
    for (auto &&[v, cs] : mtao::iterator::zip(GV, vertex_cells)) {
        coord_type c = v.coord;
        cs.insert(c);
        if (v.clamped(0)) {
            coord_type cc = v.coord;
            cc[0]--;
            cs.insert(cc);
        }
        if (v.clamped(1)) {
            coord_type cc = v.coord;
            cc[1]--;
            cs.insert(cc);
        }
        if (v.is_grid_vertex()) {
            coord_type cc = v.coord;
            cc[0]--;
            cc[1]--;
            cs.insert(cc);
        }
    }
    auto &&[hem, Es] = ret;
    if (adaptive) {

        for (int i = 0; i < E.cols(); ++i) {
            auto e = E.col(i);
            auto a = GV[e(0)];
            auto b = GV[e(1)];
            for (auto &&[i, v] : mtao::iterator::enumerate(a.mask() & b.mask())) {
                if (v) {
                    int vv = *v;
                    if (vv == 0 || vv == cell_shape()[i]) {
                        Es.emplace(Edge{ { e(0), e(1) } });
                    }
                }
            }
        }
    }


    //the cells that each single=-curve face belongs to
    mtao::map<int, std::set<coord_type>> cell_coords;

    auto ci = hem.cell_indices();
    auto vi = hem.vertex_indices();
    //for every halfedge get the cell masks
    for (int i = 0; i < hem.size(); ++i) {
        int cell = ci(i);
        if (auto it = cell_coords.find(cell); it == cell_coords.end()) {
            cell_coords[cell] = vertex_cells[vi(i)];
        } else {

            auto &coords = it->second;
            auto &&vc = vertex_cells[vi(i)];
            std::set<coord_type> gs;
            std::set_intersection(coords.begin(), coords.end(), vc.begin(), vc.end(), std::inserter(gs, gs.end()));
            coords = gs;
        }
    }

    auto cells_map = hem.cells_map();
    auto cells_halfedge_map = hem.cell_halfedges_map();
    mtao::map<coord_type, std::set<int>> coord_cells;
    for (auto &&[cc, cs] : cell_coords) {
        for (auto &&c : cs) {
            coord_cells[c].insert(cc);
        }
    }
    auto areas = hem.signed_areas(V);
    for (auto &&[c, cs] : coord_cells) {
        if (cs.size() < 2) {
            continue;
        }
        std::set<int> ches;
        for (auto &&c : cs) {
            if (c >= 0) {
                ches.emplace(cells_halfedge_map.at(c));
            }
        }
        hem.tie_nonsimple_cells(V, ches);
    }
    return ret;
}


template<>
auto CutCellEdgeGenerator<2>::compute_planar_hem(const ColVecs &V, const Edges &E, const GridDatab &interior_cell_mask) const -> std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> {

    spdlog::warn("Doing planar hem V,E,icm");
    print_gridb(interior_cell_mask);
    auto t = mtao::logging::profiler("computing planar hem", false, "profiler");
    using namespace mtao::geometry::mesh;

    bool adaptive = interior_cell_mask.empty();
    EmbeddedHalfEdgeMesh<double, 2> ehem;
    {
        //auto t = mtao::logging::timer("Making halfedge mesh");
        auto VV = V;
        VV.row(0) = V.row(1);
        VV.row(1) = V.row(0);
        std::set<Edge> E2;
        for (int i = 0; i < E.cols(); ++i) {
            auto e_ = E.col(i);
            Edge e{ { e_(0), e_(1) } };
            std::sort(e.begin(), e.end());
            E2.emplace(e);
        }
        ehem = EmbeddedHalfEdgeMesh<double, 2>::from_edges(VV, mtao::eigen::stl2eigen(E2));
        //ehem = EmbeddedHalfEdgeMesh<double,2>::from_edges(VV,E);
    }
    //auto ehem = EmbeddedHalfEdgeMesh<double,2>::from_edges(ret.vertices(),ret.cut_edges);


    {
        //auto t = mtao::logging::timer("Generating topology");
        ehem.make_topology();
    }

    std::set<Edge> boundary_edges;


    auto print_coord = [&](auto &&c) {
        std::cout << "(";
        for (auto &&v : c) {
            std::cout << v << ",";
        }
        std::cout << ")";
    };
    { // create the boundary edge data

        if (!adaptive) {
            //auto t = mtao::logging::timer("actual looping in mask");
            coord_type shape = interior_cell_mask.shape();
            for (auto &&s : shape) {
                s++;
            }
            mtao::geometry::grid::utils::multi_loop(shape, [&](const coord_type &c) {
                for (int i = 0; i < 2; ++i) {
                    auto cc = c;
                    cc[i]--;


                    bool c_valid = (!interior_cell_mask.valid_index(c)) || interior_cell_mask(c);
                    bool cc_valid = (!interior_cell_mask.valid_index(cc)) || interior_cell_mask(cc);
                    if (c_valid ^ cc_valid) {
                        Edge isarr = StaggeredGrid::edge(1 - i, c);
                        std::sort(isarr.begin(), isarr.end());
                        boundary_edges.insert(isarr);
                    }
                }
            });
        } else {
        }
    }
    return { ehem, boundary_edges };
}


}// namespace mandoline::construction
