#include "mandoline/construction/generator3.hpp"
#include <mtao/geometry/mesh/boundary_facets.h>
#include <mtao/geometry/mesh/face_normals.hpp>
#include <mtao/geometry/mesh/boundary_elements.h>
#include <mtao/geometry/trigonometry.hpp>
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/geometry/mesh/dual_edges.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/iterator/enumerate.hpp>
#include <mtao/logging/logger.hpp>
#include <mtao/reindexer.hpp>
#include "mandoline/construction/subgrid_transformer.hpp"
#include <variant>
using namespace mtao::iterator;
using namespace mtao::logging;
namespace mandoline::construction {
bool CutCellGenerator<3>::AxisHEMData::is_boundary_cell(const std::vector<int> &C) const {
    for (int i = 0; i < C.size(); ++i) {
        Edge e{ { C[i], C[(i + 1) % C.size()] } };
        std::sort(e.begin(), e.end());
        if (boundary_edges.find(e) == boundary_edges.end()) {
            return false;
        }
    }
    return true;
}


void CutCellGenerator<3>::bake_faces() {

    {
        mtao::ColVecs3d V(3, data().V().size());
        for (int i = 0; i < V.cols(); ++i) {
            V.col(i) = data().V()[i].p();
        }
        origN = mtao::geometry::mesh::face_normals(V, data().F());
    }
    m_cut_faces = data().cut_faces();

    {
        auto t = mtao::logging::profiler("caching axial primal faces", false, "profiler");
        std::vector<bool> active(data().nF());
        auto &&ti = data().triangle_intersections();
        std::transform(ti.begin(), ti.end(), active.begin(), [](auto &&ti) {
            return ti.active();
        });
        cut_cell_to_primal_map.clear();
        for (auto &&[i, f] : mtao::iterator::enumerate(m_cut_faces)) {
            cut_cell_to_primal_map[i] = f.parent_fid;
            if (active[f.parent_fid]) {
                for (auto &&[i, c] : mtao::iterator::enumerate(f.mask())) {
                    if (c) {
                        double Ni = origN.col(f.parent_fid)(i);
                        if (Ni > 0) {
                            axial_primal_faces[i].insert(smallest_ordered_edge(f.indices));
                        } else {
                            std::vector<int> IR(f.indices.rbegin(), f.indices.rend());
                            axial_primal_faces[i].insert(smallest_ordered_edge(IR));
                        }
                    }
                }
            }
        }
    }
    int face_size = m_cut_faces.size();


    using Grid2 = mtao::geometry::grid::StaggeredGrid<double, 2>;
    auto subgrid = [this](auto &&g, int d) {
        int n0 = (d + 1) % 3;
        int n1 = (d + 2) % 3;
        auto s = vertex_shape();
        std::array<int, 2> shape{ { s[n0], s[n1] } };
        mtao::Vec2d dx(g.dx()(n0), g.dx()(n1));
        mtao::Vec2d origin(g.origin()(n0), g.origin()(n1));
        return CutCellEdgeGenerator<2>(shape);
    };
    std::array<CutCellEdgeGenerator<2>, D> grids{ { subgrid(*this, 0), subgrid(*this, 1), subgrid(*this, 2) } };

    auto AC = active_cells();
    auto AE = axial_edges();
    //#pragma omp parallel for

    auto t = mtao::logging::profiler("Flagging active cells");
    for (auto &&c : AC) {
        for (int i = 0; i < D; ++i) {
            int p1 = (i + 1) % 3;
            int p2 = (i + 2) % 3;
            auto &ahdata = axis_hem_data[i];
            Edge c2{ { c[p1], c[p2] } };
            if (c[p1] < 0 || c[p1] >= cell_shape()[p1]) {
                continue;
            }
            if (c[p2] < 0 || c[p2] >= cell_shape()[p2]) {
                continue;
            }
            auto flag = [&](int idx) {
                if (auto it = ahdata.find(idx); it == ahdata.end()) {
                    ahdata[idx].active_grid_cell_mask = AxisHEMData::GridDatab::Constant(true, grids[i].cell_shape());
                }
                auto &ahd = ahdata[idx];
                ahd.active_grid_cell_mask(c2) = false;
            };
            if (c[i] >= 0) {
                flag(c[i]);
            }
            if (c[i] + 1 < vertex_shape()[i]) {
                flag(c[i] + 1);
            }
        }
    }

    auto V = all_GV();


    mtao::ColVectors<double, 2> subV(2, V.cols());
    std::vector<Vertex<3>> VV(num_vertices());
    std::vector<Vertex<2>> subVV(num_vertices());
    for (int i = 0; i < num_vertices(); ++i) {
        VV[i] = grid_vertex(i);
    }
    std::string timer_str = "Flagging active cells 0";
    for (int i = 0; i < 3; ++i) {
        int dim = i;
        auto &&ahdata = axis_hem_data[i];
        auto &&ae = AE[i];
        auto &&g = grids[i];
        //for(auto&& [dim,ahdata,ae,g_]: mtao::iterator::enumerate(axis_hem_data, AE,grids)) {
        //auto& g = g_;
        SubGridTransformer subtran(g, *this, dim, 0);
        //subV.row(0) = V.row((i+1)%3);
        //subV.row(1) = V.row((i+2)%3);

        for (int i = 0; i < num_vertices(); ++i) {
            subVV[i] = subtran.downgrade(VV[i]);
        }

        {
            auto t = mtao::logging::timer(timer_str);
            timer_str[timer_str.size() - 1] += 1;

            std::vector<int> aiinds(ae.size());
            std::transform(ae.begin(), ae.end(), aiinds.begin(), [](auto &&pr) {
                return pr.first;
            });
            auto it = aiinds.begin();

#pragma omp parallel for
            for (it = aiinds.begin(); it < aiinds.end(); it++) {
                int cidx = *it;
                auto &E = ae[cidx];
                //auto&& [cidx,E] = *it;
                //for(auto&& [cidx,E]: ae) {
                CoordType coord;
                coord[dim] = cidx;
                auto &ahd = ahdata[cidx];
                std::transform(E.begin(), E.end(), std::inserter(ahd.edges, ahd.edges.end()), [](auto s) {
                    if (s[0] > s[1]) {
                        std::swap(s[0], s[1]);
                    }
                    return s;
                });
                //ahd.edges.insert(E.begin(),E.end());
                std::cout << "AHD size: " << ahd.edges.size() << std::endl;

                std::set<Edge> bedges;
                std::tie(ahd.hem, bedges) = g.compute_planar_hem(subVV, mtao::eigen::stl2eigen(ahd.edges), ahd.active_grid_cell_mask);
                std::transform(bedges.begin(), bedges.end(), std::inserter(ahd.boundary_edges, ahd.boundary_edges.end()), [&, dim = dim](Edge e) -> Edge {
                    for (auto &&idx : e) {

                        auto c = g.vertex_unindex(idx);
                        coord[(dim + 1) % 3] = c[0];
                        coord[(dim + 2) % 3] = c[1];
                        idx = vertex_index(coord);
                    }
                    return e;
                });
            }
            for (auto &&[cidx, E] : ae) {
                CoordType coord;
                auto &ahd = ahdata[cidx];
                auto fa = ahd.hem.cell_indices().array();
                fa = (fa >= 0).select(fa + face_size, -1);
                face_size = ahd.hem.cell_indices().maxCoeff() + 1;
            }
        }
    }
    /*
               for(auto&& axis_hem: axis_hem_data) {
               for(auto&& [i,ahd]: axis_hem) {
               for(auto&& e: ahd.boundary_edges) {
               for(auto&& v: e) {
               boundary_vertices.insert(v);
               }
               }
               }
               }
               */
    compute_faces();
}


auto ring_to_edges(const std::vector<int> &F) {
    std::set<std::array<int, 2>> Eset;
    for (int i = 0; i < F.size(); ++i) {
        int a = F[i];
        int b = F[(i + 1) % F.size()];
        if (a != b) {
            if (a < b) {
                Eset.emplace(std::array<int, 2>{ { a, b } });
            } else {
                Eset.emplace(std::array<int, 2>{ { b, a } });
            }
        }
    }
    return mtao::eigen::stl2eigen(Eset);
}


void CutCellGenerator<3>::compute_faces() {
    m_faces.clear();

    {
        auto t = mtao::logging::timer("Computing input mesh faces");
        compute_faces_vertex();
    }
    for (int i = 0; i < 3; ++i) {
        auto t = mtao::logging::timer("Computing input mesh faces");
        compute_faces_axis(i);
    }
    for (auto &&[fidx, f] : m_faces) {
        auto &&I = f.indices;
        //if(f.is_mesh_face())
        //{
        //std::cout << "Mesh Face: " << fidx << std::endl;
        //} else {
        //std::cout << "Grid Face: " << fidx << std::endl;
        //}
        //for(auto&& l: I)
        //{
        //    for(auto&& v: l)
        //    {
        //        std::cout << v << ":";
        //    }
        //    std::cout << std::endl;
        //    if(l.size() > 4)
        //    {
        //        for(auto&& v: l)
        //        {
        //            std::cout << grid_info(v) << std::endl;
        //        }
        //    }
        //}
        //std::cout << std::endl;
        if (f.mask().count() == 1) {
            //auto&& loop = *I.begin();
            auto pc = possible_cells(I);
            {
                //                    if(loop.size() == 4) {
                if (pc.size() != 2) {
                    continue;
                }
                std::array<CoordType, 2> pca;
                std::copy(pc.begin(), pc.end(), pca.begin());
                //find the boundary cells axis
                int idx;
                for (idx = 0; idx < 3; ++idx) {
                    if (pca[0][idx] != pca[1][idx]) {
                        break;
                    }
                }
                if (pca[0][idx] + 1 != pca[1][idx]) {
                    std::cout << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?" << std::endl;
                }
                if (pca[0][idx] < 0) {
                    f.external_boundary = { -1, 1 };
                } else if (pca[1][idx] >= cell_shape()[idx]) {
                    f.external_boundary = { -1, 0 };
                } else if (I.size() == 1 && I.begin()->size() == 4) {//should always be a grid cell!
                    int pi = cell_index(pca[0]);
                    int ni = cell_index(pca[1]);
                    bool pa = m_active_grid_cell_mask.get(pi);
                    bool na = m_active_grid_cell_mask.get(ni);
                    bool sign = (idx%2==1);
                    if (na ^ pa) {//active inactive boundary
                        if (pa) {
                            f.external_boundary = { pi, sign };
                        } else {
                            f.external_boundary = { ni, !sign };
                        }
                    }
                }
            }
        }
    }
}
void CutCellGenerator<3>::compute_faces_vertex() {
    auto t = mtao::logging::timer("Computing input mesh faces with normals ");
    std::array<double, 3> N;
    for (auto &&[i, f] : mtao::iterator::enumerate(m_cut_faces)) {
        if (f.indices.size() > 2) {// && f.mask().count() < 2) {

            //std::cout << "Normal dots: " << f.parent_fid << ": " << origN.col(f.parent_fid).dot(data().m_triangle_intersections[f.parent_fid].N()) << std::endl;
            m_faces[i] = CutFace<D>(f, origN.col(f.parent_fid));
            mesh_face_indices.insert(i);
        } else {
            std::cout << "Cut faces should always be  real faces!" << std::endl;
        }
    }
}
void CutCellGenerator<3>::compute_faces_axis(int idx) {
    //axial half edge mesh data
    auto &ahdata = axis_hem_data[idx];

    //active axial indices
    std::vector<int> aiinds(ahdata.size());
    std::transform(ahdata.begin(), ahdata.end(), aiinds.begin(), [](auto &&pr) {
        return pr.first;
    });
    auto it = aiinds.begin();
    std::vector<mtao::map<int, CutFace<D>>> faces_vec(ahdata.size());

#pragma omp parallel for
    for (it = aiinds.begin(); it < aiinds.end(); it++) {
        int cidx = *it;
        auto &ahd = ahdata.at(cidx);
        faces_vec[std::distance(aiinds.begin(), it)] = compute_faces_axis(idx, cidx);
    }

    for (auto &&fcs : faces_vec) {
        for (auto &&[id, fs] : fcs) {
            axis_face_indices[idx].insert(id);
        }
        m_faces.insert(fcs.begin(), fcs.end());
    }
}
mtao::map<int, CutFace<3>> CutCellGenerator<3>::compute_faces_axis(int idx, int cidx) const {

    mtao::map<int, CutFace<D>> faces;
    auto &ahdata = axis_hem_data[idx];
    auto &ahd = ahdata.at(cidx);

    auto V = all_GV();

    mtao::ColVectors<double, 2> VV(2, V.cols());
    VV.row(0) = V.row((idx + 1) % 3);
    VV.row(1) = V.row((idx + 2) % 3);

    mtao::Vec3d N = mtao::Vec3d::Unit((idx + 1) % 3).cross(mtao::Vec3d::Unit((idx + 2) % 3));

    auto mesh_faces = ahd.hem.cells_multi_component_map();

    auto vols = ahd.hem.signed_areas(VV);
    for (auto &&[i, v] : mesh_faces) {
        CutFace<D> F;
        F.id = Edge{ { idx, cidx } };
        auto add_loop = [&](const std::vector<int> &v) {
            if (v.size() > 2) {
                if (axial_primal_faces[idx].find(smallest_ordered_edge(v)) == axial_primal_faces[idx].end()) {
                    F.indices.emplace(v);
                }
            }
        };
        for (auto &v : v) {// for every boundary cell figure out if this is a boundary stencil
            if (ahd.is_boundary_cell(v)) {
                auto pc = possible_cells(v);
                if (pc.empty()) {
                    continue;
                }
                std::array<CoordType, 2> pca;
                std::copy(pc.begin(), pc.end(), pca.begin());
                if (pca[0][idx] + 1 != pca[1][idx]) {
                    std::cout << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?" << std::endl;
                }
                std::array<int, 2> ind;
                ind[0] = pca[0][(idx + 1) % 3];
                ind[1] = pca[0][(idx + 2) % 3];
                if (!ahd.active_grid_cell_mask.valid_index(ind)) {
                    std::cout << "Invalid index???? " << std::endl;
                }
                if (!ahd.active_grid_cell_mask(ind)) {
                    if (vols.at(i) >= 0) {
                        add_loop(v);
                    }
                }
            } else {
                if (vols.at(i) >= 0) {
                    add_loop(v);
                }
            }
        }
        if (!F.indices.empty()) {

            F.N = N;
            F.coord_mask<D>::operator=(face_mask(F.indices));
            faces[i] = std::move(F);
            //std::cout << std::string(faces[i]) << std::endl;
        }
    }

    /*
                //debug when i was only generating quads
                bool bad = false;
                for(auto&& [fidx,f]: faces) 
                {
                    if(f.indices.begin()->size() > 4)
                    {
                        bad = true;
                    }
                }
                if(bad)
                {
                    std::cout << "A bunch of edges weren't used due to some failure. " << std::endl;
                    std::cout << " we only got " << faces.size() << "things" << std::endl;
                    std::cout << "Edges on " << idx << "," << cidx << std::endl;
                    std::cout << mtao::eigen::stl2eigen(ahd.edges) << std::endl;
                    for(auto&& [a,b]: ahd.edges)
                    {
                        std::cout << grid_info(a) << " => " << grid_info(b) << std::endl;
                    }
                    std::cout << std::endl;
                }
                */
    return faces;
}

mtao::ColVecs3i CutCellGenerator<3>::faceMap_to_faces(mtao::map<int, mtao::ColVecs3i> &fm) {

    int cols = 0;
    for (auto &&[i, f] : fm) {
        cols += f.cols();
    }
    mtao::ColVecs3i F(3, cols);

    int col = 0;
    for (auto &&[i, f] : fm) {
        F.block(0, col, f.rows(), f.cols()) = f;
        col += f.cols();
    }
    return F;
}
}// namespace mandoline::construction
