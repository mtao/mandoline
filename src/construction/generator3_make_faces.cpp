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
                            auto e = smallest_ordered_edge(f.indices);
                            //std::swap(e[0],e[1]);
                            axial_primal_faces[i].insert(e);
                        } else {
                            auto e = smallest_ordered_edge_reverse(f.indices);
                            //std::swap(e[0],e[1]);
                            axial_primal_faces[i].insert(e);
                        }
                    }
                }
            }
        }
    }
    int face_size = m_cut_faces.size();


    //{
    //    auto s = StaggeredGrid::cell_shape();
    //    spdlog::warn("cell shape: {},{},{}", s[0],s[1],s[2]);
    //}
    //{
    //    auto s = active_grid_cell_mask().shape();
    //    spdlog::warn("Interior cell mask shape: {},{},{}", s[0],s[1],s[2]);
    //}
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

    auto activeCells = active_cells();
    auto axialEdges = axial_edges();
    //#pragma omp parallel for

    auto t = mtao::logging::profiler("Flagging active cells");
    for (auto &&c : activeCells) {
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
        auto &&axialEdges_dim = axialEdges[i];
        auto &&g = grids[i];
        //for(auto&& [dim,ahdata,axialEdges_dim,g_]: mtao::iterator::enumerate(axis_hem_data, axialEdges,grids)) {
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

            std::vector<int> axial_edge_indices(axialEdges_dim.size());
            std::transform(axialEdges_dim.begin(), axialEdges_dim.end(), axial_edge_indices.begin(), [](auto &&pr) {
                return pr.first;
            });
            auto it = axial_edge_indices.begin();

#pragma omp parallel for
            for (it = axial_edge_indices.begin(); it < axial_edge_indices.end(); it++) {
                int cidx = *it;
                auto &E = axialEdges_dim[cidx];
                //auto&& [cidx,E] = *it;
                //for(auto&& [cidx,E]: axialEdges_dim) {
                coord_type coord;
                coord[dim] = cidx;
                auto &ahd = ahdata[cidx];
                std::transform(E.begin(), E.end(), std::inserter(ahd.edges, ahd.edges.end()), [](auto s) {
                    if (s[0] > s[1]) {
                        std::swap(s[0], s[1]);
                    }
                    //std::cout << s[0] << ":" << s[1] << " ";
                    return s;
                });
                //std::cout << std::endl;
                //ahd.edges.insert(E.begin(),E.end());
                //spdlog::info("AHD size {}", ahd.edges.size());
                //std::cout << "AHD size: " << ahd.edges.size() << std::endl;

                std::set<Edge> bedges;
                //if(dim == 2 && cidx == 1) {
                //    std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << std::endl;
                //    std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << std::endl;
                //    std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << std::endl;
                //    std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << std::endl;
                //}
                std::tie(ahd.hem, bedges) = g.compute_planar_hem(subVV, mtao::eigen::stl2eigen(ahd.edges), ahd.active_grid_cell_mask);
                //if(dim == 2 && cidx == 1) {
                //    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
                //    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
                //    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
                //    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
                //}
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
            for (auto &&[cidx, E] : axialEdges_dim) {
                coord_type coord;
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
    // the sign must be determined by the cell it shares, we just write to external_boundary and let generate_celll fill in the signs
    for (auto &&[fidx, f] : m_faces) {
        auto &&I = f.indices;
        if (f.mask().count() == 1) {
            //auto&& loop = *I.begin();
            // external boundary check: if this face is on a boundary
            auto pc = possible_cells(I);
            {
                //                    if(loop.size() == 4) {
                if (pc.size() != 2) {
                    continue;
                }
                std::array<coord_type, 2> pca;
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
                    f.external_boundary = { -2, 0 };
                } else if (pca[1][idx] >= cell_shape()[idx]) {
                    f.external_boundary = { -2, 0 };
                } else if (I.size() == 1 && I.begin()->size() == 4) {//should always be a grid cell!
                    bool pa = is_active_cell(pca[0]);
                    bool na = is_active_cell(pca[1]);
                                //spdlog::info("boundary region {} {} {} => {} {} {}"
                                //        , pca[0][0]
                                //        , pca[0][1]
                                //        , pca[0][2]
                                //        , pca[1][0]
                                //        , pca[1][1]
                                //        , pca[1][2]
                                //        );
                    if (na ^ pa) {//active inactive boundary
                        if (pa) {
                            if(is_valid_cell_coord(pca[0])) {
                                int pi = cell_index(pca[0]);
                                f.external_boundary = { pi, 0};
                            } else {
                                f.external_boundary = { -2, 0};
                            }
                        } else {
                            if(is_valid_cell_coord(pca[1])) {
                                int ni = cell_index(pca[1]);
                                f.external_boundary = { ni, 0};
                            } else {
                                f.external_boundary = { -2, 0};
                            }
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

    //pull the keys out so we can parallel for loop over them
    std::vector<int> axial_edge_indices(ahdata.size());
    std::transform(ahdata.begin(), ahdata.end(), axial_edge_indices.begin(), [](auto &&pr) {
        return pr.first;
    });
    auto it = axial_edge_indices.begin();
    std::vector<mtao::map<int, CutFace<D>>> faces_vec(ahdata.size());

#pragma omp parallel for
    for (it = axial_edge_indices.begin(); it < axial_edge_indices.end(); it++) {
        int cidx = *it;
        auto &ahd = ahdata.at(cidx);
        auto r = compute_faces_axis(idx, cidx);
        //spdlog::info("compute faces axis {} level {} got size {}",idx,cidx,r.size());
        faces_vec[std::distance(axial_edge_indices.begin(), it)] = std::move(r);
    }

    for (auto &&fcs : faces_vec) {
        for (auto &&[id, fs] : fcs) {
            axis_face_indices[idx].insert(id);
            //std::cout << std::string(fs) << std::endl;
        }
        m_faces.insert(fcs.begin(), fcs.end());
    }
}
mtao::map<int, CutFace<3>> CutCellGenerator<3>::compute_faces_axis(int idx, int cidx) const {
    //spdlog::error("Compute_face_axis({},{})",idx,cidx);

    mtao::map<int, CutFace<D>> faces;
    auto &ahdata = axis_hem_data[idx];
    auto &ahd = ahdata.at(cidx);

    auto V = all_GV();

    mtao::ColVectors<double, 2> VV(2, V.cols());
    VV.row(0) = V.row((idx + 1) % 3);
    VV.row(1) = V.row((idx + 2) % 3);

    mtao::Vec3d N = mtao::Vec3d::Unit((idx + 1) % 3).cross(mtao::Vec3d::Unit((idx + 2) % 3));

    auto mesh_faces = ahd.hem.cells_multi_component_map();
    spdlog::info("AHD {},{} got {} faces", idx,cidx,mesh_faces.size());

    auto vols = ahd.hem.signed_areas(VV);
    //std::cout << "Forbidden edges: ";
    //for(auto&& [a,b]: axial_primal_faces[idx]) {
    //    std::cout << a << ":" << b << " ";
    //}
    //std::cout << std::endl;
    for (auto &&[i, v] : mesh_faces) {
        CutFace<D> F;
        F.id = Edge{ { idx, cidx } };
        auto add_loop = [&](const std::vector<int> &v) {

            if (v.size() > 2) {
                //std::copy(v.begin(),v.end(),std::ostream_iterator<int>(std::cout,","));
                //std::cout << ": ";
                auto e = smallest_ordered_edge(v);
                //auto er = smallest_ordered_edge_reverse(v);
                //std::cout << "Trying edge " << e[0] << ":" << e[1] << " rev: " << er[0] << ":" << er[1] << std::endl;
                if (axial_primal_faces[idx].find(e) == axial_primal_faces[idx].end()) {
                    //std::vector<int> vv(v.size());

                    //std::copy(v.rbegin(),v.rend(),vv.begin());

                    F.indices.emplace(std::move(v));
                }
            }
        };
        //std::cout << i << ") " << vols.at(i) << ": ";
        //for (auto &v : v) {
        //    std::copy(v.begin(),v.end(),std::ostream_iterator<int>(std::cout,","));
        //    std::cout << " ";
        //}
        //std::cout << std::endl;
        for (auto &v : v) {// for every boundary cell figure out if this is a boundary stencil
            //spdlog::info("Is boundary face?: {}", ahd.is_boundary_cell(v));
            if (ahd.is_boundary_cell(v)) {
                auto pc = possible_cells(v);
                if (pc.empty()) {
                    continue;
                }
                std::array<coord_type, 2> pca;
                std::copy(pc.begin(), pc.end(), pca.begin());
                if (pca[0][idx] + 1 != pca[1][idx]) {
                    spdlog::error("SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?");
                }
                std::array<int, 2> ind;
                ind[0] = pca[0][(idx + 1) % 3];
                ind[1] = pca[0][(idx + 2) % 3];
                if (!ahd.active_grid_cell_mask.valid_index(ind)) {
                    spdlog::error("Invalid index");
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
