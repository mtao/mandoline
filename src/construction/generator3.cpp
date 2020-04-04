#include "mandoline/construction/generator3.hpp"
#include <mtao/geometry/mesh/boundary_facets.h>
#include <mtao/geometry/mesh/boundary_elements.h>
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/colvector_loop.hpp>
#include <mtao/geometry/mesh/dual_edges.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/iterator/enumerate.hpp>
#include <mtao/logging/logger.hpp>
#include "mandoline/construction/subgrid_transformer.hpp"
#include <variant>
#include "mandoline/construction/cell_collapser.hpp"
#include "mandoline/construction/adaptive_grid_factory.hpp"
using namespace mtao::iterator;
using namespace mtao::logging;


namespace mandoline::construction {


CutCellGenerator<3>::~CutCellGenerator() {}

template<>
CutCellMesh<3> CutCellEdgeGenerator<3>::generate() const {
    CutCellMesh<3> ccm = generate_edges();
    //make edge stuff
    return ccm;
}

CutCellMesh<3> CutCellGenerator<3>::generate() const {


    CutCellMesh<3> ccm = CCEG::generate();
    ccm.m_folded_faces = folded_faces;
    //extra_metadata(ccm);
    ccm.m_faces.clear();
    ccm.m_faces.resize(faces().size());
    if (adaptive) {
        if (!adaptive_grid) {
            warn() << "Adaptive grid doesn't exist, is there no room for it?";
        } else {
            ccm.m_adaptive_grid = *adaptive_grid;
            if (adaptive_grid_regions) {
                ccm.m_adaptive_grid_regions = *adaptive_grid_regions;
            }
        }
    }

    //ccm.triangulated_cut_faces.resize(faces().size());
    mtao::map<int, int> reindexer;
    int max = m_cut_faces.size();
    for (auto &&[i, ff] : mtao::iterator::enumerate(faces())) {
        auto &&[idx, f] = ff;
        reindexer[idx] = i;
    }
    for (auto &&[i, f] : m_faces) {
        ccm.m_faces[reindexer.at(i)] = f;
    }

    // reindex a full set of indices
    auto redx = [&](const std::set<int> &i) {
        std::set<int> ret;
        std::transform(i.begin(), i.end(), std::inserter(ret, ret.end()), [&](int idx) -> int { return reindexer.at(idx); });
        return ret;
    };

    for (auto &&mfi : mesh_face_indices) {
        int idx = reindexer.at(mfi);
        auto &&cutface = m_cut_faces[mfi];
        mtao::ColVecs3d B(3, cutface.size());
        for (auto &&[idx, i] : mtao::iterator::enumerate(cutface.indices)) {
            if (i < grid_vertex_size()) {
                B.col(idx) = data().get_face_bary(cutface.parent_fid, grid_vertex(i));
            } else {
                B.col(idx) = data().get_face_bary(cutface.parent_fid, crossing(i));
            }
        }
        ccm.m_mesh_cut_faces[idx] = BarycentricTriangleFace{ std::move(B), cutface.parent_fid };
    }
    //ccm.mesh_faces = redx(mesh_face_indices);
    for (auto &&[a, b] : mtao::iterator::zip(ccm.m_axial_faces, axis_face_indices)) {
        a = redx(b);
    }
    /*
           if(!performance) {
           std::cout << "Making triangulated faces" << std::endl;
           for(auto&& [i,ff]: faces_triangulated_map()) {

           if(reindexer.find(i) != reindexer.end()) {
           ccm.triangulated_cut_faces[reindexer.at(i)] = ff;
           }
           }
           }
           */

    ccm.m_cells.clear();
    ccm.m_cells.resize(cell_boundaries.size());

    for (auto &&[a, b] : mtao::iterator::zip(cell_boundaries, ccm.m_cells)) {
        b.index = a.index;
        b.region = a.region;
        std::set<int> inds;
        for (auto &&[i, j] : a) {
            int fidx = reindexer.at(i);
            b[fidx] = j;
            // pos_cells currently fails with folded faces. lets give up on them for now
            if (!ccm.is_folded_face(fidx)) {
                inds.insert(fidx);
            }
        }
        //std::cout << std::string(b) << std::endl;
        auto pos_cells = possible_cells_cell(inds, ccm.faces());
        if(pos_cells.size() == 1) {
        b.grid_cell = *pos_cells.begin();
        } else if (pos_cells.size() == 0) {
            //std::cout << "CELL: " << a.index << std::endl;
            //std::cout << "No possible cells!?!?!" << std::endl;
            //for (auto &&ind : inds) {
            //    std::cout << "Face " << ind << ": ";
            //    if (ccm.is_folded_face(ind)) {
            //        std::cout << "Folded: ";
            //    }
            //    auto &face = ccm.faces()[ind];
            //    auto pc = possible_cells(face.indices);
            //    std::cout << std::string(ccm.faces()[ind]) << ":::";
            //    for (auto &&c : possible_cells(ccm.faces()[ind].indices)) {
            //        std::cout << c[0] << ":";
            //        std::cout << c[1] << ":";
            //        std::cout << c[2] << " ";
            //    }
            //    for (auto &&c : ccm.faces()[ind].indices) {
            //        for (auto &&vi : c) {
            //            std::cout << std::string(grid_vertex(vi)) << " ";
            //        }
            //    }

            //    std::cout << std::endl;
            //}
            //std::cout << std::endl;
        } else {
            mtao::logging::warn()<< "Degenerate cell";
        }
    }
    ccm.m_origV.resize(3, origV().size());
    for (int i = 0; i < origV().size(); ++i) {
        ccm.m_origV.col(i) = origV()[i];
    }
    ccm.m_origE = data().E();
    ccm.m_origF = data().F();


    return ccm;
}


void CutCellGenerator<3>::extra_metadata(CutCellMesh<3> &mesh) const {
}


void CutCellGenerator<3>::clear() {
    cut_cell_to_primal_map.clear();
    origN = {};
    axis_hem_data = {};
    axial_primal_faces = {};
    m_newF = {};
    m_faces.clear();
    mesh_face_indices.clear();
    axis_face_indices = {};
    folded_faces.clear();
    cell_boundaries.clear();
    boundary_vertices.clear();
    boundary_faces.clear();

    CCEG::clear();
}

void CutCellGenerator<3>::bake() {


    auto t2 = mtao::logging::profiler("general bake", false, "profiler");
    CCEG::bake();


    //In order to catch odd things like vertices touching stencils we have to relearn the stencil with the faces

    update_active_grid_cell_mask();


    auto t = mtao::logging::profiler("grid bake cells", false, "profiler");
    bake_cells();
}
void CutCellGenerator<3>::update_active_grid_cell_mask() {
    return;
    GridDatab mask = GridDatab::Constant(true, StaggeredGrid::cell_shape());
    auto &old_mask = m_active_grid_cell_mask;

    for (auto &&[fidx, face] : faces()) {
        auto pc = possible_cells(face.indices);
        if (face.is_axial_face()) {
            //dont only remove things from the mask, don't add
            //therefore, if we had an external boundary before it should sitll be one
            assert(pc.size() == 2);
            if (pc.size() != 2) {
                warn() << "Axial faces should have two possible cells!: face " << fidx;
            }
            if (face.external_boundary) {
                std::array<coord_type, 2> pca;
                std::copy(pc.begin(), pc.end(), pca.begin());
                int idx;
                for (idx = 0; idx < 3; ++idx) {
                    if (pca[0][idx] != pca[1][idx]) {
                        break;
                    }
                }
                if (pca[0][idx] + 1 != pca[1][idx]) {
                    std::cout << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?" << std::endl;
                }
                bool below = std::get<1>(*face.external_boundary);
                auto choice = below ? pca[0] : pca[1];
                if (mask.valid_index(choice)) {
                    mask(choice);
                }
            } else {
                for (auto &&c : pc) {
                    if (mask.valid_index(c)) {

                        mask(c) = false;
                    }
                }
            }
        } else {
            for (auto &&c : pc) {
                if (mask.valid_index_(c)) {

                    mask(c) = false;
                }
            }
        }
    }
    old_mask = mask;
}
void CutCellGenerator<3>::bake_cells() {
    {
        auto t = mtao::logging::profiler("cell collapser", false, "profiler");
        CellCollapser cc(m_faces);
        auto V = all_GV();

        //auto [a,b] = adaptive_grid->compute_edges(adaptive_level);
        //for(auto&& [hf,bd]: cc.m_halfface_to_cell) {
        //    auto&& [face_index, edge] = hf;
        //    auto&& [a,b] = edge;
        //    auto&& [bs,sgn] = bd;
        //    spdlog::warn("{}: ({},{}) => {}({})", face_index,a,b,bs,sgn);
        //}

        cc.bake(V);
        //spdlog::error("Baked!");
        //for(auto&& [hf,bd]: cc.m_halfface_to_cell) {
        //    auto&& [face_index, edge] = hf;
        //    auto&& [a,b] = edge;
        //    auto&& [bs,sgn] = bd;
        //    spdlog::warn("{}: ({},{}) => {}({})", face_index,a,b,bs,sgn);
        //}
        folded_faces = cc.folded_faces();
        //m_normal_faces = cc.m_faces;
        boundary_vertices.clear();
        //for(int i = 0; i < StaggeredGrid::vertex_size(); ++i) {
        //    boundary_vertices.insert(i);
        //}

        //cc.remove_boundary_cells_from_vertices(boundary_vertices);
        cc.remove_boundary_cells();
        //even though max index value could be vertex_shape in size, the max value is cell_Shape
        cc.remove_grid_boundary_cells(StaggeredGrid::cell_shape());

        if constexpr (false) {
            auto Vs = vertices();
            std::map<int, double> vol;
            for (auto &&[i, f] : m_faces) {
                vol[i] = f.brep_volume(Vs);
            }
            cc.remove_boundary_cells_by_volume(vol);
        }
        auto cb = cc.cell_boundaries();
        cell_boundaries.clear();

        cell_boundaries.resize(cb.size());
        for (auto &&[a, b] : mtao::iterator::zip(cell_boundaries, cb)) {
            a.insert(b.begin(), b.end());
        }
    }
    if (adaptive) {
        auto t = mtao::logging::profiler("Adaptive grid", false, "profiler");
        assert(m_active_grid_cell_mask.shape() == cell_shape());
        auto adaptive_grid_factory = AdaptiveGridFactory(m_active_grid_cell_mask);
        adaptive_grid_factory.make_cells(adaptive_level);
        adaptive_grid = adaptive_grid_factory.create();
    }

    {
        auto &&cells = cell_boundaries;
        int max_cell_id = cells.size();
        if (adaptive) {
            //make cell names unique
            int cbsize = cells.size();
            auto &ag = *adaptive_grid;
            auto gc = ag.m_cells;
            ag.m_cells.clear();
            for (auto &&[i, c] : mtao::iterator::enumerate(gc)) {
                auto &&[j, b] = c;
                int id = i + cbsize;
                ag.m_cells[id] = b;
            }
            max_cell_id = cbsize + ag.m_cells.size();
        }
        mtao::data_structures::DisjointSet<int> cell_ds;
        {
            auto t = mtao::logging::profiler("region disjoint set construction", false, "profiler");
            for (int i = 0; i < cells.size(); ++i) {
                cell_ds.add_node(i);
            }
            for (auto &&[a, b] : m_faces) {
                if (a >= 0) {
                    cell_ds.add_node(max_cell_id + a);
                }
            }
            if (adaptive) {
                auto &ag = *adaptive_grid;
                auto grid = ag.cell_ownership_grid();
                auto ag_faces = ag.faces(grid);
                for (auto &&[c, b] : ag.cells()) {
                    cell_ds.add_node(c);
                }
                for (auto &&f : ag_faces) {
                    auto &&[a, b] = f.dual_edge;
                    if (a >= 0 && b >= 0) {
                        cell_ds.join(a, b);
                        //spdlog::info("ad{} <=> ad{}", a,b);
                    }
                }
                for (auto &&[id, f] : faces()) {
                    if (f.external_boundary) {
                        auto &[cid, s] = *f.external_boundary;
                        if (cid >= 0) {
                            cell_ds.join(max_cell_id + id, grid.get(cid));
                            //spdlog::info("c{} <=> ad{}", cid,grid.get(cid));
                        }
                    }
                }
            }

            for (auto [cid, faces] : mtao::iterator::enumerate(cells)) {
                for (auto &&[fid, s] : faces) {
                    if (fid >= 0) {
                        auto &&f = m_faces[fid];
                        if (!f.is_mesh_face()) {
                            cell_ds.join(cid, fid + max_cell_id);
                            //spdlog::info("c{} <=> f{}", cid,fid);
                        }
                    }
                }
            }

            cell_ds.reduce_all();
        }

        int min_face_idx = 0;
        {
            int min_face_x = vertex_shape()[0];
            for (auto &&[i, f] : m_faces) {
                if (f[0]) {
                    int x = *f[0];
                    if (x < min_face_x) {
                        min_face_idx = i;
                        min_face_x = x;
                    }
                }
            }
        }
        /*{
              auto minf = m_faces[min_face_idx];
              std::cout << "Min facE: ";
              for(auto&& f: minf.indices) {
              std::copy(f.begin(),f.end(),std::ostream_iterator<int>(std::cout,":"));
              std::cout << " ";
              }

              std::cout << std::endl;
              }
              */
        int outside_root = -1;
        for (auto [cid, faces] : mtao::iterator::enumerate(cells)) {
            for (auto &&[fid, s] : faces) {
                if (fid >= 0) {
                    if (fid == min_face_idx) {
                        outside_root = cell_ds.get_root(cid).data;
                        break;
                    }
                }
            }
            if (outside_root != -1) {
                break;
            }
        }
        std::map<int, int> reindexer;
        reindexer[outside_root] = 0;
        for (int i = 0; i < cells.size(); ++i) {
            int root = cell_ds.get_root(i).data;
            if (root != outside_root && reindexer.find(root) == reindexer.end()) {
                reindexer[root] = reindexer.size();
            }
        }


        std::vector<int> ret(cells.size());
        std::set<int> regions;
        for (int i = 0; i < cells.size(); ++i) {
            cells[i].index = i;
            cells[i].region = reindexer[cell_ds.get_root(i).data];
            regions.insert(cells[i].region);
        }
        if (adaptive) {

            auto &ag = *adaptive_grid;
            auto &agr = *(adaptive_grid_regions = std::map<int, int>());
            for (auto &&[cid, b] : ag.cells()) {
                agr[cid] = reindexer[cell_ds.get_root(cid).data];
            }
        }


        warn() << "Region count: " << regions.size();
        auto w = warn();
        /*
            for(auto&& r: regions) {
                w << r << "(" << std::count_if(cells.begin(),cells.end(), [&](auto&& c) {
                            return c.region == r;
                            })<< ") ";
            }
            for(auto&& r: regions) {
                std::cout << r << "==========================\n";
                for(auto&& c: cells) {
                    if(c.region == r) {
                        std::cout << c.size() << ":";
                    }
                }
                std::cout << std::endl;
                for(auto&& c: cells) {
                    if(c.region == r) {
                        for(auto&& [a,b]: c) {
                            std::cout << a << std::endl;
                            std::cout << std::string(m_faces.at(a)) << std::endl;
                        }
                        std::cout << std::endl;
                    }
                }
            }
            */
    }
    mtao::logging::debug() << "Cell count: " << cell_boundaries.size();
}
mtao::Vec3d CutCellGenerator<3>::area_normal(const std::vector<int> &F) const {
    auto V = all_GV();
    mtao::Vec3d N = mtao::Vec3d::Zero();
    for (int i = 0; i < F.size(); ++i) {
        int j = (i + 1) % F.size();
        int k = (i + 2) % F.size();
        auto x = V.col(F[i]);
        auto y = V.col(F[j]);
        auto z = V.col(F[k]);
        (z - y).cross(y - x);
        N += (z - y).cross(x - y);
    }
    return N / 2;
}
mtao::Vec3d CutCellGenerator<3>::area_normal(const std::set<std::vector<int>> &F) const {
    mtao::Vec3d N = mtao::Vec3d::Zero();
    for (auto &&f : F) {
        N += area_normal(f);
    }
    return N;
}

bool CutCellGenerator<3>::check_face_utilization() const {
    auto gv = all_GV();
    std::set<int> faces;
    for (auto &&[fidx, cs] : m_faces) {
        bool do_it = true;
        for (auto &&v : cs.indices) {

            for (int i = 0; i < 3; ++i) {
                if (cs[i]) {
                    if (axial_primal_faces[i].find(smallest_ordered_edge(v)) != axial_primal_faces[i].end()) {
                        do_it = false;
                    }
                }
            }
        }
        if (do_it) {
            faces.insert(fidx);
        }
    }
    for (auto &&cell : cell_boundaries) {
        for (auto &&[f, s] : cell) {
            faces.erase(f);
        }
    }
    for (auto &&f : faces) {
        spdlog::error("Unused face");
        std::cout << "Unused face: " << f << std::endl;
        auto face = m_faces.at(f);
        std::cout << std::string(face) << std::endl;
        std::cout << std::string(face.mask()) << std::endl;
        for (auto &&f : face.indices) {
            std::copy(f.begin(), f.end(), std::ostream_iterator<int>(std::cout, ":"));
            std::cout << " ";
        }
        std::cout << std::endl;
    }

    return faces.empty();
}

bool CutCellGenerator<3>::check_cell_containment() const {
    auto gv = all_GV();
    bool ret = true;
    for (auto &&[fidx, cs] : m_faces) {
        if (!is_in_cell(cs.indices)) {
            for (auto &&c : cs.indices) {
                std::cout << "Bad face, lacks containment " << std::endl;
                for (auto &&v : c) {
                    std::string e = std::string(GV(v));
                    std::cout << v << ":";
                }
                std::cout << std::endl;
                for (auto &&v : c) {
                    std::string e = std::string(GV(v));
                    std::cout << v << "(" << gv.col(v).transpose() << ")";
                }
                std::cout << std::endl;
                ret = false;
            }
        }
    }
    return ret;
}
auto CutCellGenerator<3>::edge_slice(int dim, int coord) const -> std::set<Edge> {
    if (auto it = axis_hem_data[dim].find(coord); it != axis_hem_data[dim].end()) {
        return it->second.edges;
    } else {
        return {};
    }
}
}// namespace mandoline::construction
