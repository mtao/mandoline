#include "mandoline/construction/generator.hpp"
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/colvector_loop.hpp>

namespace mandoline::construction {
    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate() const {
            return generate_edges();
        }


    CutCellMesh<2> CutCellGenerator<2>::generate() const {
        auto ccm = generate_faces();
        return ccm;
    }


    void CutCellGenerator<2>::add_boundary_elements(const BoundaryElements& E) {
        //BoundaryElements
        add_edges(E);
    }


    void CutCellGenerator<2>::bake_faces() {

        std::tie(hem,std::ignore) = compute_planar_hem(all_GV(),mtao::eigen::stl2eigen(edges()),m_active_grid_cell_mask,StaggeredGrid::cell_size());

        std::vector<bool> active(data().nE());// edges that sit on an axis
        auto&& ti = data().edge_intersections();
        std::transform(ti.begin(),ti.end(), active.begin(), [](auto&& ti) {
                return ti.active();
                });
        cut_cell_to_primal_map.clear();



        mtao::logging::debug() << "Number of faces: " << hem.cell_halfedges().size();

        for(auto&& [i,f]: mtao::iterator::enumerate(m_cut_faces)) {
            cut_cell_to_primal_map[i] = f.parent_fid;
            if(active[f.parent_fid]) {
                axial_primal_faces[i].insert(smallest_ordered_edge(f.indices));
            }
        }
    }
    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate_faces() const {
            // generate the edge structure
            auto ret = generate_edges();


            // make a halfedge mesh and store it
            std::tie(ret.hem,std::ignore) = compute_planar_hem(all_V(),mtao::eigen::stl2eigen(edges()),m_active_grid_cell_mask,0);




            {
                //auto t = mtao::logging::timer("Gluing halfedges to edges");
                mtao::map<std::array<int,2>,int> emap;
                for(int i = 0; i < ret.m_cut_edges.cols(); ++i) {
                    std::array<int,2> v;
                    IVecMap(v.data()) = ret.m_cut_edges.col(i);
                    emap[v] = i;
                    std::swap(v[0],v[1]);
                    emap[v] = i;
                }
                ret.halfedges_per_edge.resizeLike(ret.m_cut_edges);


                for(int i = 0; i < ret.hem.size(); ++i) {
                    std::array<int,2> v;
                    auto e = ret.hem.edge(i);
                    v[0] = e.vertex();
                    e.dual();
                    int dhe = int(e);
                    v[1] = e.vertex();
                    int idx = ret.halfedge_to_edge_index[i] = emap[v];
                    auto hpee = ret.halfedges_per_edge.col(idx);
                    auto ee = ret.m_cut_edges.col(idx);
                    if(v[0] == ee(0)) {
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
                //auto t = mtao::logging::timer("Cell grid ownership");

                auto dv = ret.dual_vertices();
                auto  cells = ret.hem.cell_halfedges();
                for(auto&& c: cells) {
                    int cell_ind = ret.hem.cell_index(c);
                    if(is_valid_grid_index(cell_ind)) {
                        continue;
                    }
                    if((cell_ind) != 0) {
                        auto ij = get_grid_cell(dv.col(cell_ind));
                        int grid_ind = StaggeredGrid::cell_index(ij);
                        ret.cell_grid_ownership[grid_ind].insert(cell_ind);
                    }
                }
            }
            //extra_metadata(ret);
            return ret;
        }


    void CutCellGenerator<2>::extra_metadata(CutCellMesh<2>& mesh) const {

    }


    void CutCellGenerator<2>::compute_faces_vertex() {
        auto t = mtao::logging::timer("Computing input mesh faces with normals ");
        std::array<double,3> N;
        for(auto&& [i,f]: mtao::iterator::enumerate(m_cut_faces)) {
            if(f.indices.size() > 2) {// && f.mask().count() < 2) {

                m_faces[i] = CutFace<D>(f,origN.col(f.parent_fid));
                mesh_face_indices.insert(i);
            } else {
                std::cout << "Cut faces should always be  real faces!" << std::endl;
            }
            }
        }

        mtao::map<int,CutFace<3>> CutCellGenerator<3>::compute_faces_axis()const {

            mtao::map<int,CutFace<D>> faces;

            auto V = all_GV();


            auto mesh_faces = hem.cells_multi_component_map();

            auto vols = hem.signed_areas(V);
            for(auto&& [i,v]: mesh_faces) {
                CutFace<D> F;
                F.id = Edge{{idx,cidx}};
                auto add_loop = [&](const std::vector<int>& v) {
                    if(v.size() > 2) {
                        if(axial_primal_faces[idx].find(smallest_ordered_edge(v)) == axial_primal_faces[idx].end()) {
                            F.indices.emplace(v);
                        }
                    }
                };
                // TODO:remove handling this crap
                /*
                for(auto& v: v) {
                    if(ahd.is_boundary_cell(v)) {
                        auto pc = possible_cells(v);//possible grid cell
                        if(pc.empty()) {
                            continue;
                        }
                        std::array<CoordType,2> pca;
                        std::copy(pc.begin(),pc.end(),pca.begin());
                        if(pca[0][idx]+1 != pca[1][idx]) {
                            std::cout << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?" << std::endl;
                        }
                        std::array<int,2> ind;
                        ind[0] = pca[0][(idx+1)%3];
                        ind[1] = pca[0][(idx+2)%3];
                        if(!ahd.active_grid_cell_mask.valid_index(ind)) {
                            std::cout << "Invalid index???? " << std::endl;
                        }
                        if(!ahd.active_grid_cell_mask(ind)) {
                            if(vols.at(i) >= 0) {
                                add_loop(v);
                            }
                        }
                    } else {
                        if(vols.at(i) >= 0) {
                            add_loop(v);
                        }
                    }
                }
                */
                if(!F.indices.empty()) {

                    F.N = N;
                    F.coord_mask<D>::operator=(face_mask(F.indices));
                    faces[i] = std::move(F);
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

    }
