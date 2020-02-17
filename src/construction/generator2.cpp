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
        mtao::logging::debug() << "Number of faces: " << hem.cell_halfedges().size();
    }
    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate_faces() const {
            auto ret = generate_edges();


            std::tie(ret.hem,std::ignore) = compute_planar_hem(all_V(),mtao::eigen::stl2eigen(edges()),m_active_grid_cell_mask,StaggeredGrid::cell_size());


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



}
