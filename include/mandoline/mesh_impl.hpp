#pragma once
#include "mandoline/mesh.hpp"
#include "mtao/geometry/grid/triangulation.hpp"

namespace mandoline {
    template <int D, typename Derived>
            auto CutCellMeshBase<D,Derived>::get_world_vertex(const Vertex& p) const -> Vec {
                auto&& g = vertex_grid();
                return g.world_coord(p.p());
                return g.origin() + (p.p().cwiseProduct(g.dx()));

            }

    template <int D, typename Derived>
        int CutCellMeshBase<D,Derived>::cut_vertex_size() const {
            return m_cut_vertices.size();
        }
    template <int D, typename Derived>
        int CutCellMeshBase<D,Derived>::vertex_size() const {
            return StaggeredGrid::vertex_size() + cut_vertex_size();
        }
    template <int D, typename Derived>
        int CutCellMeshBase<D,Derived>::edge_size() const {
            return StaggeredGrid::edge_size() + cut_edge_size();

        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::vertex(int idx) const -> Vec {
            return get_world_vertex(masked_vertex(idx));
        }

    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::masked_vertex(int idx) const -> Vertex {
            assert(idx >= 0);
            if(is_grid_vertex(idx) ) {
                return StaggeredGrid::vertex_unindex(idx);
            } else {
                assert(idx - StaggeredGrid::vertex_size() < cut_vertex_size());
                return cut_vertex(idx - StaggeredGrid::vertex_size());
            }
        }

    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::grid_boundary(bool dirichlet_boundary) const -> SpMat {
            auto B = SpMat{};
            //auto B = StaggeredGrid::boundary<D>();
            //TODO: implement dirichlet boundary
            return B;
        }

    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::edge(int idx) const -> Edge{

            if(is_grid_edge(idx)) {
                return grid_edge(idx);
            } else {
                Edge e;
                IVecMap(e.data()) = cut_edge(idx - StaggeredGrid::edge_size());
                return e;
            }
        }

    template <int D, typename Derived>
        int CutCellMeshBase<D,Derived>::grid_edge_type(int idx) const {
            return StaggeredGrid::form_type<1>(idx);
        }

    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::edge_volumes() const -> VecX{
            VecX V(edge_size());
            auto offsets = StaggeredGrid::template offsets<1>();
            for(int i = 0; i < D; ++i) {
                V.segment(offsets[i],StaggeredGrid::template staggered_size<1>(i)).array() = dx()(i);
            }
            auto ceV = V.bottomRows(cut_edge_size());
            for(int i = 0; i < cut_edge_size(); ++i) {
                auto e = cut_edge(i);
                ceV(i) = (vertex(e(1)) - vertex(e(0))).norm();
            }
            return V;
        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::vertices() const -> ColVecs {
            return mtao::eigen::hstack(StaggeredGrid::vertices(),cut_vertices_colvecs());
        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::grid_space_vertices() const -> ColVecs {
            return mtao::eigen::hstack(StaggeredGrid::local_vertices(),cut_vertices_colvecs());
        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::grid_space_cut_vertices_colvecs() const -> ColVecs {
            ColVecs V(3,cut_vertex_size());
            for(auto&& [i,v]: mtao::iterator::enumerate(cut_vertices())) {
                V.col(i) = v.p();
            }
            return V;

        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::cut_vertices_colvecs() const -> ColVecs {
            auto&& g = vertex_grid();
            return g.world_coord(grid_space_cut_vertices_colvecs());
            /*
            ColVecs V(3,cut_vertex_size());
            for(auto&& [i,v]: mtao::iterator::enumerate(cut_vertices())) {
                V.col(i) = get_world_vertex(v);
            }
            return V;
            */

        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::grid_edge(coord_type edge_coord, int type) const -> Edge {
            int fidx = StaggeredGrid::vertex_index(edge_coord);
            edge_coord[type]++;
            int nidx = StaggeredGrid::vertex_index(edge_coord);
            return Edge{fidx,nidx};
        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::grid_edge(int idx) const -> Edge {
            int et = StaggeredGrid::edge_type(idx);
            auto coord = StaggeredGrid::template staggered_unindex<1>(idx,et);
            return grid_edge(coord,et);
        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::grid_dual_edge(int idx) const -> Edge {
            int et = StaggeredGrid::edge_type(idx);
            auto coord = StaggeredGrid::template staggered_unindex<D-1>(idx,et);
            int fidx = StaggeredGrid::template staggered_index<D>(coord,et);
            coord[et]--;
            int nidx = StaggeredGrid::template staggered_index<D>(coord,et);
            return Edge{{et,nidx}};
        }
    template <int D, typename Derived>
        auto CutCellMeshBase<D,Derived>::edges() const -> Edges {
            return mtao::eigen::hstack(mtao::geometry::grid::GridTriangulator<StaggeredGrid>(vertex_grid()).edges(),cut_edges());
        }
}
