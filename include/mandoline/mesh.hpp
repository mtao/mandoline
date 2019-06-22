#pragma once
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop
#include <mtao/geometry/mesh/halfedge.hpp>
#include <array>
#include <map>
#include <vector>
#include "mandoline/cell.hpp"
#include "cutmesh.pb.h"
#include "mandoline/adaptive_grid.hpp"
#include <set>
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/eigen/stack.h>
#include "mandoline/barycentric_triangle_face.hpp"
#include "mandoline/vertex.hpp"

namespace mandoline {
    namespace construction {
    template <int D>
        class CutCellEdgeGenerator;
    template <int D>
        class CutCellGenerator;
    }

    template <int D, typename Derived_>
        struct CutCellMeshBase: public mtao::geometry::grid::StaggeredGrid<double,D> {
            friend class construction::CutCellEdgeGenerator<D>;
            friend class construction::CutCellGenerator<D>;
            using Derived = Derived_;
            Derived& derived() { return *static_cast<Derived*>(this); }
            const Derived& derived() const { return *static_cast<const Derived*>(this); }
            using StaggeredGrid = mtao::geometry::grid::StaggeredGrid<double,D>;
            using StaggeredGrid::vertex_grid;
            using coord_type = typename StaggeredGrid::coord_type;
            using GridData = mtao::geometry::grid::GridDataD<double,D>;
            using GridDatab = mtao::geometry::grid::GridDataD<bool,D>;
            using ColVecs = mtao::ColVectors<double,D>;
            using VecX = mtao::VectorX<double>;
            using Edges = mtao::ColVectors<int,2>;
            using Edge = std::array<int,2>;
            using Vec = mtao::Vector<double,D>;
            using VecRef = Eigen::Ref<Vec>;
            using VecCRef = Eigen::Ref<const Vec>;
            using VecMap = Eigen::Map<Vec>;
            using VecCMap = Eigen::Map<const Vec>;
            using IVec = mtao::Vector<int,D>;
            using IVecRef = Eigen::Ref<IVec>;
            using IVecCRef = Eigen::Ref<const IVec>;
            using IVecMap = Eigen::Map<IVec>;
            using IVecCMap = Eigen::Map<const IVec>;
            using SpMat = Eigen::SparseMatrix<double>;
            using Vertex = ::mandoline::Vertex<D>;
            CutCellMeshBase() = default;
            CutCellMeshBase(const CutCellMeshBase&) = default;
            CutCellMeshBase(CutCellMeshBase&&) = default;
            CutCellMeshBase& operator=(const CutCellMeshBase&) = default;
            CutCellMeshBase& operator=(CutCellMeshBase&&) = default;
            CutCellMeshBase(const StaggeredGrid& g, const std::vector<Vertex>& v = {}): StaggeredGrid(g), m_cut_vertices(v), m_active_grid_cell_mask(GridDatab::Constant(true,g.cell_shape())) {
            }
            using StaggeredGrid::dx;

            ColVecs vertices() const;
            ColVecs dual_vertices() const;
            ColVecs cut_vertices_colvecs() const;

            Vec vertex(int vertex) const;
            Vertex masked_vertex(int vertex) const;//in grid space, with masking
            const Vertex& cut_vertex(int vertex) const {
                return m_cut_vertices[vertex];
            }
            //Vec dual_vertex(int vertex) const;
            int vertex_size() const;
            int cut_vertex_size() const;
            int edge_size() const;


            VecX volumes() const { return derived().volumes(); }
            VecX dual_edge_volumes() const { return derived().dual_edge_volumes(); }
            SpMat grid_boundary(bool dirichlet_boundary) const;
            VecX edge_volumes() const;

            Edge grid_edge(int idx) const;
            int grid_edge_type(int idx) const;
            Edge grid_edge(coord_type edge_coord, int type) const;
            Edge grid_dual_edge(int idx) const;
            Edge edge(int idx) const;
            Edge dual_edge(int idx) const;
            Edges edges() const;

            auto cut_edge(int idx) const { return m_cut_edges.col(idx); }
            int cut_edge_size() const { return m_cut_edges.cols(); }

            auto grid_vertices() const { return StaggeredGrid::vertices(); }
            auto grid_vertex(int idx) const { return StaggeredGrid::vertex(idx); }
            auto vertex(const coord_type& coord) const { return StaggeredGrid::vertex(coord); }
            bool is_grid_vertex(int index) const { return index < StaggeredGrid::template form_size<0>(); }
            bool is_grid_edge(int index) const { return index < StaggeredGrid::template form_size<1>(); }
            bool is_grid_face(int index) const { return index < StaggeredGrid::template form_size<2>(); }

            int cell_index(const VecCRef& v) const { return derived().cell_index(v); }

            Eigen::SparseMatrix<double> boundary(bool dirichlet_boundary) const;


            Vec get_world_vertex(const Vertex& p) const;

            const std::vector<Vertex>& cut_vertices() const {return m_cut_vertices;}

            const GridDatab& active_grid_cell_mask() const { return m_active_grid_cell_mask;}
            const Edges& cut_edges() const {return m_cut_edges;}

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            protected:
            std::vector<Vertex> m_cut_vertices;

            GridDatab m_active_grid_cell_mask;
            Edges m_cut_edges;
        };
    template <int D>
        struct CutCellMesh;
}

#include "mandoline/mesh_impl.hpp"


