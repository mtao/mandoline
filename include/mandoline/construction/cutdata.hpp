#pragma once
#ifdef _OPENMP
# include <omp.h>
#endif
#include <mtao/types.hpp>
#include <mtao/geometry/mesh/boundary_elements.h>
#include <mtao/geometry/mesh/boundary_facets.h>
#include "mandoline/construction/facet_intersections.hpp"
#include <mtao/eigen/stack.h>
#include <mtao/logging/timer.hpp>
#include <mtao/logging/profiler.hpp>


namespace mandoline::construction {
    template <int D>
        class CutCellEdgeGenerator;
    template <int D, typename Indexer_  = typename EdgeIntersections<D>::SGType::Indexer>
        struct CutData: public Indexer_ {
            friend class CutCellEdgeGenerator<D>;
            using Edge = std::array<int,2>;
            using Edges = mtao::ColVectors<int,2>;
            using Vec = mtao::Vec3d;
            using Faces = mtao::ColVectors<int,3>;
            using Triangles = mtao::ColVectors<int,3>;
            using ColVecs = mtao::ColVectors<double,D>;
            using EdgeIsect = EdgeIntersection<D>;
            using TriIsect = TriangleIntersection<D>;
            using EdgeIsects = EdgeIntersections<D>;
            using TriIsects = TriangleIntersections<D>;
            using SGType = typename EdgeIntersections<D>::SGType;
            using VType = Vertex<D>;
            using VPtrEdge = std::array<const VType*,2>;

            using Indexer = Indexer_;

            CutData() = default;
            CutData(CutData&&) = default;
            CutData& operator=(CutData&&) = default;
            CutData(const Indexer& indexer): Indexer(indexer) {}
            CutData(const Indexer& indexer, const mtao::vector<VType>& V, const Faces& F = {});
            CutData(const Indexer& indexer, const mtao::vector<VType>& V, const Edges& E,const Faces& F = {},  const Faces& FEA = {});

            void set_topology(const Faces& F = {});
            void set_topology(const Edges& E,const Faces& F = {},  const Faces& FEA = {}) ;

            void bake(const std::optional<SGType>& grid= {}, bool fuse=true);
            void clear();
            void reset_intersections() ;





            const std::vector<Crossing<D>>& crossings() const;
            std::set<Edge> stl_edges() const ;
            mtao::ColVectors<int,2> edges() const ;
            std::vector<std::vector<int>> faces() const;
            const std::vector<CutMeshFace<D>>& cut_faces() const { return m_cut_faces; }


            void update_vertices(const mtao::vector<VType>& V);
            void set_vertices(const mtao::vector<VType>& V) { m_V = V; }



            void add_edge_isect(int idx, double t) ;


            ColVecs cut_vertices() const;

            //this includes the grid vertices
            int cut_vertex_size() const;
            int nV() const { return m_V.size(); }
            int nE() const { return m_E.cols(); }
            int nF() const { return m_F.cols(); }
            int nFE() const { return m_FE.cols(); }
            auto V(int idx) const { return m_V[idx]; }
            auto E(int idx) const { return m_E.col(idx); }
            auto F(int idx) const { return m_F.col(idx); }
            auto FE(int idx) const { return m_FE.col(idx); }
            const auto& V() const { return m_V; }
            const auto& E() const { return m_E; }
            const auto& F() const { return m_F; }
            const auto& FE() const { return m_FE; }
            Faces F_asCrossings() const;

            std::map<std::vector<int>,int> triangle_index_ownership() const;
            std::map<Edge,int> edge_index_ownership() const;
            mtao::Vec3d get_bary(int face_index, const Crossing<D>& crossing) const;


            int index(const VType* v) const { return m_vertex_indexer.at(v); }
            int grid_size() const { return Indexer::size(); }
            int grid_index(const VType& v) const { return Indexer::index(v.coord); }
            int grid_index(const VType* v) const { return grid_index(*v); }


            Eigen::SparseMatrix<double> barycentric_map() const;

            //TODO: figure out what i wanted from these
            void clean_edges();
            void clean_triangles();
            private:
            std::vector<const EdgeIntersection<D>*> flat_edge_intersections() const ;
            std::vector<const TriangleIntersection<D>*> flat_triangle_intersections() const ;
            std::vector<Crossing<D>> compute_crossings(bool fuse=true) const ;
            std::vector<Crossing<D>> vertex_crossings() const ;
            std::vector<Crossing<D>> edge_crossings() const ;
            std::vector<Crossing<D>> face_crossings() const ;

            template <typename IsectType>
                std::vector<Crossing<D>> facet_crossings(const std::vector<const IsectType*>& isects) const ;
            std::map<const VType*, int> gv_index_map(const std::vector<Crossing<D>>& crossings) const ;

            template <int U>
                mtao::ColVectors<int,U> facets(const std::map<const VType*,int>& gv_idx_map, const std::vector<std::array<const VType*,U>>& facets) const ;

            std::vector<Edge> gvedge2edges(const std::vector<VPtrEdge>& gvedges, std::map<const VType*,int>& indexer) const ;
            mtao::ColVectors<int,2> edge_edges() const ;
            std::map<VPtrEdge,int> edge_indices() const ;
            std::set<Edge> edge_edges(const std::map<const VType*,int>& gv_idx_map) const ;

            std::map<VPtrEdge, std::set<VPtrEdge>> edge_ownership() const ;

            std::set<Edge> face_edges() const ;
            std::set<Edge> face_edges(const std::map<const VType*,int>& gv_idx_map) const ;

            mtao::ColVectors<int,2> edges(const std::map<const VType*,int>& gv_idx_map) const ;
            std::set<Edge> stl_edges(const std::map<const VType*,int>& gv_idx_map) const ;

            private:


            mtao::vector<Vertex<D>> m_V;
            Edges m_E;
            Faces m_F;
            Faces m_FE;
            std::vector<EdgeIntersections<D>> m_edge_intersections;
            std::vector<TriangleIntersections<D>> m_triangle_intersections;
            //needs to be baked
            std::vector<Crossing<D>> m_crossings;
            std::map<const VType*,int> m_vertex_indexer;
            mtao::vector<CutMeshFace<D>> m_cut_faces;

        };


}
#include "mandoline/construction/cutdata_impl.hpp"
