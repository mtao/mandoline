#pragma once
#ifdef MTAO_OPENMP
#include <omp.h>
#endif
#include <balsa/eigen/types.hpp>
#include <mtao/geometry/mesh/boundary_elements.h>
#include <mtao/geometry/mesh/boundary_facets.h>
#include "mandoline/construction/facet_intersections.hpp"
#include <balsa/eigen/stack.hpp>
#include <mtao/logging/timer.hpp>
#include <mtao/logging/profiler.hpp>


namespace mandoline::construction {
template<int D>
class CutCellEdgeGenerator;
template<int D, typename Indexer_ = typename EdgeIntersections<D>::SGType::GridType::Indexer>
struct CutData : public Indexer_ {
    friend class CutCellEdgeGenerator<D>;
    using Edge = std::array<int, 2>;
    using Edges = balsa::eigen::ColVectors<int, 2>;
    using Vec = balsa::eigen::Vec3d;
    using Faces = balsa::eigen::ColVectors<int, 3>;
    using Triangles = balsa::eigen::ColVectors<int, 3>;
    using ColVecs = balsa::eigen::ColVectors<double, D>;
    using EdgeIsect = EdgeIntersection<D>;
    using TriIsect = TriangleIntersection<D>;
    using EdgeIsects = EdgeIntersections<D>;
    using TriIsects = TriangleIntersections<D>;
    using SGType = typename EdgeIntersections<D>::SGType;
    using GType = typename SGType::GridType;
    using VType = Vertex<D>;
    using VPtrEdge = std::array<const VType *, 2>;

    using Indexer = Indexer_;

    CutData() = default;
    CutData(CutData &&) = default;
    CutData &operator=(CutData &&) = default;
    CutData(const Indexer &indexer) : Indexer(indexer) {}
    CutData(const Indexer &indexer, const mtao::vector<VType> &V, const Faces &F = {});
    CutData(const Indexer &indexer, const mtao::vector<VType> &V, const Edges &E, const Faces &F = {}, const Faces &FEA = {});

    void set_topology(const Faces &F = {});
    void set_topology(const Edges &E, const Faces &F = {}, const Faces &FEA = {});
    void update_topology_masks();

    void bake(const std::optional<SGType> &grid = {}, bool fuse = true, bool filter_external = false);
    void clear();// clear intersections, useful if vertices changed position
    void reset();// reset internal data
    void reset_topology();// reset internal data
    void reset_intersections();


    const std::vector<Crossing<D>> &crossings() const;
    std::set<Edge> stl_edges() const;
    balsa::eigen::ColVectors<int, 2> edges() const;
    std::vector<std::vector<int>> faces() const;
    const std::vector<CutMeshFace<D>> &cut_faces() const { return m_cut_faces; }
    const std::vector<CutMeshEdge<D>> &cut_edges() const { return m_cut_edges; }


    void update_vertices(const mtao::vector<VType> &V);
    void update_grid(const Indexer &indexer);
    template<typename S>
    void update_grid(const mtao::geometry::grid::StaggeredGrid<S, D> &sg) { update_grid(sg.vertex_grid()); }


    void set_vertices(const mtao::vector<VType> &V) { m_V = V; }
    void set_vertices(mtao::vector<VType> &&V) { m_V = std::move(V); }


    void add_edge_isect(int idx, double t);


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
    const auto &V() const { return m_V; }
    const auto &E() const { return m_E; }
    const auto &F() const { return m_F; }
    const auto &FE() const { return m_FE; }
    Faces F_asCrossings() const;

    std::map<std::vector<int>, int> triangle_index_ownership() const;
    std::map<Edge, int> edge_index_ownership() const;
    // get the barycentric coordinates of a crossing along an edge
    balsa::eigen::Vec3d get_face_bary(int face_index, const Crossing<D> &crossing) const;
    // get the coordinate over an edge
    double get_edge_coord(int edge_index, const Crossing<D> &crossing) const;


    int index(const VType *v) const { return m_vertex_indexer.at(v); }
    int grid_size() const { return Indexer::size(); }
    int grid_index(const VType &v) const { return Indexer::index(v.coord); }
    int grid_index(const VType *v) const { return grid_index(*v); }


    Eigen::SparseMatrix<double> barycentric_map() const;

    //TODO: figure out what i wanted from these
    void clean_edges();
    void clean_triangles();

    const std::vector<EdgeIntersections<D>> &edge_intersections() const { return m_edge_intersections; }
    const std::vector<TriangleIntersections<D>> &triangle_intersections() const { return m_triangle_intersections; }

  private:
    std::vector<const EdgeIntersection<D> *> flat_edge_intersections() const;
    std::vector<const TriangleIntersection<D> *> flat_triangle_intersections() const;
    /// fuse: identifies crossings that are actually grid vertices and fuses them
    // this should pretty much always be run, but some inputs can be really bad FP-wise anyway
    std::vector<Crossing<D>> compute_crossings(bool fuse = true, bool filter_external = false) const;
    std::vector<Crossing<D>> vertex_crossings() const;
    std::vector<Crossing<D>> edge_crossings() const;
    std::vector<Crossing<D>> face_crossings() const;

    template<typename IsectType>
    std::vector<Crossing<D>> facet_crossings(const std::vector<const IsectType *> &isects) const;
    std::map<const VType *, int> gv_index_map(const std::vector<Crossing<D>> &crossings) const;

    template<int U>
    balsa::eigen::ColVectors<int, U> facets(const std::map<const VType *, int> &gv_idx_map, const std::vector<std::array<const VType *, U>> &facets) const;

    std::vector<Edge> gvedge2edges(const std::vector<VPtrEdge> &gvedges, std::map<const VType *, int> &indexer) const;
    balsa::eigen::ColVectors<int, 2> edge_edges() const;
    std::map<VPtrEdge, int> edge_indices() const;
    std::set<Edge> edge_edges(const std::map<const VType *, int> &gv_idx_map) const;

    std::map<VPtrEdge, std::set<VPtrEdge>> edge_ownership() const;

    std::set<Edge> face_edges() const;
    std::set<Edge> face_edges(const std::map<const VType *, int> &gv_idx_map) const;

    balsa::eigen::ColVectors<int, 2> edges(const std::map<const VType *, int> &gv_idx_map) const;
    std::set<Edge> stl_edges(const std::map<const VType *, int> &gv_idx_map) const;

  private:
    mtao::vector<Vertex<D>> m_V;
    Edges m_E;
    Faces m_F;
    Faces m_FE;
    std::vector<EdgeIntersections<D>> m_edge_intersections;
    std::vector<TriangleIntersections<D>> m_triangle_intersections;
    //needs to be baked
    std::vector<Crossing<D>> m_crossings;
    std::map<const VType *, int> m_vertex_indexer;
    mtao::vector<CutMeshEdge<D>> m_cut_edges;
    mtao::vector<CutMeshFace<D>> m_cut_faces;
};


}// namespace mandoline::construction
#include "mandoline/construction/cutdata_impl.hpp"
