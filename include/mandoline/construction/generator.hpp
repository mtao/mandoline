#pragma once
#include <vector>
#include <tuple>
#include <mtao/types.h>
#include <bitset>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/functional.hpp>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/geometry/grid/grid_data.hpp>
#include <mtao/logging/timer.hpp>
#include <map>
#include <set>
#include "mandoline/construction/cutdata.hpp"
#include "mandoline/cutface.hpp"
#include <iterator>
#include <mtao/geometry/mesh/halfedge.hpp>
#include <mtao/type_utils.h>

#define USE_INTERIOR_MASK

namespace mandoline {
template<int D>
class CutCellMesh;

}

namespace mandoline::construction {

//template<int D, typename IndexContainerType>
//std::optional<std::tuple<int, bool>> make_boundary_pair(const std::vector<Vertex<D>> &V, const CoordMaskedGeometry<D, IndexContainerType> &boundary_facet);

std::array<int, 2> smallest_ordered_edge(const std::vector<int> &v);
std::array<int, 2> smallest_ordered_edge_reverse(const std::vector<int> &v);
//ASSUMES SIMPLICIAL INPUTS
template<int D>
class CutCellEdgeGenerator : public mtao::geometry::grid::StaggeredGrid<double, D> {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    using Vec = mtao::Vector<double, D>;
    using Veci = mtao::Vector<int, D>;
    using ColVecs = mtao::ColVectors<double, D>;
    using Edge = std::array<int, 2>;
    using Edges = mtao::ColVectors<int, 2>;
    using BoundaryElements = mtao::ColVectors<int, D>;
    using VecVector = mtao::vector<Vec>;
    using StaggeredGrid = mtao::geometry::grid::StaggeredGrid<double, D>;
    using GridType = mtao::geometry::grid::GridD<double, D>;
    using StaggeredGrid::grid;
    using StaggeredGrid::vertex_grid;
    using StaggeredGrid::vertex_shape;
    using StaggeredGrid::vertex_unindex;
    using StaggeredGrid::cell_grid;
    using StaggeredGrid::origin;
    using StaggeredGrid::dx;
    using StaggeredGrid::vertex_index;
    using StaggeredGrid::cell_index;
    using StaggeredGrid::vertex;
    using coord_type = std::array<int, D>;
    //vertex that lives in the grid as a grid cell + local offset
    using VType = Vertex<D>;
    using EdgeIntersectionType = EdgeIntersection<D>;
    using CrossingType = Crossing<D>;
    using CoordMaskedEdgeType = CoordMaskedEdge<D>;
    using GridDatab = mtao::geometry::grid::GridDataD<bool, D>;



    // GLOBAL OPTIONS THAT CAN BE SET BETWEEN SETTING AND BAKING!
    bool toss_external_facets = false;

    //returns a new mesh and a set of boundary vertices
    std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> compute_planar_hem(const std::vector<VType> &GV, const ColVecs &V, const Edges &E, const GridDatab &interior_cell_mask) const;
    // auxilliary call for if we haven't created t
    std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> compute_planar_hem(const std::vector<VType> &GV, const Edges &E, const GridDatab &interior_cell_mask) const;
    // internal call that has no awareness of grid vertices / grid pruning
    std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> compute_planar_hem(const ColVecs &V, const Edges &E, const GridDatab &interior_cell_mask) const;
    std::tuple<mtao::geometry::mesh::HalfEdgeMesh, std::set<Edge>> compute_planar_hem(const std::map<std::array<int, 2>, std::tuple<int, bool>> &tangent_map, const ColVecs &T, const Edges &E, const GridDatab &interior_cell_mask) const;


    virtual void bake();
    virtual void bake_vertices();
    void bake_active_grid_cell_mask();
    virtual void bake_edges();
    virtual void bake_faces() {}
    virtual void bake_cells() {}


    virtual void clear();

    CutCellMesh<D> generate_vertices() const;// generate the initial mesh object
    CutCellMesh<D> generate_edges() const;// generate the edges on top of vertices
    CutCellMesh<D> generate_faces() const;// generate the faces on top of edges
    CutCellMesh<D> generate() const;

    //This threshold is to fuse vertices near grid vertices. for the results of float to double conversions 1e-6 seemed reasonable
    CutCellEdgeGenerator(const VecVector &V, const StaggeredGrid &grid, std::optional<double> threshold = 1e-6);
    CutCellEdgeGenerator(const ColVecs &V, const StaggeredGrid &grid, std::optional<double> threshold = 1e-6);
    CutCellEdgeGenerator(const StaggeredGrid &grid);
    CutCellEdgeGenerator() = default;
    CutCellEdgeGenerator(CutCellEdgeGenerator &&) = default;
    CutCellEdgeGenerator &operator=(CutCellEdgeGenerator &&) = default;
    //CutCellGenerator(const VecVector& V, const Vec& dx = Vec::Ones());

    void filter_external_cells(bool value = true) { m_filter_external_cells = value; }

    template<typename Derived>
    VType get_vertex(const Eigen::MatrixBase<Derived> &p) const {
        auto &&g = vertex_grid();
        Vec local = p.cwiseQuotient(g.dx()) - g.origin().cwiseQuotient(g.dx());
        return VType::from_vertex(local);
    }
    Vec get_world_vertex(const VType &p) const {
        auto &&g = vertex_grid();
        return g.origin() + (p.p().cwiseProduct(g.dx()));
    }

    size_t grid_vertex_size() const { return StaggeredGrid::vertex_size(); }

    bool is_grid_vertex(int idx) const { return idx < grid_vertex_size(); }
    bool is_orig_vertex(int idx) const { return !is_grid_vertex(idx) && idx < grid_vertex_size() + origV().size(); }
    bool is_new_vertex(int idx) const { return !is_orig_vertex(idx) && idx < num_vertices(); }


    void reset(const mtao::vector<VType> &grid_vertices);

    void add_boundary_elements(const BoundaryElements &E);
    void set_boundary_elements(const BoundaryElements &E);
    void add_edges(const Edges &E);
    static std::set<EdgeIntersectionType> cell_edge_intersections(const std::set<coord_type> &cells);
    template<typename Func>
    static void per_cell_vertex_looper(Func &&f, const coord_type &c) {
        for (int i = 0; i < (2 << D); ++i) {
            coord_type cc = c;
            std::bitset<D> bs(i);
            for (int j = 0; j < D; ++j) {
                cc[j] += bs[j] ? 1 : 0;
            }
            f(cc, bs);
        }
    }
    template<typename Func>
    static void per_dual_cell_vertex_looper(Func &&f, const coord_type &c) {
        for (int i = 0; i < (2 << D); ++i) {
            coord_type cc = c;
            std::bitset<D> bs(i);
            for (int j = 0; j < D; ++j) {
                cc[j] -= bs[j] ? 1 : 0;
            }
            f(cc, bs);
        }
    }
    template<typename Func>
    static void per_boundary_cell_vertex_looper(int N, Func &&f, const coord_type &c) {
        for (int i = 0; i < (2 << (D - 1)); ++i) {
            std::bitset<D> bs(i);
            if (bs[N]) {//if bit is 1 we move on
                continue;
            }
            coord_type cc = c;
            Vec v = Vec::Zero();
            for (int j = 0; j < D; ++j) {
                cc[j] += bs[j];
            }
            f(cc, bs);
        }
    }


    coord_mask<D> get_mask(int idx) const {
        int gvs = grid_vertex_size();
        using CM = coord_mask<D>;
        if (idx < gvs) {
            return CM{ StaggeredGrid::template staggered_unindex<0, 0>(idx) };
        } else {
            return crossing(idx).mask();
        }
    }

    // returns crossing according to the world index size. 
    auto &&crossing(int idx) const { return m_crossings[idx - grid_vertex_size()]; }

    size_t total_vertex_size() const { return m_crossings.size() + grid_vertex_size(); }
    auto &&crossings() const { return m_crossings; }

    size_t new_vertex_offset() const { return grid_vertex_size() + origV().size(); }

    std::set<coord_type> active_cells() const;
    std::array<mtao::map<coord_type, std::set<int>>, D> crossing_indices() const;

    auto &&origE() const { return data().E(); }
    auto origE(int idx) const { return data().E(idx); }


    const auto &origV() const { return m_origV; }
    const Vec &origV(int i) const { return m_origV[i]; }

    const auto &newV() const { return m_newV; }
    const Vec &newV(int i) const { return m_newV[i]; }

    std::tuple<coord_type, std::array<double, D>> grid_coord(const Vec &v) const;
    ColVecs compact_vertices() const;
    ColVecs V() const;// the non-grid velocities (from crossings)
    ColVecs all_V() const;// all of the vertices grid
    ColVecs all_GV() const;// all vertices in index/grid space
    VecVector stl_V() const;

    bool valid() const;
    bool valid_grid() const;

    std::vector<Vertex<D>> all_vertices() const;
    size_t Vsize() const { return newV().size(); }
    Vec V(int i) const;
    VType GV(int i) const;//vertices with grid masks
    mtao::map<int, int> vertex_reindexer() const;
    Vec vertex(int i) const;
    size_t num_vertices() const;

    Vertex<D> grid_vertex(int idx) const {//vertices with grid mask
        if (idx < grid_vertex_size()) {
            return vertex_unindex(idx);
        } else {
            return crossing(idx).vertex();
        }
    }
    std::string grid_info(int idx) const {
        if (idx < grid_vertex_size()) {
            std::stringstream ss;
            ss << "Grid vertex(";
            auto v = vertex_unindex(idx);
            std::copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, ","));
            ss << ")";
            return ss.str();
        } else {
            return std::string(crossing(idx));
        }
    }

    bool is_valid_cell_coord(const coord_type &c) const {
        return StaggeredGrid::cell_grid().valid_index(c);
    }
    bool is_valid_grid_index(int idx) const {
        return idx >= 0 && idx < StaggeredGrid::cell_size();
    }
    template<typename Derived>
    coord_type get_grid_cell(const Eigen::MatrixBase<Derived> &p) const {
        return std::get<0>(StaggeredGrid::coord(p));
    }
    void update_vertices_from_intersections();
    const mtao::map<Edge, int> &origEMap() const { return m_origEMap; }
    auto &&grid_vertices() const { return data().V(); }
    const CutData<D> &data() const { return m_data; }
    void update_vertices(const ColVecs &V, const std::optional<double> &threshold = -1);
    void update_vertices(const VecVector &V, const std::optional<double> &threshold = -1);
    void update_grid(const StaggeredGrid &indexer);
    std::set<Edge> edges() const;
    std::set<Edge> edge_slice(int dim, int slice) const;
    std::array<mtao::map<int, std::set<Edge>>, D> axial_edges() const;

    std::set<coord_type> possible_cells(const std::vector<int> &face) const;
    std::set<coord_type> possible_cells(const std::set<std::vector<int>> &face) const;
    std::set<coord_type> possible_cells_cell(const std::set<int> &faces, const std::vector<CutFace<D>> &) const;
    coord_mask<D> face_mask(const std::vector<int> &face) const;
    coord_mask<D> face_mask(const std::set<std::vector<int>> &face) const;
    bool is_in_cell(const std::vector<int> &face) const;
    bool is_in_cell(const std::set<std::vector<int>> &face) const;

    static VecVector colvecs_to_vecvector(const ColVecs &V);

  protected:
    // helper for assigning boundary information to

    using crossing_store_type = mtao::map<coord_type, std::set<EdgeCrossing<D>>>;
    std::array<crossing_store_type, D> get_per_axis_crossing_indices() const;
    std::array<crossing_store_type, D> get_per_axis_crossing_indices(const mtao::vector<Crossing<D>> &isects) const;

    //CutVertexMetas cut_vertex_metas() const;

    void extra_metadata(CutCellMesh<D> &mesh) const;


    //Note that grid-edge and grid-vertex values belong to multiple cells
    mtao::map<coord_type, std::set<int>> vertex_ownership() const;

    const GridDatab &active_grid_cell_mask() const {
        return m_active_grid_cell_mask;
    }
    bool is_active_cell(const coord_type &c) const;

  protected:

    CutData<D> m_data;
    mtao::map<Edge, int> m_origEMap;
    VecVector m_origV;

    bool m_filter_external_cells = false;


    //baked by bake_vertices
    VecVector m_newV;
    std::array<crossing_store_type, D> m_per_axis_crossings;
    std::vector<CrossingType> m_crossings;
    //baked by bake_edges
    Edges cut_edges;
    std::vector<CoordMaskedEdge<D>> m_grid_edges;
    std::vector<CutMeshEdge<D>> m_cut_edges;
    std::vector<CutMeshFace<D>> m_cut_faces;
    std::conditional_t<D == 2, mtao::geometry::mesh::HalfEdgeMesh, mtao::types::empty_t> m_hem;

    GridDatab m_active_grid_cell_mask;
};

template<int D>
class CutCellGenerator {};

}// namespace mandoline::construction

#include "mandoline/construction/generator_impl.hpp"
