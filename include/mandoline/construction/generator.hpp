#pragma once
#include <vector>
#include <tuple>
#include <mtao/types.h>
#include <bitset>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/functional.hpp>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/logging/timer.hpp>
#include <map>
#include <set>
#include "mandoline/mesh2.hpp"
#include "mandoline/mesh3.hpp"
#include "mandoline/construction/cutdata.hpp"
#include "mandoline/cutface.hpp"
#include <iterator>
#include "mandoline/adaptive_grid.hpp"

#define USE_INTERIOR_MASK

namespace mandoline::construction {

    //ASSUMES SIMPLICIAL INPUTS
    template <int D> 
        class CutCellEdgeGenerator: public mtao::geometry::grid::StaggeredGrid<double,D> {
            public:
                constexpr static double threshold_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
                using Vec = mtao::Vector<double,D>;
                using Veci = mtao::Vector<int,D>;
                using ColVecs = mtao::ColVectors<double,D>;
                using Edge = std::array<int,2>;
                using Edges = mtao::ColVectors<int,2>;
                using BoundaryElements = mtao::ColVectors<int,D>;
                using VecVector = mtao::vector<Vec>;
                using StaggeredGrid = mtao::geometry::grid::StaggeredGrid<double,D>;
                using GridType = mtao::geometry::grid::GridD<double,D>;
                using StaggeredGrid::grid;
                using StaggeredGrid::vertex_grid;
                using StaggeredGrid::vertex_shape;
                using StaggeredGrid::vertex_unindex;
                using StaggeredGrid::cell_grid;
                using StaggeredGrid::origin;
                using StaggeredGrid::dx;
                using StaggeredGrid::shape;
                using StaggeredGrid::vertex_index;
                using StaggeredGrid::cell_index;
                using StaggeredGrid::cell_size;
                using StaggeredGrid::vertex;
                using CoordType = std::array<int,D>;
                //vertex that lives in the grid as a grid cell + local offset
                using VType = Vertex<D>;
                using EdgeIntersectionType = EdgeIntersection<D>;
                using CrossingType = Crossing<D>;
                using CoordMaskedEdgeType = CoordMaskedEdge<D>;
                using GridDatab = mtao::geometry::grid::GridDataD<bool,D>;

                struct GridVertex: public CoordType {
                    GridVertex(const StaggeredGrid& sg, const CoordType& c): CoordType(c), index(sg.vertex_index(c)) {}
                    GridVertex(const StaggeredGrid& sg, int index): CoordType(sg.template form_unindex<0>(index)), index(index) {}
                    int index = -1;
                };

                bool performance = false;


                //returns a new mesh and a set of boundary vertices
                std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>> compute_planar_hem(const ColVecs& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const;
                //returns a new mesh and a set of boundary vertices
                std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>> compute_planar_hem(const std::vector<VType>& GV, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const;
                std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>> compute_planar_hem(const std::vector<VType>& GV, const ColVecs& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const;


                virtual void bake();
                virtual void bake_vertices();
                void bake_active_grid_cell_mask();
                virtual void bake_edges();
                virtual void bake_faces() {}
                virtual void bake_cells() {}
                
                virtual void clear();

                CutCellMesh<D> generate_edges() const;
                CutCellMesh<D> generate_faces() const;
                CutCellMesh<D> generate() const;

                //This threshold is to fuse vertices near grid vertices. for the results of float to double conversions 1e-6 seemed reasonable
                CutCellEdgeGenerator(const VecVector& V, const StaggeredGrid& grid, std::optional<double> threshold  = 1e-6);
                CutCellEdgeGenerator(const ColVecs& V, const StaggeredGrid& grid, std::optional<double> threshold  = 1e-6);
                CutCellEdgeGenerator(const StaggeredGrid& grid);
                CutCellEdgeGenerator() = default;
                CutCellEdgeGenerator(CutCellEdgeGenerator&&) = default;
                CutCellEdgeGenerator& operator=(CutCellEdgeGenerator&&) = default;
                //CutCellGenerator(const VecVector& V, const Vec& dx = Vec::Ones());




                template <typename Derived>
                    VType get_vertex(const Eigen::MatrixBase<Derived>& p) const {
                        auto&& g = vertex_grid();
                        Vec local = p.cwiseQuotient(g.dx()) - g.origin().cwiseQuotient(g.dx());
                        return VType::from_vertex(local);
                    }
                Vec get_world_vertex(const VType& p) const {
                    auto&& g = vertex_grid();
                    return g.origin() + (p.p().cwiseProduct(g.dx()));

                }
                Vec get_world_vertex(const GridVertex& p) const {
                    auto&& g = vertex_grid();
                    return g.origin() + (mtao::eigen::stl2eigen(p).cwiseProduct(g.dx()));
                }

                size_t grid_vertex_size() const { return StaggeredGrid::vertex_size(); }

                bool is_grid_vertex(int idx) const { return idx < grid_vertex_size(); }
                bool is_orig_vertex(int idx) const { return !is_grid_vertex(idx) && idx < grid_vertex_size() + origV().size(); }
                bool is_new_vertex(int idx) const { return !is_orig_vertex(idx) && idx < num_vertices(); }




                void reset(const mtao::vector<VType>& grid_vertices);

                void add_boundary_elements(const BoundaryElements& E);
                void add_edges(const Edges& E);
                static std::set<EdgeIntersectionType> cell_edge_intersections(const std::set<CoordType>& cells) ;
                template <typename Func>
                    static void per_cell_vertex_looper(Func&& f, const CoordType& c) {
                        for(int i = 0; i < (2 << D); ++i) {
                            CoordType cc = c;
                            std::bitset<D> bs(i);
                            for(int j = 0; j < D; ++j) {
                                cc[j] += bs[j]?1:0;
                            }
                            f(cc,bs);

                        }
                    }
                template <typename Func>
                    static void per_dual_cell_vertex_looper(Func&& f, const CoordType& c) {
                        for(int i = 0; i < (2 << D); ++i) {
                            CoordType cc = c;
                            std::bitset<D> bs(i);
                            for(int j = 0; j < D; ++j) {
                                cc[j] -= bs[j]?1:0;
                            }
                            f(cc,bs);

                        }
                    }
                template <typename Func>
                    static void per_boundary_cell_vertex_looper(int N, Func&& f, const CoordType& c) {
                        for(int i = 0; i < (2 << (D-1)); ++i) {
                            std::bitset<D> bs(i);
                            if(bs[N]) {//if bit is 1 we move on
                                continue;
                            }
                            CoordType cc = c;
                            Vec v = Vec::Zero();
                            for(int j = 0; j < D; ++j) {
                                cc[j] += bs[j];
                            }
                            f(cc,bs);

                        }
                    }


                coord_mask<D> get_mask(int idx) const {
                    int gvs = grid_vertex_size();
                    using CM = coord_mask<D>;
                    if(idx < gvs) {
                        return CM{StaggeredGrid::template staggered_unindex<0,0>(idx)};
                    } else {
                        return crossing(idx).mask();
                    }
                }

                auto&& crossing(int idx) const { return m_crossings[idx - grid_vertex_size()]; }

                size_t total_vertex_size() const { return m_crossings.size() + grid_vertex_size(); }
                auto&& crossings() const { return m_crossings; }

                size_t new_vertex_offset() const { return grid_vertex_size() + origV().size();}

                std::set<CoordType> active_cells() const;
                std::array<mtao::map<CoordType,std::set<int>>,D> crossing_indices() const;

                auto&& origE() const { return data().E(); }
                auto origE(int idx) const { return data().E(idx); }




                const auto& origV() const { return m_origV; }
                const Vec& origV(int i) const { return m_origV[i]; }

                const auto& newV() const { return m_newV; }
                const Vec& newV(int i) const { return m_newV[i]; }

                std::tuple<CoordType,std::array<double,D>> grid_coord(const Vec& v) const;
                ColVecs compact_vertices() const;
                ColVecs V() const;
                ColVecs all_V() const;
                ColVecs all_GV() const;
                VecVector stl_V() const;

                size_t Vsize() const { return newV().size(); }
                Vec V(int i) const;
                VType GV(int i) const;
                mtao::map<int,int> vertex_reindexer() const;
                Vec vertex(int i) const;
                size_t num_vertices() const;

                Vertex<D> grid_vertex(int idx) const {
                    if(idx < grid_vertex_size()) {
                        return vertex_unindex(idx);
                    } else {
                        return crossing(idx).vertex();
                    }
                }
                std::string grid_info(int idx) const {
                    if(idx < grid_vertex_size()) {
                        std::stringstream ss;
                        ss << "Grid vertex(";
                        auto v = vertex_unindex(idx);
                        std::copy(v.begin(),v.end(),std::ostream_iterator<int>(ss,","));
                        ss << ")";
                        return ss.str();
                    } else {
                        return std::string(crossing(idx));
                    }
                }

                bool is_valid_grid_index(int idx) const {
                    return idx >= 0 && idx < StaggeredGrid::cell_size();
                }
                template <typename Derived>
                    CoordType get_grid_cell(const Eigen::MatrixBase<Derived>& p) const {
                        return std::get<0>(StaggeredGrid::coord(p));
                    }
                void update_vertices_from_intersections();
                const mtao::map<Edge, int>& origEMap() const { return m_origEMap; }
                auto&& grid_vertices() const { return data().V(); }
                const CutData<D>& data() const { return m_data; }
                void update_vertices(const ColVecs& V, const std::optional<double>& threshold = -1);
                void update_vertices(const VecVector& V, const std::optional<double>& threshold = -1);
                void update_grid(const StaggeredGrid& indexer);
                std::set<Edge> edges() const;
                std::set<Edge> edge_slice(int dim, int slice) const;
                std::array<mtao::map<int,std::set<Edge>>,D> axial_edges() const;

                std::set<CoordType> possible_cells(const std::vector<int>& face)const ;
                std::set<CoordType> possible_cells(const std::set<std::vector<int>>& face)const ;
                std::set<CoordType> possible_cells_cell(const std::set<int>& faces, const std::vector<CutFace<D>>&)const ;
                coord_mask<D> face_mask(const std::vector<int>& face)const ;
                coord_mask<D> face_mask(const std::set<std::vector<int>>& face)const ;
                bool is_in_cell(const std::vector<int>& face)const ;
                bool is_in_cell(const std::set<std::vector<int>>& face)const ;

                static VecVector colvecs_to_vecvector(const ColVecs& V);
            private:
                using crossing_store_type = mtao::map<CoordType,std::set<EdgeCrossing<D>>>;
                std::array<crossing_store_type,D> get_per_axis_crossing_indices() const;
                std::array<crossing_store_type,D> get_per_axis_crossing_indices(const mtao::vector<Crossing<D>>& isects) const;

                //CutVertexMetas cut_vertex_metas() const;

                void extra_metadata(CutCellMesh<D>& mesh) const;


                //Note that grid-edge and grid-vertex values belong to multiple cells
                mtao::map<CoordType,std::set<int>> vertex_ownership()const ;

                const GridDatab& active_grid_cell_mask() const {
                    return m_active_grid_cell_mask;
                }

            protected:
                CutData<D> m_data;
                mtao::map<Edge, int> m_origEMap;
                VecVector m_origV;
                //baked by bake_vertices
                VecVector m_newV;
                std::array<crossing_store_type,D> m_per_axis_crossings;
                std::vector<CrossingType> m_crossings;
                //baked by bake_edges
                std::set<CoordMaskedEdgeType> m_edges;
                Edges cut_edges;
                std::vector<CutMeshFace<D>> m_cut_faces;

                GridDatab m_active_grid_cell_mask;

        };

    template <int D>
        class CutCellGenerator;
    template <>
        class CutCellGenerator<2>: public CutCellEdgeGenerator<2> {
            public:
                using CCEG = CutCellEdgeGenerator<2>;
                using BoundaryElements = typename CCEG::BoundaryElements;
                using CCEG::CCEG;
                using CoordType = typename CCEG::CoordType;
                CutCellGenerator() = default;
                CutCellGenerator(CutCellGenerator&&) = default;
                CutCellGenerator& operator=(CutCellGenerator&&) = default;
                CutCellMesh<2> generate() const ;
                void add_boundary_elements(const BoundaryElements& E);
                void bake_faces() override;
                void extra_metadata(CutCellMesh<2>& mesh) const;
                mtao::geometry::mesh::HalfEdgeMesh hem;
        };
    template <>
        class CutCellGenerator<3>: public CutCellEdgeGenerator<3> {
            public:
                using CCEG = CutCellEdgeGenerator<3>;
                using BoundaryElements = typename CCEG::BoundaryElements;
                using CoordType = typename CCEG::CoordType;
                using CCEG::add_boundary_elements;
                using CCEG::CCEG;
                CutCellGenerator() = default;
                CutCellGenerator(CutCellGenerator&&) = default;
                CutCellGenerator& operator=(CutCellGenerator&&) = default;
                using CCEG::shape;

                size_t Vsize() const {
                    return CCEG::Vsize();
                }
                ~CutCellGenerator();

                const mtao::map<int,CutFace<D>>& faces()const { return m_faces; }

                void compute_faces();
                void compute_faces_vertex();
                void compute_faces_axis(int idx);
                mtao::map<int,CutFace<D>> compute_faces_axis(int idx, int cidx)const;

                void bake_faces() override;
                void bake_cells() override;
                void update_active_grid_cell_mask();
                bool adaptive = true;
                std::optional<int> adaptive_level = 0;
                std::optional<AdaptiveGrid> adaptive_grid;
                std::optional<std::map<int,int>> adaptive_grid_regions;
                void set_region_map(const std::map<int,std::set<int>>& region_vertices);
                void bake() override;
                void clear() override;

                Edge smallest_ordered_edge(const std::vector<int>& v) const;



                static mtao::ColVecs3i faceMap_to_faces(mtao::map<int,mtao::ColVecs3i>& fm);

                void update_vertices_from_intersections();
                CutCellMesh<3> generate() const ;
                //void make_faces(const CutCellMesh<3>& ccm ) const;

                void extra_metadata(CutCellMesh<3>& mesh) const;
                /*
                   const auto& newV() const { return m_newV; }
                   const Vec& newV(int i) const { return m_newV[i]; }
                   auto&& stl_face_cutE() const { return m_face_cutedges; }
                   const BoundaryElements& origF() const { return m_origF;}
                   auto origF(int i) const { return m_origF.col(i);}
                   const BoundaryElements& origFE() const { return m_cell_edge_map;}
                   auto origFE(int i) const { return m_cell_edge_map.col(i);}
                   */

                //baked by bake_faces
                struct AxisHEMData {
                    using GridDatab = mtao::geometry::grid::GridDataD<bool,2>;
                    GridDatab active_grid_cell_mask;
                    std::set<Edge> edges;
                    std::set<Edge> boundary_edges;
                    mtao::geometry::mesh::HalfEdgeMesh hem;
                    bool is_boundary_cell(const std::vector<int>& verts) const;

                };
                bool check_cell_containment() const;
                bool check_face_utilization() const;

                mtao::Vec3d area_normal(const std::vector<int>& F) const;
                mtao::Vec3d area_normal(const std::set<std::vector<int>>& F) const;
                std::set<Edge> edge_slice(int dim, int slice) const;
                mtao::map<int,int> cut_cell_to_primal_map;
                mtao::ColVecs3d origN;
                std::array<mtao::map<int,AxisHEMData>,3> axis_hem_data;
                std::array<std::set<Edge>,3> axial_primal_faces;
                BoundaryElements m_newF;

                mtao::map<int,CutFace<D>> m_faces;
                std::set<int> mesh_face_indices;
                std::array<std::set<int>,3> axis_face_indices;
                std::set<int> folded_faces;
                std::set<Edge> adaptive_edges;
                std::set<Edge> adaptive_bedges;


                std::vector<CutCell> cell_boundaries;
                std::set<int> boundary_vertices;
                std::set<int> boundary_faces;
        };

    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate_faces() const;
    template <>
        CutCellMesh<2> CutCellEdgeGenerator<2>::generate() const;
    template <>
        CutCellMesh<3> CutCellEdgeGenerator<3>::generate() const;
    template <>
        auto CutCellEdgeGenerator<2>::compute_planar_hem(const ColVecs& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const-> std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>>;
    template <>
        auto CutCellEdgeGenerator<2>::compute_planar_hem(const std::vector<VType>& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const-> std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>>;
    template <>
        auto CutCellEdgeGenerator<2>::compute_planar_hem(const std::vector<VType>& GV, const ColVecs& V, const Edges& E, const GridDatab& interior_cell_mask, int cell_size) const-> std::tuple<mtao::geometry::mesh::HalfEdgeMesh,std::set<Edge>>;
}




#include "mandoline/construction/generator_impl.hpp"

