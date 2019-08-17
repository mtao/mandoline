#pragma once
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/geometry/grid/grid_data.hpp>
#include <set>
#include "cutmesh.pb.h"
#include <mtao/eigen/stl2eigen.hpp>
namespace mandoline {
    template <int D>
    class CutCellMesh;
    namespace construction {
    template <int D>
    class CutCellGenerator;
    }

    class AdaptiveGrid: public mtao::geometry::grid::StaggeredGrid<double,3> {
        public:
            using Base = mtao::geometry::grid::StaggeredGrid<double,3>;
            using GridData3i = mtao::geometry::grid::GridDataD<int,3>;
            using coord_type = std::array<int,3>;
            using Edge = std::array<int,2>;
            struct Cell: public std::tuple<coord_type,int> {
                using Parent = std::tuple<coord_type,int>;
                using Parent::Parent;
                using Parent::operator=;
                const coord_type& corner() const { return std::get<0>(*this); }
                coord_type vertex(int a, int b, int c) const;
                int width() const { return std::get<1>(*this); }


                mtao::Vec3d center() const { return mtao::eigen::stl2eigen(corner()).array().cast<double>() + width()/2.0; }

                void  serialize(CutMeshProto::Cube&) const;
                static Cell from_proto(const CutMeshProto::Cube&);
            };
            struct Square: public std::tuple<coord_type,int,int> {
                using Parent = std::tuple<coord_type,int,int>;
                using Parent::Parent;
                Square(Square&&) = default;
                Square& operator=(Square&&) = default;
                Square(const Square&) = default;
                Square& operator=(const Square&) = default;
                using Parent::operator=;
                void  serialize(CutMeshProto::Square&) const;
                static Square from_proto(const CutMeshProto::Square&);
                const coord_type& corner() const { return std::get<0>(*this); }
                coord_type vertex(int a, int b) const;
                int width() const { return std::get<2>(*this); }
                int dimension() const { return std::get<1>(*this); }
                int axis() const { return std::get<1>(*this); }
            };
            struct Face: public Square {
                Face(Square&& s, const Edge& e): Square(std::move(s)), dual_edge(e) {}
                Face() = default;
                Face(Face&&) = default;
                Face& operator=(Face&&) = default;
                Face(const Face&) = default;
                Face& operator=(const Face&) = default;
                const Square& geometry() const { return *this; }
                void  serialize(CutMeshProto::Face&) const;
                static Face from_proto(const CutMeshProto::Face&);
                Edge dual_edge = {{-1,-1}};
                bool has_negative() const { return dual_edge[0] != -1; }
                bool has_positive() const { return dual_edge[1] != -1; }
                bool full_edge() const { return has_positive() && has_negative(); }
            };

            friend class CutCellMesh<3>;
            friend class construction::CutCellGenerator<3>;



            AdaptiveGrid() = default;
            AdaptiveGrid(const AdaptiveGrid&) = default;
            AdaptiveGrid(AdaptiveGrid&&) = default;
            AdaptiveGrid& operator=(const AdaptiveGrid&) = default;
            AdaptiveGrid& operator=(AdaptiveGrid&&) = default;
            AdaptiveGrid(const Base& b, const std::map<int,Cell>& cells = {}): Base(b), m_cells(cells) {make_faces();}
            std::array<int,4> face(const Cell& c, int axis, bool sign) const;
            std::array<int,4> face(int idx, int axis, bool sign) const;
            mtao::ColVecs3i triangulated(int idx) const;
            mtao::ColVecs3i triangulated(const Cell& c) const;
            mtao::ColVecs3d boundary_centroids() const;
            void cell_centroids(mtao::ColVecs3d&) const;
            GridData3i grid() const;
            //const std::vector<Edge>& boundary() const { return m_boundary; }
            std::vector<Edge> boundary_pairs() const;
            const std::map<int,Cell>& cells() const { return m_cells; }
            const Cell& cell(int idx) const { return m_cells.at(idx); }
            mtao::VecXd dual_edge_lengths() const;
            mtao::VecXd face_volumes() const;
            mtao::VecXd cell_volumes() const;
            std::vector<Eigen::Triplet<double>> boundary_triplets(int min_edge_index) const;
            int num_faces() const;
            int num_cells() const;
            std::vector<Face> faces(const GridData3i& grid) const;
            std::vector<Face> faces() const { return m_faces; }
            void make_faces();
            static GridData3i grid_from_cells(const coord_type& shape, const std::map<int,Cell>& cells);
            std::vector<Eigen::Triplet<double>> grid_face_projection(int min_edge_index) const;
            std::vector<Eigen::Triplet<double>> grid_cell_projection() const;


            inline bool is_valid_edge(const Edge& e) const {
                auto [a,b] = e;
                return a != b && m_cells.find(a) != m_cells.end() && m_cells.find(b) != m_cells.end();
            }

        private:

            //std::vector<Edge> m_boundary;//Beware of -1!
            std::vector<Face> m_faces;
            std::vector<Edge> m_edges;
            std::map<int,Cell> m_cells;
    };

    class AdaptiveGridFactory {
        public:
            constexpr static int logwidth = 1;
            constexpr static int width = 1 << logwidth;
            using Edge = std::array<int,2>;
            using GridData3i = mtao::geometry::grid::GridDataD<int,3>;
            using ActiveMask = std::bitset<1 << logwidth * 3>;
            using GridData3 = mtao::geometry::grid::GridDataD<bool,3>;
            using GridData3b = mtao::geometry::grid::GridDataD<ActiveMask,3>;
            using coord_type = std::array<int,3>;
            using Cell = AdaptiveGrid::Cell;
            using AxialBEdgeMap = std::array<std::map<int,std::set<Edge>>,3> ;

            AdaptiveGridFactory(const GridData3& mask);



            std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> compute_axial_edges(const std::optional<int>& max_level = {}) const ;
            std::tuple<std::set<Edge>,AxialBEdgeMap> compute_edges(const std::optional<int>& max_level = {}) const ;

            void make_cells(const std::optional<int>& max_level = {});

            using Indexer = mtao::geometry::grid::indexing::OrderedIndexer<3>;

            const static Indexer cmask_indexer;
            const static Indexer vmask_indexer;
            const static std::array<Indexer,3> mask_edge_indexers;
            const static std::array<coord_type,3> mask_edge_shapes;

            Indexer vertex_indexer;

            static int cmask_index(const coord_type& c) { return cmask_indexer.index(c); }
            static int cmask_index(int a, int b, int c) { return cmask_indexer.index(a,b,c); }
            static int vmask_index(const coord_type& c) { return vmask_indexer.index(c); }
            static int vmask_index(int a, int b, int c) { return vmask_indexer.index(a,b,c); }
            static int emask_index(int dim, const coord_type& c) { return mask_edge_indexers[dim].index(c); }
            static int emask_index(int dim, int a, int b, int c) { return mask_edge_indexers[dim].index(a,b,c); }
            int vertex_index(const coord_type& c) const { return vertex_indexer.index(c); }
            int vertex_index(int a, int b, int c) const { return vertex_indexer.index(a,b,c); }


            void make_edges(const std::optional<int>& max_level = {}) ;
            void make_cells(const ActiveMask& mask, int level, const coord_type& coord) ;
            void make_cells(const GridData3& mask, int level);

            void add_cell(const coord_type& corner, int jump);

            GridData3i grid_from_cells(const std::map<int,Cell>& cells) const;

            std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> get_edges(const std::array<GridData3,3>& edge_masks, int level, const coord_type& offset = {}) const;
            Edge get_edge(const coord_type& start, int jump, int dim) const;

            int get_jump(int level) const { return  1 << (level*logwidth); }
            int get_parent_jump(int level) const { return  1 << ((level+1)*logwidth); }

            std::array<coord_type,3> make_edge_shapes(const coord_type& coord) const;

            std::array<GridData3,3> empty_edge_masks(const coord_type& shape,bool value=false) const;
            std::array<GridData3,3> empty_edge_masks(int level, bool value=false) const;

            std::array<GridData3,3> make_edge_mask(const GridData3& mask, int level) const;
            std::array<GridData3,3> make_edge_mask(const GridData3b& mask, int level) const;
            std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> make_edges(const GridData3& mask, int level) const;
            std::tuple<std::array<std::set<Edge>,3>,AxialBEdgeMap> make_edges(const GridData3b& mask, int level) const;


            AdaptiveGrid create() const;


        public:
            mtao::ColVecs2i edges ;
            mtao::ColVecs2i boundary_edges ;
            std::vector<GridData3b> levels;
            std::vector<GridData3> levels_mask;
            GridData3 original;
            std::map<int,Cell> cells;

        private:
            template <typename GridAccessor>
                void set_edge_masks(const coord_type& shape, const coord_type& corner, std::array<GridData3,3>& edge_mask , const GridAccessor& accessor) const;
    };

}
