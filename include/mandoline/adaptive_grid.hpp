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
            using Vec = typename Base::Vec;
            struct Cell: public std::tuple<coord_type,int> {
                using Parent = std::tuple<coord_type,int>;
                using Parent::Parent;
                using Parent::operator=;
                const coord_type& corner() const { return std::get<0>(*this); }
                coord_type vertex(int a, int b, int c) const;
                int width() const { return std::get<1>(*this); }


                mtao::Vec3d center() const { return mtao::eigen::stl2eigen(corner()).array().cast<double>() + width()/2.0; }

                void  serialize(protobuf::Cube&) const;
                static Cell from_proto(const protobuf::Cube&);
                //local coordinates
                bool is_inside(const Vec& p) const;
            };
            struct Square: public std::tuple<coord_type,int,int> {
                using Parent = std::tuple<coord_type,int,int>;
                using Parent::Parent;
                Square(Square&&) = default;
                Square& operator=(Square&&) = default;
                Square(const Square&) = default;
                Square& operator=(const Square&) = default;
                using Parent::operator=;
                void  serialize(protobuf::Square&) const;
                static Square from_proto(const protobuf::Square&);
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
                void  serialize(protobuf::Face&) const;
                static Face from_proto(const protobuf::Face&);
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
            GridData3i cell_ownership_grid() const;
            //const std::vector<Edge>& boundary() const { return m_boundary; }
            std::vector<Edge> boundary_pairs() const;
            const std::map<int,Cell>& cells() const { return m_cells; }
            const Cell& cell(int idx) const { return m_cells.at(idx); }
            mtao::VecXd dual_edge_lengths() const;
            mtao::VecXd face_volumes() const;
            mtao::VecXd cell_volumes() const;
            std::vector<Eigen::Triplet<double>> boundary_triplets(int min_edge_index) const;
            int num_edges() const;
            int num_faces() const;
            int num_cells() const;
            mtao::ColVecs2i edges() const;
            std::vector<Face> faces(const GridData3i& grid) const;
            const std::vector<Face>& faces() const { return m_faces; }
            void make_faces();
            static GridData3i grid_from_cells(const coord_type& shape, const std::map<int,Cell>& cells);
            std::vector<Eigen::Triplet<double>> grid_face_projection(int min_edge_index) const;
            std::vector<Eigen::Triplet<double>> grid_cell_projection() const;

            int get_cell_index(const Vec& p) const;

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


}
