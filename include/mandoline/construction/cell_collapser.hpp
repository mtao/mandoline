#pragma once

#include <array>
#include <tuple>
#include <map>
#include <set>
#include <vector>
#include <mtao/data_structures/disjoint_set.hpp>
#include <mtao/geometry/mesh/triangle_fan.hpp>
#include <mtao/geometry/trigonometry.hpp>
#include "mandoline/cutface.hpp"
#include <mtao/logging/logger.hpp>
#include <mtao/logging/profiler.hpp>


namespace mandoline::construction {
    struct CellCollapser {
        public:
            using Edge = std::array<int,2>;
            using HalfFace = std::tuple<int,Edge>;
            using CoordType = std::array<int,3>;
            //CellCollapser(const mtao::map<int,std::set<std::vector<int>>>& faces);
            CellCollapser(const mtao::map<int,CutFace<3>>& faces);
            std::tuple<std::set<std::vector<int>>,std::set<Edge>> clean_unmanifold(const std::vector<int>&);
            mtao::map<Edge,std::set<const HalfFace*>> collect_halffaces() const;

            int dual_cell(const HalfFace& hf) const;
            int cell(const HalfFace& hf) const;

            template <typename Derived>
                void merge(const Eigen::MatrixBase<Derived>& V);
            template <typename Derived>
                void bake(const Eigen::MatrixBase<Derived>& V);

            std::vector<std::set<int>> cell_faces() const;
            const auto& faces() const { return m_faces; }
            const auto& face(int idx) const { return m_faces.at(idx); }

            const auto& cell_boundaries() const { return m_cell_boundaries; }
            void remove_boundary_cells_from_vertices(const std::set<int>& boundary_vertices);
            void remove_boundary_cells_from_faces(const std::set<int>& boundary_faces);
            std::set<int> folded_faces() const; //faces that have their duals as well
            //private:
            mtao::map<HalfFace,std::tuple<int,bool>> m_halfface_to_cell;
            mtao::map<int,CutFace<3>> m_faces;
            std::map<int,std::set<Edge>> flap_edges;
            mtao::map<int,std::set<CoordType>> face_cell_possibilities;
            std::vector<mtao::map<int,bool>> m_cell_boundaries;
            std::set<Edge> invalid_edges;
            mtao::data_structures::DisjointSet<int> cell_ds;

    };

    template <typename Derived>
        void CellCollapser::bake(const Eigen::MatrixBase<Derived>& V) {
            merge(V);

        }
    template <typename Derived>
        void CellCollapser::merge(const Eigen::MatrixBase<Derived>& V) {
            auto t2 = mtao::logging::profiler("cell collapser merge",false,"profiler");
            for(auto [e,hfs]: collect_halffaces()) {
                //auto t3 = mtao::logging::profiler("cell collapser edge fan processessing",false,"profiler");
                if(hfs.empty()) {
                    continue;
                }
                auto [a,b] = e;
                if(a > b) {
                    continue;
                }
                auto va = V.col(a);
                auto vb = V.col(b);
                mtao::Vec3d t = (va - vb).normalized();
                mtao::Matrix<double,3,2> uv;
                auto u = uv.col(0);
                auto v = uv.col(1);
                for(int i = 0; i < 3; ++i) {
                    if(t((i+1)%3) != 0) {
                        u = t.cross(mtao::Vec3d::Unit(i));
                        break;
                    }
                }
                v = u.cross(t);
                mtao::map<double,HalfFace> index_map;
                for(auto&& hfp: hfs) {
                    auto [fidx,e] = *hfp;

                    auto&& F = m_faces.at(fidx);
                    auto& N = F.N;
                    bool sign = std::get<1>(m_halfface_to_cell.at(*hfp));

                    mtao::Vec2d p = uv.transpose() * (sign?1:-1) * N;
                    double ang = mtao::geometry::trigonometry::angle(p)(0);
                    index_map[ang] = *hfp;
                }
                auto it = index_map.begin();
                auto it1 = it;
                it1++;
                for(; it != index_map.end(); ++it, ++it1) {
                    if(it1 == index_map.end()) {
                        it1 = index_map.begin();
                    }
                    cell_ds.join(cell(it->second),dual_cell(it1->second));
                }

            }
            cell_ds.reduce_all();
            mtao::map<int,int> reindexer;
            for(auto&& i: cell_ds.root_indices()) {
                reindexer[cell_ds.node(i).data] = reindexer.size();
            }
            m_cell_boundaries.resize(reindexer.size());
            for(auto&& [hf,cb]: m_halfface_to_cell) {
                auto& [c,b] = cb;
                int root = cell_ds.get_root(c).data;
                int cell = reindexer[root];
                c = cell;
                m_cell_boundaries[cell][std::get<0>(hf)] = b;
            }

            if(!face_cell_possibilities.empty()) { 
                /*
                   mtao::map<CoordType,std::tuple<int,bool>> map;
                   for(auto&& [c,fs]: cell_boundaries) {
                   for(auto&& [f,b]: fs) {

                   }
                   }
                   */
            }



            m_cell_boundaries.erase(std::remove_if(m_cell_boundaries.begin(),m_cell_boundaries.end(),[](auto&& m) {
                        return m.size() < 4;
                        }), m_cell_boundaries.end());
        }
}
