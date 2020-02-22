#pragma once

#include <array>
#include <tuple>
#include <numeric>
#include <map>
#include <set>
#include <vector>
#include <mtao/data_structures/disjoint_set.hpp>
#include <mtao/geometry/mesh/triangle_fan.hpp>
#include <mtao/geometry/trigonometry.hpp>


namespace mandoline::construction {
    // collapsers a collection of edges into faces
    struct FaceCollapser {
        public:
            using Edge = std::array<int,2>;
            using CoordType = std::array<int,3>;
            //FaceCollapser(const std::map<int,std::set<std::vector<int>>>& faces);
            FaceCollapser(const std::set<Edge>& edges);
            std::map<int,std::set<int>> collect_edges() const;

            int dual_face(Edge e) const;
            int face(const Edge& e) const;


            // merge edges into faces
            template <typename Derived>
                void merge(const Eigen::MatrixBase<Derived>& V);
            // bake the face structure. for now it is only merge, but could be more?
            template <typename Derived>
                void bake(const Eigen::MatrixBase<Derived>& V);

            // cleans up the face identities to start from 0
            void finalize();

            std::map<int,std::vector<int>> faces() const;
            std::map<int,std::map<int,int>> face_adjacency_map() const;
            std::map<int,std::map<Edge,Edge>> face_dual_adjacency_map() const;

            // indicate that a loop of vertices, interpreted as directed edges [i,i+1] are on the outside
            void set_edges_for_removal(const std::vector<int>& boundary_loop);
            // indicate that a single directed edges in on the outside
            void set_edge_for_removal(const Edge& e);
            private:
            // a directed edge maps to a face identity and whether this is the same order as the input
            std::map<Edge,std::tuple<int,bool>> m_edge_to_face;
            // structure for combining face identities
            mtao::data_structures::DisjointSet<int> face_ds;
            // face to face map
            std::map<Edge,Edge> dual_edge_graph;

    };

    template <typename Derived>
        void FaceCollapser::bake(const Eigen::MatrixBase<Derived>& V) {
            merge(V);

        }
    template <typename Derived>
        void FaceCollapser::merge(const Eigen::MatrixBase<Derived>& V) {
            for(auto [a,bs]: collect_edges()) {

                std::vector<int> indices(bs.begin(),bs.end());
                auto va = V.col(a);
                mtao::ColVecs2d D(2,bs.size());
                for(auto [i,j]: mtao::iterator::enumerate(indices)) {
                    D.col(i) = V.col(j) - va;
                }
                std::vector<char> quadrants(bs.size());
                constexpr static std::array<int,4> __quadrants{{4,1,3,2}};
                for(int i = 0; i < bs.size(); ++i) {
                    auto b = D.col(i);
                    // ++ +- -+ -- => 1 4 2 3
                    quadrants[i] = __quadrants[2 * std::signbit(b.y()) + std::signbit(b.x())];
                }
                // sort by quadrant and then by cross product volume
                auto comp = [&](int ai, int bi) -> bool {
                    const char qa = quadrants[ai];
                    const char qb = quadrants[bi];
                    if(qa == qb) {
                        auto a = D.col(ai);
                        auto b = D.col(bi);
                        return b.x() * a.y() < a.x() * b.y();
                    } else {
                        return qa < qb;
                    }
                };
                std::vector<int> ordered_indices(bs.size());
                std::iota(ordered_indices.begin(),ordered_indices.end(),0);
                std::sort(ordered_indices.begin(),ordered_indices.end(),comp);
                std::transform(ordered_indices.begin(),ordered_indices.end(),ordered_indices.begin(),[&](int idx) -> int { return indices[idx]; });
                auto it = ordered_indices.begin();
                auto it1 = it;
                it1++;
                for(; it != ordered_indices.end(); ++it, ++it1) {
                    if(it1 == ordered_indices.end()) {
                        it1 = ordered_indices.begin();
                    }
                    Edge e{{*it1,a}};
                    Edge ne{{a,*it}};
                    int triangle_indexx = face(e);
                    int ntriangle_indexx = face(ne);
                    face_ds.join(triangle_indexx,ntriangle_indexx);
                    dual_edge_graph[e] = ne;
                }

            }
            finalize();
        }
}
