#pragma once

#include <array>
#include <tuple>
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
            // takes a pair [ quadrant, vertex index ]
            auto comp = [&](const std::array<int,2>& a, const std::array<int,2>& b) -> bool {
                if(a[0]==b[0]) {
                    auto va = V.col(a[1]);
                    auto vb = V.col(b[1]);
                    return 
                } else {
                    return a[0] < b[0];
                }
            };
            for(auto [a,bs]: collect_edges()) {

                std::map<double,int> index_map;
                auto va = V.col(a);
                for(auto&& b: bs) {
                    auto vb = V.col(b);

                    mtao::Vec2d p = vb - va;
                    double ang = mtao::geometry::trigonometry::angle(p)(0);
                    index_map[ang] = b;
                }
                auto it = index_map.begin();
                auto it1 = it;
                it1++;
                for(; it != index_map.end(); ++it, ++it1) {
                    if(it1 == index_map.end()) {
                        it1 = index_map.begin();
                    }
                    Edge e{{it1->second,a}};
                    Edge ne{{a,it->second}};
                    int triangle_indexx = face(e);
                    int ntriangle_indexx = face(ne);
                    face_ds.join(triangle_indexx,ntriangle_indexx);
                    dual_edge_graph[e] = ne;
                }

            }
            finalize();
        }
}
