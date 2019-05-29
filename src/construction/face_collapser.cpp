#include "mandoline/construction/face_collapser.hpp"
#include <algorithm>
#include <iterator>

#include <mtao/iterator/enumerate.hpp>

namespace mandoline::construction {
    FaceCollapser::FaceCollapser(const std::set<Edge>& edges) {
        for(auto&& [cid,e]: mtao::iterator::enumerate(edges)) {
            face_ds.add_node(2*cid);
            face_ds.add_node(2*cid+1);
            {
                m_edge_to_face[e] = std::make_tuple(2*cid,true);
            }
            Edge e2 = e;
            std::swap(e2[0],e2[1]);
            {
                m_edge_to_face[e2] = std::make_tuple(2*cid+1,false);
            }
        }
    }

    auto FaceCollapser::collect_edges() const -> std::map<int,std::set<int>> {
        std::map<int,std::set<int>> ret;
        for(auto&& [e,pr]: m_edge_to_face) {
            ret[e[0]].insert(e[1]);
        }
        return ret;
    }
    int FaceCollapser::dual_face(Edge e) const {
        std::swap(e[0],e[1]);

        return face(e);
    }
    int FaceCollapser::face(const Edge& e) const {
        return std::get<0>(m_edge_to_face.at(e));
    }


    void FaceCollapser::set_edges_for_removal(const std::vector<int>& boundary_loop) {
        face_ds.add_node(-1);
        for(int i = 0; i < boundary_loop.size(); ++i) {
            int j = (i+1)%boundary_loop.size();
            Edge e{{i,j}};
            if(m_edge_to_face.find(e) != m_edge_to_face.end()) {
                face_ds.join(std::get<0>(m_edge_to_face[e]),-1);
            }
        }
    }

    void FaceCollapser::set_edge_for_removal(const Edge& e) {
        face_ds.add_node(-1);

        if(m_edge_to_face.find(e) != m_edge_to_face.end()) {
            face_ds.join(std::get<0>(m_edge_to_face[e]),-1);
        }
    }

    std::map<int,std::map<int,int>> FaceCollapser::face_adjacency_map() const {
        std::map<int,std::map<int,int>> ret;

        for(auto&& [e,pr]: m_edge_to_face) {
            auto&& [a,b] = e;
            auto&& [c,s] = pr;
            ret[c][a] = b;
        }
        return ret;
    }
    std::map<int,std::vector<int>> FaceCollapser::faces() const {
        std::map<int,std::vector<int>> ret;

        auto deg = dual_edge_graph;
        auto inc = [&](auto&& it) {
            return deg.find(it->second);
        };
        std::set<Edge> seen_edges;
        for(auto&& [a,b]: dual_edge_graph) {
            seen_edges.insert(a);
        }

        for(auto it = deg.begin(); it != deg.end(); ++it) {
            if(seen_edges.find(it->first) == seen_edges.end()) {
                continue;
            }
            int myface = face(it->first);
            if(myface < 0) {
                seen_edges.erase(it->first);
                continue;
            }
            auto it1 = it;
            auto it2 = it;

            std::vector<int> face;
            face.reserve(deg.size());
            face.push_back(it1->first[0]);
            seen_edges.erase(it1->first);

            it1 = inc(it1);
            it2 = inc(it2);
            it2 = inc(it2);
            for(; it1 != it && it1 != it2; it1=inc(it1), it2=inc(inc(it2))) {
                seen_edges.erase(it1->first);
                if(myface >= 0) {
                    auto e0 = it1->first;
                    auto e1 = it1->second;
                    auto& [a,b] = e0;
                    auto& [c,d] = e1;
                    face.push_back(a);
                }
            }
            if(it1 != it) {
                assert(it1 != it2);
            }
            ret[myface] = std::move(face);
        }
        return ret;

    }

    void FaceCollapser::finalize() {
        face_ds.reduce_all();
        std::map<int,int> reindexer;

        int null_root = -1;
        if(face_ds.has_node(-1)) {
            null_root = face_ds.get_root(-1).data;
        }
        for(auto&& i: face_ds.root_indices()) {
            if(i != null_root) {
                reindexer[face_ds.node(i).data] = reindexer.size();
            }
        }
        reindexer[null_root] = -1;

        for(auto&& [e,pr]: m_edge_to_face) {
            auto&& [a,b] = e;
            auto&& [c,s] = pr;

            int root = face_ds.get_root(c).data;
            c = reindexer[root];
        }



    }
}
