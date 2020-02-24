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


            // merge edge face ids 
            template <typename Derived>
                void unify_boundary_loops(const Eigen::MatrixBase<Derived>& V);
            // bake the face structure. for now it is only merge, but could be more?
            template <typename Derived>
                void bake(const Eigen::MatrixBase<Derived>& V, bool nonsimple_faces=true);



            template <typename Derived>
                bool is_inside(const Eigen::MatrixBase<Derived>& V, const std::vector<int>& a, const std::vector<int>& b) const;

            template <typename Derived>
                void merge_faces(const Eigen::MatrixBase<Derived>& V);
            template <typename Derived>
                std::map<int,typename Derived::Scalar> volumes(const Eigen::MatrixBase<Derived>& V) const;

            // cleans up the face identities to start from 0
            void finalize();

            // returns faces with the assumption that there are no holes.
            // when there are holes, this code will generate every boundary loop in the domain, some of which may have positive or negative volume.
            std::map<int,std::vector<int>> faces_no_holes() const;
            std::map<int,std::set<std::vector<int>>> faces() const;
            template <typename Derived>
            void faces(const Eigen::MatrixBase<Derived>& V) const;
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
        void FaceCollapser::bake(const Eigen::MatrixBase<Derived>& V, bool nonsimple_faces) {
            unify_boundary_loops(V);

            if(nonsimple_faces) {
                merge_faces(V);
            }

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

    template <typename Derived>
    std::map<int,std::vector<std::vector<int>>> FaceCollapser::faces(const Eigen::MatrixBase<Derived>& V) const {
        std::map<int,std::vector<std::vector<int>>> ret;

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
    template <typename Derived>
        void FaceCollapser::merge_faces(const Eigen::MatrixBase<Derived>& V) {
            auto vols = volumes(V);
            std::map<int,std::set<int>> halfedge_partial_ordering;
            auto faces = faces();
            std::set<int> outer_hes;
            std::set<int> interior_hes;
            mtao::data_structures::DisjointSet<int> ds;
            ds.add_node(-1);
            /*
            for(auto&& he: cell_halfedges) {
                auto area = signed_area(V,edge(he));
                if(area  > 0) {
                    outer_hes.insert(he);
                } else if(area  < 0) {
                    interior_hes.insert(he);
                } else {

                    //assert(area != 0);
                    int ci = cell_index(he);
                    ds.add_node(ci);
                    ds.join(-1,ci);

                }
            }

            //the vertices that make each cell
            std::map<int,std::vector<int>> CM;
            std::map<int,std::set<int>> CMs;
            for(auto&& he: cell_halfedges) {
                auto&& v = CM[he] = cell_he(he);
                CMs[he] = std::set<int>(v.begin(),v.end());
            }

            for(auto&& he: cell_halfedges) {
                int cell = cell_index(he);
                ds.add_node(cell);
                ds.add_node(cell_index(dual_index(he)));
            }
            for(auto&& ohe: outer_hes) {
                int ocell = cell_index(ohe);
                //int docell = cell_index(dual_index(ohe));
                for(auto&& ihe: interior_hes) {
                    int icell = cell_index(ihe);
                    int dicell = cell_index(dual_index(ihe));
                    if(ocell == icell || dicell == ocell) {
                        continue;
                    } else {
                        if( is_inside(V,edge(ohe),CM[ihe])) {
                            auto&& potential_parent = CMs[ohe];
                            auto&& potential_child = CMs[ihe];
                            //if there's some random edge connecting these two we need to find it
                            if(std::includes(potential_parent.begin(),potential_parent.end(),potential_child.begin(),potential_child.end())) {
                                auto&& child = CM[ihe];
                                auto&& parent = CM[ohe];
                                std::set<std::array<int,2>> edges;
                                for(int i = 0; i < parent.size(); ++i) {
                                    std::array<int,2> e;
                                    e[0] = parent[i];
                                    e[1] = parent[(i+1)%parent.size()];
                                    std::sort(e.begin(),e.end());
                                    edges.emplace(e);
                                }
                                for(int i = 0; i < child.size(); ++i) {
                                    std::array<int,2> e;
                                    e[0] = child[i];
                                    e[1] = child[(i+1)%child.size()];
                                    std::sort(e.begin(),e.end());
                                    if(edges.find(e) == edges.end()) {
                                        auto v = (V.col(e[0]) + V.col(e[1]))/2;
                                        if(is_inside(V,edge(ohe),v)) {
                                            halfedge_partial_ordering[ohe].insert(ihe);
                                            break;
                                        }
                                    }
                                }
                            } else {
                                halfedge_partial_ordering[ohe].insert(ihe);
                            }
                        }
                    }
                }

            }

            std::map<int,std::set<int>> ordered_partial_ordering;
            std::set<int> seen_cells;
            for(auto&& [h, hs]: halfedge_partial_ordering) {
                ordered_partial_ordering[hs.size()].insert(h);
            }
            for(auto&& pr: ordered_partial_ordering) {
                auto&& hs = std::get<1>(pr);
                for(auto&& h: hs) {
                    int cell = cell_index(h);
                    auto&& children = halfedge_partial_ordering[h];

                    for(auto&& che: children) {
                        if(seen_cells.find(che) == seen_cells.end()) {
                            int ci = cell_index(che);
                            ds.join(cell,ci);
                            seen_cells.emplace(che);
                        }
                    }
                }
            }
            ds.reduce_all();
            {

                auto ci = cell_indices();
                for(int i = 0; i < size(); ++i) {
                    int c = ci(i);
                    if(ds.has_node(c)) {
                        int root = ds.get_root(ci(i)).data;
                        cell_indices()(i) = root;
                    }
                }
            }
            */
        }
    template <typename Derived>
        std::map<int,typename Derived::Scalar> FaceCollapser::volumes(const Eigen::MatrixBase<Derived>& V) const {
        }

    template <typename Derived>
        std::map<int,typename Derived::Scalar> FaceCollapser::volumes(const Eigen::MatrixBase<Derived>& V) const {
        }
    template <typename Derived>
        bool FaceCollapser::is_inside(const Eigen::MatrixBase<Derived>& V, const std::vector<int>& a, const std::vector<int>& b) const {
        }
}
