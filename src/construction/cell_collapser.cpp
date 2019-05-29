#include "mandoline/construction/cell_collapser.hpp"
#include <algorithm>
#include <iterator>

#include <mtao/iterator/enumerate.hpp>
namespace mandoline::construction {

    CellCollapser::CellCollapser(const mtao::map<int,CutFace<3>>& faces) {
        for(auto&& [cid,pr]: mtao::iterator::enumerate(faces)) {
            auto&& [fid,cutface] = pr;
            auto&& CS2 = cutface.indices;
            auto&& N = cutface.N;
            CutFace<3> CS = cutface;
            std::set<Edge> my_flap_edges;
            CS.indices.clear();
            for(auto&& C2: CS2) {
                std::set<std::vector<int>> Cs;
                std::tie(Cs,my_flap_edges) = clean_unmanifold(C2);
                for(auto&& C: Cs) {
                    if(C.size() < 3) {
                        continue;
                    } else {
                        CS.indices.emplace(std::move(C));
                    }
                }
            }
            if(!CS.indices.empty()) {
                auto&& CSS = m_faces[fid] = std::move(CS);
                cell_ds.add_node(2*cid);
                cell_ds.add_node(2*cid+1);

                auto add_halfface = [&](Edge e) {
                    {
                        HalfFace hf{fid,e};
                        m_halfface_to_cell[hf] = std::make_tuple(2*cid,true);
                    }
                    std::swap(e[0],e[1]);
                    {
                        HalfFace hf{fid,e};
                        m_halfface_to_cell[hf] = std::make_tuple(2*cid+1,false);
                    }
                };
                for(auto&& C: CSS.indices) {
                    for(int i = 0; i < C.size(); ++i) {
                        add_halfface(Edge{{C[i],C[(i+1)%C.size()]}});
                    }
                }
                if(!my_flap_edges.empty()) {
                    continue;
                    for(auto&& e: my_flap_edges) {
                        auto [a,b] = e;
                        add_halfface(e);
                    }
                }
            }
        }
    }

    auto CellCollapser::clean_unmanifold(const std::vector<int>& C) -> std::tuple<std::set<std::vector<int>>,std::set<Edge>>{
        //auto t = mtao::logging::profiler("clean unmanifold",false,"profiler");
        //M is true if we want to keep a vertex
        //std::vector<bool> M(C.size(),true);
        std::set<Edge> E;
        for(int i = 0; i < C.size(); ++i) {
            Edge e{{C[i],C[(i+1)%C.size()]}};
            E.emplace(e);
        }
        std::set<Edge> flaps;
        //If an edge is duplicated we remove the ends
        for(int i = 0; i < C.size(); ++i) {
            int j = (i+1)%C.size();
            Edge e{{C[j],C[i]}};
            if(E.find(e) != E.end()) {
                if(e[0] > e[1]) {
                    flaps.insert(e);
                }
                invalid_edges.insert(e);
                std::swap(e[0],e[1]);
                if(e[0] > e[1]) {
                    flaps.insert(e);
                }
                invalid_edges.insert(e);
            }
        }
        if(invalid_edges.size() == 0) {
            return {{C},flaps};
        }
        mtao::map<int,int> next;
        for(int i = 0; i < C.size(); ++i) {
            int j = (i+1)%C.size();
            Edge e{{C[i],C[j]}};
            if(invalid_edges.find(e) == invalid_edges.end()) {
                next[C[i]] = C[j];
            }
        }
        /*
           for(int i = 0; i < C.size(); ++i) {
           int j = (i+1)%C.size();
           if(M[i] || M[j]) {
           next[C[i]] = C[j];
           }

           }
           */
        if(next.size() < C.size()) {
            std::set<std::vector<int>> ret;
            std::set<int> seen;
            for(auto&& [start, n]: next) {
                if(seen.find(start) != seen.end()) {
                    continue;
                }
                std::vector<int> NV;
                NV.push_back(start);
                seen.insert(start);
                int cur = n;
                while(seen.find(cur) == seen.end()) {
                    NV.push_back(cur);
                    seen.insert(cur);
                    cur = next[cur];

                }
                ret.insert(NV);
            }

            return {ret,flaps};
        } else {
            return {{C},flaps};
        }
    }


    std::set<int> CellCollapser::folded_faces() const {
        std::set<int> folded;
        for(auto&& [hf, pr]: m_halfface_to_cell) {
            auto&& [cell,sgn] = pr;
            if(dual_cell(hf) == cell) {
                folded.insert(std::get<0>(hf));
            }
        }
        return folded;
    }


    auto CellCollapser::collect_halffaces() const -> mtao::map<Edge,std::set<const HalfFace*>> {
        mtao::map<Edge,std::set<const HalfFace*>> ret;
        for(auto&& [hf,pr]: m_halfface_to_cell) {
            auto&& e = std::get<1>(hf);
            if(invalid_edges.find(e) == invalid_edges.end()) {
                ret[e].insert(&hf);
            }
        }
        return ret;
    }
    int CellCollapser::dual_cell(const HalfFace& hf) const {
        HalfFace h = hf;
        auto& e = std::get<1>(h);
        std::swap(e[0],e[1]);

        return cell(h);
    }
    int CellCollapser::cell(const HalfFace& hf) const {
        return std::get<0>(m_halfface_to_cell.at(hf));
    }


    std::vector<std::set<int>> CellCollapser::cell_faces() const {
        std::vector<std::set<int>> ret(m_cell_boundaries.size());
        std::transform(m_cell_boundaries.begin(),m_cell_boundaries.end(), ret.begin(), [](auto&& V) {
                std::set<int> ret;
                std::transform(V.begin(),V.end(), std::inserter(ret,ret.end()), [](auto&& v) {
                        return std::get<0>(v);
                        });

                return ret;
                });

        return ret;

    }

    void CellCollapser::remove_boundary_cells_from_faces(const std::set<int>& boundary_faces) {
        m_cell_boundaries.erase(std::remove_if(m_cell_boundaries.begin(),m_cell_boundaries.end(),[&](auto&& m) -> bool {

                    for(auto&& [i,b]: m) {
                    if(boundary_faces.find(i) == boundary_faces.end()) {
                    return false;
                    }
                    }
                    return true;
                    }), m_cell_boundaries.end());
    }
    void CellCollapser::remove_boundary_cells_from_vertices(const std::set<int>& boundary_vertices) {


        std::set<int> boundary_faces;
        for(auto&& [i,F]: m_faces) {
            auto&& N = F.N;
            auto&& Fs = F.indices;
            bool is_boundary = true;
            for(auto&& F: Fs) {
                for(auto&& v: F) {
                if(boundary_vertices.find(v) == boundary_vertices.end()) {
                    is_boundary = false;
                    break;
                }
            }
            if(!is_boundary) {
                break;
            }
        }
        if(is_boundary) {
            boundary_faces.insert(i);
        }

    }
    remove_boundary_cells_from_faces(boundary_faces);
}
}
