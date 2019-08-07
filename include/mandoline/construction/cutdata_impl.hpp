#pragma once
#include "mandoline/construction/cutdata.hpp"
#include <map>
#include <mtao/iterator/enumerate.hpp>
#include <set>
#include <tuple>
#include <iostream>
#include <thread>

namespace mandoline::construction {

    template <int D, typename Indexer>
        CutData<D,Indexer>::CutData(const Indexer& indexer, const mtao::vector<VType>& V, const Faces& F): Indexer(indexer), m_V(V)  {
            set_topology(F);
        }
    template <int D, typename Indexer>
        CutData<D,Indexer>::CutData(const Indexer& indexer, const mtao::vector<VType>& V, const Edges& E,const Faces& F,  const Faces& FEA): Indexer(indexer), m_V(V)  {
            set_topology(E,F,FEA);

        }
    template <int D, typename Indexer>
        void CutData<D,Indexer>::set_topology(const Faces& F) {
            auto E = mtao::geometry::mesh::boundary_facets(F);
            set_topology(E,F);
        }

    template <int D, typename Indexer>
        void CutData<D,Indexer>::clean_edges() {
            std::vector<std::set<std::tuple<double,size_t>>> edges(nE());
            //for(auto&& [idx,eisect]: mtao::iterator::enumerate(m_edge_intersections)) {
                //edges[eisect.edge_index].emplace(eisect.edge_coord,idx);
                //}
            for(auto&& eisect: m_edge_intersections) {
                eisect.clear();
            }
            /*
            m_edge_intersections.clear();
            m_edge_intersections.reserve(m_E.size());
            for(int i = 0; i < m_E.cols(); ++i) {
                m_edge_intersections.emplace_back(m_V,m_E,i);
            }
            */
        }
    template <int D, typename Indexer>
        void CutData<D,Indexer>::clean_triangles() {
            for(auto&& tisect: m_triangle_intersections) {
                tisect.clear();
            }

            /*
            m_triangle_intersections.clear();
            m_triangle_intersections.reserve(m_F.size());
            for(int i = 0; i < m_F.cols(); ++i) {
                m_triangle_intersections.emplace_back(m_V,m_F,m_E,m_FE, m_edge_intersections,i);
            }
            */
        }


    template <int D, typename Indexer>
        auto CutData<D,Indexer>::F_asCrossings() const -> Faces {
            return m_F.unaryExpr([&](int index) -> int {
                    return m_crossings[index].index;
                    });
        }

    template <int D, typename Indexer>
        std::vector<std::vector<int>> CutData<D,Indexer>::faces() const {
            std::vector<std::vector<int>> F(m_cut_faces.size());
            std::transform(m_cut_faces.begin(),m_cut_faces.end(),F.begin(),[](const CutMeshFace<D>& cf) {
                    return cf.indices;
                    });
            return F;
        }



    template <int D, typename Indexer>
        std::map<std::vector<int>,int> CutData<D,Indexer>::triangle_index_ownership() const {
            return {};
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edge_index_ownership() const -> std::map<Edge,int>{
            return {};
        }



    template <int D, typename Indexer>
        void CutData<D,Indexer>::bake(const std::optional<SGType>& grid, bool fuse) {
            auto t = mtao::logging::profiler("mesh bake",false,"profiler");
            {
                auto t = mtao::logging::profiler("mesh vertex bake",false,"profiler");
                {

                    //auto t = mtao::logging::timer("data bake edges");
                    auto t = mtao::logging::profiler("mesh edge intersections",false,"profiler");
                    int i;
#pragma omp parallel for
                    for (i=0; i<m_edge_intersections.size(); i++) {
                        m_edge_intersections[i].bake(grid);
                    }
                }
                {
                    auto t = mtao::logging::profiler("mesh face intersections",false,"profiler");
                    //auto t = mtao::logging::timer("data bake faces");
                    int i;
#pragma omp parallel for
                    for (i=0; i<m_triangle_intersections.size(); i++) {
                        m_triangle_intersections[i].bake(grid);
                    }
                    //#pragma omp parallel
                }
                {
                    //auto t = mtao::logging::timer("Baking crossings");
                    auto t = mtao::logging::profiler("mesh vertex unification",false,"profiler");
                    m_crossings = compute_crossings();
                    m_vertex_indexer = gv_index_map(m_crossings);
                }
            }

            using namespace mtao::geometry::mesh;
            using Edge = std::array<int,2>;

            std::map<int,int> cut_to_primal;
            int cut_face_index = 0;

            {
                auto tt = mtao::logging::profiler("mesh face bake",false,"profiler");
                std::mutex face_mutex;
                auto add_face = [&](const coord_mask<D>& m, const std::vector<int>& F, int fidx) {
                    std::scoped_lock lock(face_mutex);
                    cut_to_primal[m_cut_faces.size()] = fidx;
                    m_cut_faces.emplace_back(m,F,fidx);

                };

                int i = 0;
#pragma omp parllel for
                for (i=0; i<m_triangle_intersections.size(); i++) {
                    auto&& FI = m_triangle_intersections[i];
                    coord_mask<D> fi_mask = FI.mask();
                    bool can_skip = true;
                    if(FI.intersections.empty()) {
                        for(auto&& e: FI.edge_isects) {
                            if(!e->intersections.empty()) {
                                can_skip = false;
                                break;
                            }
                        }

                        if(can_skip) {
                            std::vector<int> f(3);
                            auto fi = F(FI.triangle_index);
                            for(auto&& [i,v]: mtao::iterator::enumerate(f)) {
                                v = m_crossings[fi(i)].index;
                            }
                            //std::transform(FI.vptr_tri.begin(),FI.vptr_tri.end(),f.begin(),[&](auto&& v) {
                            //        return index(v);
                            //        });
                            add_face(fi_mask,f,FI.triangle_index);
                        }

                    }  else {
                        can_skip = false;
                    }
                    if(!can_skip){
                        auto Fs = FI.faces(m_vertex_indexer);
                        for(auto&& F: Fs) {
                            add_face(fi_mask,F,FI.triangle_index);
                        }
                    }
                    //std::cout << "Face intersections: " << FI.face_index << std::endl;

                }
            }

        }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::get_bary(int face_index, const Crossing<D>& crossing) const -> mtao::Vec3d{
            auto& FI = m_triangle_intersections[face_index];
            return std::visit(
                    [&](auto&& v) -> mtao::Vec3d {
                    using T = typename std::decay_t<decltype(v)>;
                    if constexpr(std::is_same_v<T,VType const*>) {
                    //try to find it on the vertices
                    for(auto&& [i,ptr]: mtao::iterator::enumerate(FI.vptr_tri)) {
                    if(ptr == v) {
                    return mtao::Vec3d::Unit(i);
                    }
                    }
                    } else if constexpr(std::is_same_v<T,EdgeIsect const*>) {
                    auto&& EI = m_edge_intersections[v->edge_index];
                    return FI.edge_bary(EI,v->edge_coord);
                    } else if constexpr(std::is_same_v<T,TriIsect const*>) {
                    return v->bary_coord;
                    }
                    return FI.get_bary(*v);
                    }
                    ,crossing.vv);
        }


    template <int D, typename Indexer>
        auto CutData<D,Indexer>::cut_vertices() const -> ColVecs {
            int gsize = grid_size();
            int size = cut_vertex_size();
            if(size == gsize) {
                return {};
            }
            ColVecs ret(D,size - gsize);
            for(auto&& c: crossings()) {
                if(c.index >= gsize) {
                    ret.col(c.index - gsize) = c.point();
                }
            }
            return ret;
        }



    //this includes the grid vertices
    template <int D, typename Indexer>
        int CutData<D,Indexer>::cut_vertex_size() const {
            int size = grid_size();
            for(auto&& c: crossings()) {
                size = std::max<int>(size,c.index+1);
            }
            return size;
        }
    template <int D, typename Indexer>
        void CutData<D,Indexer>::set_topology(const Edges& E,const Faces& F,  const Faces& FEA) {
            m_E = E;
            m_F = F;
            m_FE = FEA;
            reset_intersections();
        }

    template <int D, typename Indexer>
        void CutData<D,Indexer>::update_vertices(const mtao::vector<VType>& V) {
            assert(V.size() == m_V.size());
            std::copy(V.begin(),V.end(),m_V.begin());
            clear();
        }
    template <int D, typename Indexer>
        void CutData<D,Indexer>::update_grid(const Indexer& I) {
            Indexer::operator=(I);
        }

    template <int D, typename Indexer>
        void CutData<D,Indexer>::clear() {
            m_crossings.clear();
            m_vertex_indexer.clear();
            m_cut_faces.clear();
            clean_edges();
            clean_triangles();
        }
    template <int D, typename Indexer>
        void CutData<D,Indexer>::reset_intersections() {
            m_edge_intersections.clear();
            m_triangle_intersections.clear();
            m_edge_intersections.reserve(m_E.size());
            m_triangle_intersections.reserve(m_F.size());
            for(int i = 0; i < m_E.cols(); ++i) {
                m_edge_intersections.emplace_back(m_V,m_E,i);
            }

            if(m_FE.cols() == 0) {
                m_FE = mtao::geometry::mesh::boundary_elements(m_F,m_E);
            }
            for(int i = 0; i < m_F.cols(); ++i) {
                m_triangle_intersections.emplace_back(m_V,m_F,m_E,m_FE, m_edge_intersections,i);
            }

        }

    template <int D, typename Indexer>
        void CutData<D,Indexer>::add_edge_isect(int idx, double t) {
            auto eis = m_edge_intersections[idx];
            eis.intersections.emplace_back(eis.from_coord(t));
        }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::flat_edge_intersections() const -> std::vector<const EdgeIntersection<D>*>{
            std::vector<const EdgeIntersection<D>*> ret;
            for(auto&& eis: m_edge_intersections) {
                std::transform(eis.intersections.begin(),eis.intersections.end(), std::back_inserter(ret), [](const EdgeIntersection<D>& v) {
                        return &v;
                        });
            }
            return ret;
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::flat_triangle_intersections() const -> std::vector<const TriangleIntersection<D>*>{
            std::vector<const TriangleIntersection<D>*> ret;
            for(auto&& fis: m_triangle_intersections) {
                std::transform(fis.intersections.begin(),fis.intersections.end(), std::back_inserter(ret), [](const TriangleIntersection<D>& v) {
                        return &v;
                        });
            }
            return ret;
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::crossings() const -> const std::vector<Crossing<D>>& { return m_crossings; }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::compute_crossings(bool fuse) const -> std::vector<Crossing<D>>{

                auto t2 = mtao::logging::timer("creating crossings");
            auto V = vertex_crossings();
            auto E = edge_crossings();
            auto F = face_crossings();
            std::vector<Crossing<D>> ret;
            ret.resize(V.size() + E.size() + F.size());

            std::map<VType, int> index_map;
            std::vector<VType> gvs;

            int eoff = V.size();
            int foff = V.size() + E.size();




                {
                auto t = mtao::logging::timer("copying crossings");
                std::copy(V.begin(),V.end(),ret.begin());
                std::copy(E.begin(),E.end(),ret.begin()+eoff);
                std::copy(F.begin(),F.end(),ret.begin()+foff);
                }


                int count = grid_size();
                auto t = mtao::logging::timer("indexing crossings");

                if(fuse) {
                    int i;
                    //for(auto&& [i,v]: mtao::iterator::enumerate(V.size()+E.size(),F)) {
                    for(i = 0; i < ret.size(); ++i) {
                        auto& p = ret[i];
                        auto& gv = p.vertex();
                        if(gv.is_grid_vertex()) {
                            p.index = index_map[gv] = grid_index(gv);
                        } else {
                            {
                                if(auto it = index_map.find(gv); it != index_map.end()) {
                                    p.index = it->second;
                                } else {
                                    p.index = index_map[gv] = count++;
                                }
                            }

                        }
                    }
                } else {
                    int i;
                    for(i = 0; i < ret.size(); ++i) {
                        auto& p = ret[i];
                        p.index = i;
                    }
                }
            return ret;

        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::vertex_crossings() const -> std::vector<Crossing<D>> {
            std::vector<Crossing<D>> ret;
            auto e = mtao::iterator::enumerate(m_V);
            int grid_size = this->grid_size();
            std::transform(e.begin(),e.end(), std::back_inserter(ret), 
                    [grid_size](auto&& pr) {
                    auto&& [i,iptr] = pr;
                    return Crossing<D>(&iptr,i+grid_size);
                    });
            return ret;
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edge_crossings() const -> std::vector<Crossing<D>> {
            return facet_crossings(flat_edge_intersections());
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::face_crossings() const -> std::vector<Crossing<D>>{
            return facet_crossings(flat_triangle_intersections());
        }
    template <int D, typename Indexer>
        template <typename IsectType>
        auto CutData<D,Indexer>::facet_crossings(const std::vector<const IsectType*>& isects) const -> std::vector<Crossing<D>>{
            std::vector<Crossing<D>> ret;
            auto e = mtao::iterator::enumerate(isects);
            std::transform(e.begin(),e.end(), std::back_inserter(ret), 
                    [](auto&& pr) {
                    auto&& [i,iptr] = pr;
                    return Crossing<D>(iptr,i);
                    });
            return ret;
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::gv_index_map(const std::vector<Crossing<D>>& crossings) const -> std::map<const VType*, int>{
            std::map<const VType*, int> ret;
            for(auto&& c: crossings) {
                ret[&c.vertex()] = c.index;
            }
            return ret;
        }

    template <int D, typename Indexer>
        template <int U>
        auto CutData<D,Indexer>::facets(const std::map<const VType*,int>& gv_idx_map, const std::vector<std::array<const VType*,U>>& facets) const -> mtao::ColVectors<int,U> {
            auto E = mtao::eigen::stl2eigen(facets);
            return E.unaryExpr([&](const VType* ptr) -> int {
                    return gv_idx_map.at(ptr);
                    });

        }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::gvedge2edges(const std::vector<VPtrEdge>& gvedges, std::map<const VType*,int>& indexer) const -> std::vector<Edge> {
            auto GVE = gvedges;
            std::vector<Edge> E;
            E.reserve(GVE);
            for(auto&& gve: GVE) {
                Edge e;
                std::transform(gve.begin(),gve.end(),e.begin(),[&](const VType* p) {
                        return indexer.at(p);
                        });
                if(e[0] != e[1]) {
                    E.push_back(e);
                }
            }
            return E;
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edge_edges() const -> mtao::ColVectors<int,2>{
            return mtao::eigen::stl2eigen(edge_edges(m_vertex_indexer));
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edge_indices() const -> std::map<VPtrEdge,int> {
            std::map<VPtrEdge,int> ret;
            for(auto&& e: m_edge_intersections) {
                for(auto gve: e.vptr_edges()) {
                    ret[gve] = e.edge_index;
                    std::swap(gve[0],gve[1]);
                    ret[gve] = e.edge_index;
                }
            }
            return ret;
        }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edge_edges(const std::map<const VType*,int>& gv_idx_map) const -> std::set<Edge>{
            auto t = mtao::logging::timer("data pulling edge edges");
            std::set<Edge> edges;



            std::vector<std::set<Edge>> per_edges(m_edge_intersections.size());
            int i = 0;
#ifdef MTAO_OPENMP
#pragma omp parallel for
#endif//MTAO_OPENMP
            for(i=0; i < m_edge_intersections.size(); ++i) {
                auto gv = m_edge_intersections[i].vptr_edges();
                auto&& e = per_edges[i];
                for(auto&& ve: gv) {
                    Edge e{{index(ve[0]),index(ve[1])}};
                    if(e[0] != e[1]) {
                        std::sort(e.begin(),e.end());
                        per_edges[i].insert(e);

                    }
                }
            }

            for(auto&& e:  per_edges) {
                edges.insert(e.begin(),e.end());
            }
            return edges;
        }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edge_ownership() const -> std::map<VPtrEdge, std::set<VPtrEdge>>{
            std::map<VPtrEdge, std::set<VPtrEdge>> ret;
            for(auto&& ei: m_edge_intersections) {
                auto gv = ei.vptr_edges();
                auto& edges = ret[ei.vptr_edge];
                edges.insert(edges.end(),gv.begin(),gv.end());
            }
            return ret;
        }

    template <int D, typename Indexer>
        auto CutData<D,Indexer>::face_edges() const -> std::set<Edge>{
            return face_edges(m_vertex_indexer);
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::face_edges(const std::map<const VType*,int>& gv_idx_map) const -> std::set<Edge> {
            auto t = mtao::logging::timer("data pulling face edges");
            std::set<Edge> edges;

            int i = 0;
#ifdef MTAO_OPENMP
#pragma omp parallel for
#endif//MTAO_OPENMP
            for(i=0; i < m_triangle_intersections.size(); ++i) {
                auto& FT = m_triangle_intersections[i];
                //if(!FT.is_cut()) {
                //    continue;
                //}
                auto gv = FT.vptr_edges();
                for(auto&& ve: gv) {
                    Edge e{{index(ve[0]),index(ve[1])}};
                    if(e[0] != e[1]) {
                        std::sort(e.begin(),e.end());
#ifdef MTAO_OPENMP
#pragma omp critical
#endif//MTAO_OPENMP
                        edges.insert(e);

                    }
                }
            }



            return edges;
        }


    template <int D, typename Indexer>
        auto CutData<D,Indexer>::stl_edges(const std::map<const VType*,int>& gv_idx_map) const -> std::set<Edge> {
            auto t = mtao::logging::timer("data pulling edges");
            std::set<Edge> EE = edge_edges(gv_idx_map), EF = face_edges(gv_idx_map);
            EE.insert(EF.begin(),EF.end());
            return EE;
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edges(const std::map<const VType*,int>& gv_idx_map) const -> mtao::ColVectors<int,2> {
            return mtao::eigen::stl2eigen(stl_edges(gv_idx_map));
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::edges() const -> mtao::ColVectors<int,2>{
            return edges(m_vertex_indexer);
        }
    template <int D, typename Indexer>
        auto CutData<D,Indexer>::stl_edges() const -> std::set<Edge>{
            return stl_edges(m_vertex_indexer);
        }

    template <int D, typename Indexer>

        Eigen::SparseMatrix<double> CutData<D,Indexer>::barycentric_map() const {

            Eigen::SparseMatrix<double> M(cut_vertex_size(),nV());
            std::vector<Eigen::Triplet<double>> trips;
            auto CV = cut_vertices();
            for(auto&& c: crossings()) {              
                int row = c.index;
                std::visit(
                        [&](auto&& v) {
                        using T = typename std::decay_t<decltype(v)>;
                        if constexpr(std::is_same_v<T,VType const*>) {
                        int col = std::distance(&m_V[0],v);
                            trips.emplace_back(row,col,1);
                        } else if constexpr(std::is_same_v<T,EdgeIsect const*>) {
                            auto e = E(v->edge_index);
                            double t = v->edge_coord;
                            trips.emplace_back(row,e(0),1-t);
                            trips.emplace_back(row,e(1),t);
                        } else if constexpr(std::is_same_v<T,TriIsect const*>) {
                            auto f = F(v->triangle_index);
                            auto&& bc = v->bary_coord;
                            for(int i = 0; i < 3; ++i) {
                                trips.emplace_back(row,f(i),bc(i));
                            }
                        }
                    }
                    ,c.vv);
            }
            M.setFromTriplets(trips.begin(),trips.end());
        return M;
    }
}
