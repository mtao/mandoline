#pragma once
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/iterable.hpp>
#include <mtao/iterator/zip.hpp>
#include <mtao/iterator/interval.hpp>
#include "mandoline/construction/vertex_types.hpp"
#include <type_traits>

namespace mandoline::construction {

    template <int D, typename CrossingContainer>
        std::map<const Vertex<D>*, int> crossing_indexer(const CrossingContainer& C) {
            std::map<const Vertex<D>*, int> M;
            if constexpr(std::is_same_v<typename CrossingContainer::value_type,Crossing<D>>) {
                for(auto&& c: C) {
                    M[c.vertex_ptr()] = c.index;
                }
                return M;
            } else {
                for(auto&& c: C) {
                    auto m = crossing_indexer(c);
                    M.insert(m.begin(),m.end());
                }
            }
        }
    template <int D, typename CrossingContainer>
        std::map<int, const Vertex<D>*> inverse_crossing_indexer(const CrossingContainer& C) {
            std::map<int,const Vertex<D>*> M;
            if constexpr(std::is_same_v<typename CrossingContainer::value_type,Crossing<D>>) {
                for(auto&& c: C) {
                    M[c.index] = c.vertex_ptr();
                }
                return M;
            } else {
                for(auto&& c: C) {
                    auto m = crossing_indexer(c);
                    M.insert(m.begin(),m.end());
                }
            }
        }

    template <int D, typename CrossingContainer>
        void populate_crossing_indices(CrossingContainer& C) {
            int index = 0;
            for(auto&& c: C) {
                c.index = index++;
            }
        }

    template <int D>
        Vertex<D>::Vertex(const CoordType& c): coord(c), clamped_indices((1<<D)-1) {}

    template <int D>
        template <typename Derived>
        Vertex<D>::Vertex(const CoordType& c, const Eigen::MatrixBase<Derived>& q): coord(c), quot(q), clamped_indices(bs_from_quot(q)) {}


    template <int D>
        template <typename Derived>
        Vertex<D>::Vertex(const CoordType& c, const Eigen::MatrixBase<Derived>& q, const std::bitset<D> & clamped_indices): coord(c), quot(q), clamped_indices(clamped_indices) {}

    template <int D>
        template <typename Derived>
        Vertex<D> Vertex<D>::from_vertex(const Eigen::MatrixBase<Derived>& p) {
            Vertex ret{{},p};
            ret.repair();
            ret.clamped_indices = ret.bs_from_quot(ret.quot);
            return ret;
        }

    template <int D>
        template <typename Derived>
        std::bitset<D> Vertex<D>::bs_from_quot(const Eigen::MatrixBase<Derived>& q) {
            std::bitset<D> bs;
            for(int i = 0; i < D; ++i) {
                bs[i] = (q(i) == 0);
            }
            return bs;
        }



    template <int D>
        bool Vertex<D>::operator==(const Vertex& o) const {
            return quot == o.quot && coord == o.coord;
        }

    template <int D>
        bool Vertex<D>::operator!=(const Vertex& o) const {
            return quot != o.quot || coord != o.coord;
        }

    template <int D>
        bool Vertex<D>::operator<(const Vertex& o) const {
            //If on the same line we can truely compare them
            if(coord < o.coord) {
                return true;
            } else if(coord == o.coord) {
                for(int i = 0; i < D; ++i) {
                    auto& a = quot(i);
                    auto& b = o.quot(i);
                    if(a < b) {
                        return true;
                    } else if(a == b) {
                        continue;
                    } else {
                        return false;
                    }
                }
            }
            return false;
        }

    template <int D>
        bool Vertex<D>::approx(const Vertex& o) const {
            for(int i = 0; i < D; ++i) {
                auto& q = quot(i);
                auto& c = coord[i];
                auto& oq = o.quot(i);
                auto& oc = o.coord[i];
                if(q < .5) {
                    if(oc == c) {
                        if(std::abs(q-oq) > threshold_epsilon) {
                            return false;
                        }
                    } else if(oc +1 == c) {
                        if(std::abs(q-(oq+1)) > threshold_epsilon) {
                            return false;
                        }
                    } else {
                        return false;
                    }
                } else {
                    if(oc == c) {
                        if(std::abs(q-oq) > threshold_epsilon) {
                            return false;
                        }
                    } else if(oc-1 == c) {
                        if(std::abs(oq-(q+1)) > threshold_epsilon) {
                            return false;
                        }
                    } else {
                        return false;
                    }
                }
            }
            return true;
        }



    template <int D>
        auto Vertex<D>::p() const -> Vec {
            return Eigen::Map<const mtao::Vector<int,D>>(coord.begin()).template cast<double>() + quot; 
        }

    template <int D>
        auto Vertex<D>::optional_index() const -> OptInd {
            OptInd ret;
            for(int i = 0; i < D; ++i) {
                if(clamped_indices[i]) {
                    ret[i] = coord[i];
                }
            }
            return ret;
        }

    template <int D>
        auto Vertex<D>::mask() const -> MaskType {
            MaskType m;
            for(int i = 0; i < D; ++i) {
                if(quot(i) == 0) {
                    m[i] = coord[i];
                }
            }
            return m;
        }

    template <int D>
        Vertex<D>::operator std::string() const {
            std::stringstream ss;
            ss << "(" <<mtao::eigen::stl2eigen(coord).transpose() << "{" << quot.transpose() << "}"<< "[";
            for(int i = 0; i < clamped_indices.size(); ++i) {
                ss << clamped_indices[i];
            }
            ss << "])";
            return ss.str();
        }



    template <int D>
        bool Vertex<D>::clamped(int index) const {
            return clamped_indices[index];
        }

    template <int D>
        size_t Vertex<D>::clamped_count() const { 
            return clamped_indices.count(); 
        }

    template <int D>
        bool Vertex<D>::is_grid_vertex() const {
            return clamped_indices.all(); 
        }

    template <int D>
        auto Vertex<D>::possible_cells() const -> std::set<CoordType> {
            std::set<CoordType> ret;
            for(int i = 0; i < (2 << D); ++i) {
                std::bitset<D> bs(i);
                if((bs&clamped_indices) == bs) {

                    CoordType cc = coord;
                    for(int j = 0; j < D; ++j) {
                        cc[j] -= bs[j]?1:0;
                    }
                    ret.emplace(cc);

                }
            }
            return ret;
        }



    template <int D>
        void Vertex<D>::repair() {
            auto rcoord = mtao::eigen::stl2eigen(coord);
            Vec qfloor = quot.array().floor();
            rcoord += qfloor.template cast<int>();
            quot = quot - qfloor;

        }

    template <int D>
        void Vertex<D>::apply_thresholding() {
            apply_thresholding(threshold_epsilon);
        }

    template <int D>
        void Vertex<D>::apply_thresholding(double thresh) {
            for(int i = 0; i < D; ++i) {
                auto& q = quot(i);
                auto& c = coord[i];
                if(q < .5) {
                    if(q < thresh) {
                        clamped_indices[i] = true;
                        q = 0;
                    }
                } else {
                    if(q > 1-thresh) {
                        clamped_indices[i] = true;
                        q = 0;
                        c++;
                    }
                }
            }
        }


    template <int D>
        auto Vertex<D>::operator+(const Vertex& o) const -> Vertex<D> {
            using namespace mtao::eigen;
            Vertex ret;
            auto mcm = stl2eigen(coord);
            auto ocm = stl2eigen(o.coord);
            auto rcm = stl2eigen(ret.coord);
            rcm = mcm + ocm;
            ret.clamped_indices = clamped_indices & o.clamped_indices;
            //ret.clamped_indices = 0;
            ret.quot = o.quot + quot;
            ret.repair();
            return ret;
        }

    template <int D>
        auto Vertex<D>::operator*(double val) const -> Vertex<D> {
            using namespace mtao::eigen;
            auto mcm = stl2eigen(coord);
            Vertex ret = from_vertex(val * mcm.template cast<double>());
            ret.quot += val * quot;
            ret.repair();
            ret.clamped_indices = 0;
            return ret;
        }

    template <int D>
        auto Vertex<D>::lerp(const Vertex& other, double t) const -> Vertex<D> {

            Vertex ret = (*this) * (1-t) + other * t;
            for(int i = 0; i < D; ++i) {
                if(coord[i] == other.coord[i] && clamped_indices[i] && other.clamped_indices[i]) {
                    ret.clamped_indices[i] = true;
                }
            }

            return ret;
        }


    template <int D>
        EdgeIntersection<D>::EdgeIntersection(const VType& v, double coord, int index): VType(v), edge_coord(coord), edge_index(index) {
            assert(coord > 0 && coord < 1);
        }

    template <int D>
        bool EdgeIntersection<D>::operator<(const EdgeIntersection& o) const {
            //If on the same line we can truely compare them
            if(edge_coord >= 0 && edge_index == o.edge_index) {
                return edge_coord < o.edge_coord;
            } else {
                //Thsi is so these intersections can belong on an edge
                return edge_index < o.edge_index 
                    && VType::operator<(o);
            }
        }

    template <int D>
        TriangleIntersection<D>::TriangleIntersection(const VType& is, const mtao::Vec3d& v, int triangle_index): VType(is), bary_coord(v), triangle_index(triangle_index) {
            assert(bary_coord.minCoeff() > 0);
        }


    template <int D>
        bool Crossing<D>::operator<(const Crossing& o) const {
            if(vv.index() != o.vv.index()) {
                return vv.index() < o.vv.index();
            } else {
                return std::visit([&](auto&& v) {
                        using T = std::decay_t<decltype(v)>;
                        if constexpr(std::is_same_v<T,VType const*>) {
                        return  *v < *std::get<const VType*>(o.vv);
                        } else if constexpr(std::is_same_v<T,EdgeIsect const*>) {
                        return  *v < *std::get<const EdgeIsect*>(o.vv);
                        } else if constexpr(std::is_same_v<T,TriangleIsect const*>) {
                        return  *v < *std::get<const TriangleIsect*>(o.vv);
                        }
                        },vv);
            }

        }

    template <int D>
        bool Crossing<D>::is_grid_vertex() const {
            return std::holds_alternative<VPtr>(vv);
        }

    template <int D>
        bool Crossing<D>::is_edge_intersection() const {
            return std::holds_alternative<EPtr>(vv);
        }

    template <int D>
        bool Crossing<D>::is_triangle_intersection() const {
            return std::holds_alternative<FPtr>(vv);
        }

    template <int D>
        auto Crossing<D>::vertex_ptr() const -> const VType* {
            return std::visit([](auto&& val) -> const VType*{
                    return static_cast<const VType*>(val);
                    }, vv);
        }

    template <int D>
        auto Crossing<D>::vertex() const -> const VType&{
            return *vertex_ptr();
        }

    template <int D>
        auto Crossing<D>::mask() const -> typename VType::MaskType { 
            return vertex().mask(); 
        }

    template <int D>
        auto Crossing<D>::point() const -> Vec {
            return vertex().p();
        }

    template <int D>
        Crossing<D>::operator std::string() const {
            std::stringstream ss;
            std::visit([&](auto&& v) {
                    using T = std::decay_t<decltype(v)>;
                    if constexpr(std::is_same_v<T,VType const*>) {
                    ss << "V";
                    } else if constexpr(std::is_same_v<T,EdgeIsect const*>) {
                    ss << "E" << "{" << v->edge_index << "}";
                    } else if constexpr(std::is_same_v<T,TriangleIsect const*>) {
                    ss << "F" << "{" << v->triangle_index<< "}";
                    } else {
                    ss << "???";
                    }
                    },vv);
            if(index >= 0) {
                ss << "[" << index << "]" << std::string(vertex());
            } else {
                ss << "[?]" << std::string(vertex());
            }
            return ss.str();
        }

    template <int D>
        EdgeCrossing<D>::operator std::string() const {
            std::stringstream ss;
            ss << "[" << index << "_" << index << ":" << edge_coord() <<"]" << std::string(vertex());
            return ss.str();
        }

    template <int D>
        int EdgeCrossing<D>::find_unbound(const VType& v) {
            //We are not allowed to have grid vertices
            assert(v.clamped_indices.count() == D-1);
            for(int i = 0; i < D; ++i) {
                if(!v.clamped(i)) {
                    return i;
                }
            }
        }

    template <int D>
        EdgeCrossing<D>::EdgeCrossing(const VType* v, int i): val(v), axis(find_unbound(*v)), index(i) {}

    template <int D>
        EdgeCrossing<D>::EdgeCrossing(const VType& v, int i): EdgeCrossing(&v,i) {}

    template <int D>
        EdgeCrossing<D>::EdgeCrossing(const Crossing<D>& v, int i): EdgeCrossing(&v.vertex(),i) {}

    template <int D>
        EdgeCrossing<D>::EdgeCrossing(const Crossing<D>& v): EdgeCrossing(&v.vertex(),v.index) {}

    template <int D>
        bool EdgeCrossing<D>::operator<(const EdgeCrossing& o) const {
            ////If on the same line we can truely compare them, by their position along the line
            //assert(index == o.index);
            return edge_coord() < o.edge_coord();
        }

    template <int D>
        double EdgeCrossing<D>::edge_coord() const {
            return vertex().quot(axis);
        }

    template <int D>
        auto EdgeCrossing<D>::vertex() const -> const VType& {
            return *val;
        }

    template <int D>
        auto EdgeCrossing<D>::mask() const -> typename VType::MaskType {
            return vertex().mask(); 
        }

    template <int D>
        auto EdgeCrossing<D>::point() const -> Vec {
            return vertex().p();
        }

}

