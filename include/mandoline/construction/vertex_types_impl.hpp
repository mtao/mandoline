#pragma once
#include <balsa/eigen/stl2eigen.hpp>
#include <mtao/iterator/zip.hpp>
#include <mtao/iterator/interval.hpp>
#include "mandoline/construction/vertex_types.hpp"
#include <type_traits>

namespace mandoline::construction {

template<int D, typename CrossingContainer>
std::map<const Vertex<D> *, int> crossing_indexer(const CrossingContainer &C) {
    std::map<const Vertex<D> *, int> M;
    if constexpr (std::is_same_v<typename CrossingContainer::value_type, Crossing<D>>) {
        for (auto &&c : C) {
            M[c.vertex_ptr()] = c.index;
        }
        return M;
    } else {
        for (auto &&c : C) {
            auto m = crossing_indexer(c);
            M.insert(m.begin(), m.end());
        }
    }
}
template<int D>
std::map<int, const Vertex<D> *> inverse_crossing_indexer_from_indexer(const std::map<const Vertex<D> *, int>  &C) {
    std::map<int, const Vertex<D> *> M;
    for(auto&& [a,b]: C) {
        M.emplace(b,a);
    }
    return M;
}
template<int D, typename CrossingContainer>
std::map<int, const Vertex<D> *> inverse_crossing_indexer(const CrossingContainer &C) {
    std::map<int, const Vertex<D> *> M;
    if constexpr (std::is_same_v<typename CrossingContainer::value_type, Crossing<D>>) {
        for (auto &&c : C) {
            M[c.index] = c.vertex_ptr();
        }
        return M;
    } else {
        for (auto &&c : C) {
            auto m = crossing_indexer(c);
            M.insert(m.begin(), m.end());
        }
    }
}

template<int D, typename CrossingContainer>
void populate_crossing_indices(CrossingContainer &C) {
    int index = 0;
    for (auto &&c : C) {
        c.index = index++;
    }
}


template<int D>
EdgeIntersection<D>::EdgeIntersection(const VType &v, double coord, int index) : VType(v), edge_coord(coord), edge_index(index) {
    assert(coord > 0 && coord < 1);
}

template<int D>
bool EdgeIntersection<D>::operator<(const EdgeIntersection &o) const {
    //If on the same line we can truely compare them
    if (edge_coord >= 0 && edge_index == o.edge_index) {
        return edge_coord < o.edge_coord;
    } else {
        //Thsi is so these intersections can belong on an edge
        return edge_index < o.edge_index
               && VType::operator<(o);
    }
}

template<int D>
TriangleIntersection<D>::TriangleIntersection(const VType &is, const balsa::eigen::Vec3d &v, int triangle_index) : VType(is), bary_coord(v), triangle_index(triangle_index) {
    assert(bary_coord.minCoeff() > 0);
}


template<int D>
bool Crossing<D>::operator<(const Crossing &o) const {
    if (vv.index() != o.vv.index()) {
        return vv.index() < o.vv.index();
    } else {
        return std::visit([&](auto &&v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, VType const *>) {
                return *v < *std::get<const VType *>(o.vv);
            } else if constexpr (std::is_same_v<T, EdgeIsect const *>) {
                return *v < *std::get<const EdgeIsect *>(o.vv);
            } else if constexpr (std::is_same_v<T, TriangleIsect const *>) {
                return *v < *std::get<const TriangleIsect *>(o.vv);
            }
        },
                          vv);
    }
}

template<int D>
bool Crossing<D>::is_grid_vertex() const {
    return std::holds_alternative<VPtr>(vv);
}

template<int D>
bool Crossing<D>::is_edge_intersection() const {
    return std::holds_alternative<EPtr>(vv);
}

template<int D>
bool Crossing<D>::is_triangle_intersection() const {
    return std::holds_alternative<FPtr>(vv);
}

template<int D>
auto Crossing<D>::vertex_ptr() const -> const VType * {
    return std::visit([](auto &&val) -> const VType * {
        return static_cast<const VType *>(val);
    },
                      vv);
}

template<int D>
auto Crossing<D>::vertex() const -> const VType & {
    return *vertex_ptr();
}

template<int D>
auto Crossing<D>::mask() const -> typename VType::MaskType {
    return vertex().mask();
}

template<int D>
auto Crossing<D>::point() const -> Vec {
    return vertex().p();
}

template<int D>
Crossing<D>::operator std::string() const {
    std::stringstream ss;
    std::visit([&](auto &&v) {
        using T = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<T, VType const *>) {
            ss << "V";
        } else if constexpr (std::is_same_v<T, EdgeIsect const *>) {
            ss << "E"
               << "{" << v->edge_index << "}";
        } else if constexpr (std::is_same_v<T, TriangleIsect const *>) {
            ss << "F"
               << "{" << v->triangle_index << "}";
        } else {
            ss << "???";
        }
    },
               vv);
    if (index >= 0) {
        ss << "[" << index << "]" << std::string(vertex());
    } else {
        ss << "[?]" << std::string(vertex());
    }
    return ss.str();
}

template<int D>
EdgeCrossing<D>::operator std::string() const {
    std::stringstream ss;
    ss << "[" << index << "_" << index << ":" << edge_coord() << "]" << std::string(vertex());
    return ss.str();
}

template<int D>
int EdgeCrossing<D>::find_unbound(const VType &v) {
    //We are not allowed to have grid vertices
    assert(v.clamped_indices.count() == D - 1);
    for (int i = 0; i < D; ++i) {
        if (!v.clamped(i)) {
            return i;
        }
    }
    return -1;
}

template<int D>
EdgeCrossing<D>::EdgeCrossing(const VType *v, int i) : val(v), axis(find_unbound(*v)), index(i) {}

template<int D>
EdgeCrossing<D>::EdgeCrossing(const VType &v, int i) : EdgeCrossing(&v, i) {}

template<int D>
EdgeCrossing<D>::EdgeCrossing(const Crossing<D> &v, int i) : EdgeCrossing(&v.vertex(), i) {}

template<int D>
EdgeCrossing<D>::EdgeCrossing(const Crossing<D> &v) : EdgeCrossing(&v.vertex(), v.index) {}

template<int D>
bool EdgeCrossing<D>::operator<(const EdgeCrossing &o) const {
    ////If on the same line we can truely compare them, by their position along the line
    //assert(index == o.index);
    return edge_coord() < o.edge_coord();
}

template<int D>
double EdgeCrossing<D>::edge_coord() const {
    return vertex().quot(axis);
}

template<int D>
auto EdgeCrossing<D>::vertex() const -> const VType & {
    return *val;
}

template<int D>
auto EdgeCrossing<D>::mask() const -> typename VType::MaskType {
    return vertex().mask();
}

template<int D>
auto EdgeCrossing<D>::point() const -> Vec {
    return vertex().p();
}

}// namespace mandoline::construction
