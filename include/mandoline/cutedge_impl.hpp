#pragma once
#include "mandoline/cutedge.hpp"

namespace mandoline {

template<int D>
bool CutEdgeBase<D>::is_mesh_edge() const {
    return std::holds_alternative<int>(id);
}
template<int D>
bool CutEdgeBase<D>::is_axial_edge() const {
    return std::holds_alternative<std::array<int, 2>>(id);
}
template<int D>
int CutEdgeBase<D>::as_edge_id() const {
    return std::get<int>(id);
}
template<int D>
const std::array<int, 2> &CutEdgeBase<D>::as_axial_id() const {
    return std::get<std::array<int, 2>>(id);
}
template<int D>
int CutEdgeBase<D>::as_axial_axis() const {
    return std::get<std::array<int, 2>>(id)[0];
}
template<int D>
int CutEdgeBase<D>::as_axial_coord() const {
    return std::get<std::array<int, 2>>(id)[1];
}
template<int D>
CutEdgeBase<D>::operator std::string() const {
    std::stringstream ss;
    ss << "{";

    std::visit([&](auto &&f) {
        using T = std::decay_t<decltype(f)>;
        if constexpr (std::is_same_v<T, int>) {
            ss << "(" << f << ")";
        } else {
            auto [a, b] = f;
            ss << "(" << a << ":" << b << ")";
        }
    },
               id);
    ss << "[";
    std::copy(indices.begin(), indices.end(), std::ostream_iterator<int>(ss, ","));
    ss << "],";

    /*
    if (external_boundary) {
        auto [b, s] = *external_boundary;
        if (s) {
            std::cout << "{+" << b << "}";
        } else {
            std::cout << "{-" << b << "}";
        }
    }
    */
    ss << "}";

    return ss.str();
}
template<int D>
template<typename Derived>
void CutEdgeBase<D>::update_mask(const std::vector<Vertex<D>> &V, const mtao::geometry::grid::indexing::IndexerBase<D, Derived> &indexer) {
    int offset = indexer.size();
    auto get_mask = [&](int idx) -> coord_mask<D> {
        if (idx >= offset) {
            return V[idx - offset].mask();
        } else {
            return indexer.unindex(idx);
        }
    };
    coord_mask<D> &me = *this;
    me = get_mask(indices.begin()->front());
    for (auto &&ind : indices) {
        for (auto &&i : ind) {
            auto mask = get_mask(i);
            me &= mask;
        }
    }
}


}// namespace mandoline
