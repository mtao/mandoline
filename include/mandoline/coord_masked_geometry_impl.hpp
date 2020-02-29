#pragma once
#include "mandoline/coord_masked_geometry.hpp"
namespace mandoline {

template<int D, typename IndexContainerType>
CoordMaskedGeometry<D, IndexContainerType>::CoordMaskedGeometry(const coord_mask<D> &mask) : coord_mask<D>(mask) {}
template<int D, typename IndexContainerType>
CoordMaskedGeometry<D, IndexContainerType>::CoordMaskedGeometry(const coord_mask<D> &mask, const IndexContainerType &indices) : coord_mask<D>(mask), indices(indices) {}
template<int D, typename IndexContainerType>
CoordMaskedGeometry<D, IndexContainerType>::CoordMaskedGeometry(const coord_mask<D> &mask, IndexContainerType &&indices) : coord_mask<D>(mask), indices(indices) {}
template<int D, typename IndexContainerType>
template<typename Func>
CoordMaskedGeometry<D, IndexContainerType>::CoordMaskedGeometry(const IndexContainerType &indices, Func &&f) : coord_mask<D>(get_container_mask(indices, f)), indices(indices) {}
template<int D, typename IndexContainerType>
template<typename Func>
CoordMaskedGeometry<D, IndexContainerType>::CoordMaskedGeometry(IndexContainerType &&indices, Func &&f) : coord_mask<D>(get_container_mask(indices, f)), indices(indices) {}

template<int D, typename IndexContainerType>
const coord_mask<D> &CoordMaskedGeometry<D, IndexContainerType>::mask() const { return *this; }

//Func is an object that retuns the mask of an index
template<int D, typename IndexContainerType>
template<typename Func>
coord_mask<D> CoordMaskedGeometry<D, IndexContainerType>::get_container_mask(const IndexContainerType &C, Func &&f) {
    using S = std::decay_t<typename IndexContainerType::value_type>;
    coord_mask<D> mask;
    auto it = C.begin();
    if constexpr (std::is_integral_v<S>) {
        mask = f(*it);
    } else {
        mask = get_container_mask(*it, f);
    }
    ++it;
    for (; it != C.end(); ++it) {
        if constexpr (std::is_integral_v<S>) {
            mask &= f(*it);
        } else {
            mask &= CoordMaskedGeometry<D, S>::get_container_mask(*it, f);
        }
    }
    return mask;
}

template<int D, typename IndexContainerType>
coord_mask<D> CoordMaskedGeometry<D, IndexContainerType>::get_container_mask(const IndexContainerType &C, const mtao::vector<Vertex<D>> &vertices) {
    return get_container_mask(C, [&](size_t idx) -> coord_mask<D> {
        return vertices.at(idx).mask();
    });
}
template<int D, typename IndexContainerType>
template<typename Func>
auto CoordMaskedGeometry<D, IndexContainerType>::get_min_coord(const IndexContainerType &C, Func &&f) -> coord_type {
    using S = std::decay_t<typename IndexContainerType::value_type>;
    coord_type coord;
    auto it = C.begin();
    if constexpr (std::is_integral_v<S>) {
        coord = f(*it);
    } else {
        coord = CoordMaskedGeometry<D,S>::get_min_coord(*it, f);
    }
    ++it;
    for (; it != C.end(); ++it) {
        if constexpr (std::is_integral_v<S>) {
            auto v = f(*it);
            for (int j = 0; j < D; ++j) {
                coord[j] = std::min(coord[j], v[j]);
            }
        } else {
            auto v = CoordMaskedGeometry<D, S>::get_min_coord(*it, f);
            for (int j = 0; j < D; ++j) {
                coord[j] = std::min(coord[j], v[j]);
            }
        }
    }
    return coord;
}

template<int D, typename IndexContainerType>
auto CoordMaskedGeometry<D, IndexContainerType>::get_min_coord(const mtao::vector<Vertex<D>> &verts) const -> coord_type {
    return get_min_coord(indices, verts);
    /*
            using S = std::decay_t<typename IndexContainerType::value_type>;
            coord_mask<D> mask;
            auto it = C.begin();
            if constexpr(std::is_integral_v<S>) {
                mask = verts[*it];
            } else {
                mask = get_container_mask(*it,verts);
            }
            ++it;
            for(;it != C.end(); ++it) {
                if constexpr(std::is_integral_v<S>) {
                    auto v = verts[*it].coord:
                    for(int j = 0; j < D; ++j) {
                        mask[j] = std::min(mask[j],v[j]);
                    }
                } else {
                    auto v = CoordMaskedGeometry<D,S>::get_container_mask(*it,verts);
                    for(int j = 0; j < D; ++j) {
                        mask[j] = std::min(mask[j],v[j]);
                    }
                }
            }
            return mask;
            */
}
template<int D, typename IndexContainerType>
auto CoordMaskedGeometry<D, IndexContainerType>::get_min_coord(const IndexContainerType &indices, const mtao::vector<Vertex<D>> &vertices) -> coord_type {

    return get_min_coord(indices, [&](size_t idx) -> coord_type {
        return vertices.at(idx).coord;
    });
}

template<int D, typename IndexContainerType>
auto CoordMaskedGeometry<D, IndexContainerType>::possible_cells(const mtao::vector<Vertex<D>> &vertices) const -> std::set<coord_type> {
    // pick out the bottom left corner to help the mask have a reference
    return mask().possible_cells(get_min_coord(indices, vertices));
}

template<int D, typename IndexContainerType>
bool CoordMaskedGeometry<D, IndexContainerType>::operator<(const CoordMaskedGeometry &other) const {
    return indices < other.indices;
}
}// namespace mandoline
