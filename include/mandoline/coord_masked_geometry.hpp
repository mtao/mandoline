#pragma once
#include <type_traits>
#include <iterator>
#include  "mandoline/coord_mask.hpp"
#include <array>
#include <variant>
#include <vector>
#include <set>
#include "cutmesh.pb.h"
#include <mtao/geometry/volume.h>
#include <mtao/geometry/centroid.hpp>
#include "mandoline/proto_util.hpp"


namespace mandoline {
    template <int D, typename IndexContainerType_>
        struct CoordMaskedGeometry: public coord_mask<D> {
            using IndexContainerType = IndexContainerType_;
            CoordMaskedGeometry() = default;
            CoordMaskedGeometry(const CoordMaskedGeometry&) = default;
            CoordMaskedGeometry(CoordMaskedGeometry&&) = default;
            CoordMaskedGeometry& operator=(const CoordMaskedGeometry&) = default;
            CoordMaskedGeometry& operator=(CoordMaskedGeometry&&) = default;
            CoordMaskedGeometry(const coord_mask<D>& mask): coord_mask<D>(mask) {}
            CoordMaskedGeometry(const coord_mask<D>& mask, const IndexContainerType& indices): coord_mask<D>(mask), indices(indices) {}
            CoordMaskedGeometry(const coord_mask<D>& mask, IndexContainerType&& indices): coord_mask<D>(mask), indices(indices) {}
            template <typename Func>
                CoordMaskedGeometry(const IndexContainerType& indices, Func&& f): coord_mask<D>(get_container_mask(indices,f)), indices(indices) {}
            template <typename Func>
                CoordMaskedGeometry(IndexContainerType&& indices, Func&& f): coord_mask<D>(get_container_mask(indices,f)), indices(indices) {}

            const coord_mask<D>& mask() const { return *this; }

            //Func is an object that retuns the mask of an index
            template <typename Func>
                static coord_mask<D> get_container_mask(const IndexContainerType& C, Func&& f) {
                    using S = std::decay_t<typename IndexContainerType::value_type>;
                    coord_mask<D> mask;
                    auto it = C.begin();
                    if constexpr(std::is_integral_v<S>) {
                        mask = f(*it);
                    } else {
                        mask = get_container_mask(*it,f);
                    }
                    ++it;
                    for(;it != C.end(); ++it) {
                        if constexpr(std::is_integral_v<S>) {
                            mask &= f(*it);
                        } else {
                            mask &= CoordMaskedGeometry<D,S>::get_container_mask(*it,f);
                        }
                    }
                    return mask;
                }

            bool operator<(const CoordMaskedGeometry& other) const {
                return indices < other.indices;
            }
            IndexContainerType indices;
        };

    template <int D>
        using CoordMaskedEdge = CoordMaskedGeometry<D,std::array<int,2>>;
    template <int D>
        using CoordMaskedTriangle = CoordMaskedGeometry<D,std::array<int,3>>;
    template <int D>
        using CoordMaskedSimplePolygon = CoordMaskedGeometry<D,std::vector<int>>;
    template <int D>
        using CoordMaskedPolygon = CoordMaskedGeometry<D,std::set<std::vector<int>>>;
    template <int D>
        using CoordMaskedPointSet = CoordMaskedGeometry<D,std::set<int>>;
}
