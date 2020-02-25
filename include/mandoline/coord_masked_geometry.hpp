#pragma once
#include <type_traits>
#include <iterator>
#include "mandoline/coord_mask.hpp"
#include <array>
#include <variant>
#include <vector>
#include <set>
#include "cutmesh.pb.h"
#include <mtao/geometry/volume.h>
#include <mtao/geometry/centroid.hpp>
#include "mandoline/proto_util.hpp"


namespace mandoline {
template<int D, typename IndexContainerType_>
struct CoordMaskedGeometry : public coord_mask<D> {
    using IndexContainerType = IndexContainerType_;
    using coord_type = typename coord_mask<D>::coord_type;
    CoordMaskedGeometry() = default;
    CoordMaskedGeometry(const CoordMaskedGeometry &) = default;
    CoordMaskedGeometry(CoordMaskedGeometry &&) = default;
    CoordMaskedGeometry &operator=(const CoordMaskedGeometry &) = default;
    CoordMaskedGeometry &operator=(CoordMaskedGeometry &&) = default;
    CoordMaskedGeometry(const coord_mask<D> &mask);
    CoordMaskedGeometry(const coord_mask<D> &mask, const IndexContainerType &indices);
    CoordMaskedGeometry(const coord_mask<D> &mask, IndexContainerType &&indices);
    template<typename Func>
    CoordMaskedGeometry(const IndexContainerType &indices, Func &&f);
    template<typename Func>
    CoordMaskedGeometry(IndexContainerType &&indices, Func &&f);

    const coord_mask<D> &mask() const;

    //Func is an object that retuns the mask of an index
    template<typename Func>
    static coord_mask<D> get_container_mask(const IndexContainerType &C, Func &&f);
    static coord_mask<D> get_container_mask(const IndexContainerType &C, const mtao::vector<Vertex<D>> &vertices);
    template<typename Func>
    static coord_type get_min_coord(const IndexContainerType &C, Func &&f);


    //coord_mask<D> get_container_mask(const IndexContainerType& C, const mtao::vector<Vertex<D>>& vertices);
    static coord_type get_min_coord(const IndexContainerType &C, const mtao::vector<Vertex<D>> &vertices);
    coord_type get_min_coord(const mtao::vector<Vertex<D>> &vertices) const;

    std::set<coord_type> possible_cells(const mtao::vector<Vertex<D>> &vertices) const;

    bool operator<(const CoordMaskedGeometry &other) const;


    IndexContainerType indices;
};

template<int D>
using CoordMaskedEdge = CoordMaskedGeometry<D, std::array<int, 2>>;
template<int D>
using CoordMaskedTriangle = CoordMaskedGeometry<D, std::array<int, 3>>;
template<int D>
using CoordMaskedSimplePolygon = CoordMaskedGeometry<D, std::vector<int>>;
template<int D>
using CoordMaskedPolygon = CoordMaskedGeometry<D, std::set<std::vector<int>>>;
template<int D>
using CoordMaskedPointSet = CoordMaskedGeometry<D, std::set<int>>;
}// namespace mandoline
#include "mandoline/coord_masked_geometry_impl.hpp"
