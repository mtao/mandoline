#pragma once
#include "mandoline/coord_masked_geometry.hpp"
#include "mandoline/vertex.hpp"
#include <mtao/geometry/grid/grid.h>
#include <igl/solid_angle.h>

namespace mandoline {

template<int D>
struct CutMeshEdge : public CoordMaskedEdge<D> {

    using Base = CoordMaskedEdge<D>;
    using Base::mask;
    using Base::indices;
    using Base::Base;


    CutMeshEdge(const coord_mask<D> &m, const std::array<int, 2> &indices, int pid) : Base(m, indices), parent_eid(pid) {}
    CutMeshEdge() = default;
    CutMeshEdge(const CutMeshEdge &) = default;
    CutMeshEdge(CutMeshEdge &&) = default;
    CutMeshEdge &operator=(const CutMeshEdge &) = default;
    CutMeshEdge &operator=(CutMeshEdge &&) = default;
    size_t size() const { return indices.size(); }

    int parent_eid = -1;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template<int D>
struct CutEdgeBase : public CoordMaskedEdge<D> {

    using Base = CoordMaskedEdge<D>;
    using Base::mask;
    using Base::indices;
    using Base::Base;
    using Vec = mtao::Vector<double, D>;

    using IDType = std::variant<int, std::array<int, 2>>;


    CutEdgeBase(const coord_mask<D> &m, const std::array<int, 2> &indices, const IDType &id) : Base(m, indices), id(id) {}
    CutEdgeBase(const coord_mask<D> &m, const std::array<int, 2> &indices, const std::array<int, 2> &id) : Base(m, indices), id(id) {}

    CutEdgeBase(const CutMeshEdge<D> &E) : CutEdgeBase(E.mask(), E.indices, E.parent_eid) {}
    CutEdgeBase() = default;
    CutEdgeBase(const CutEdgeBase &) = default;
    CutEdgeBase(CutEdgeBase &&) = default;
    CutEdgeBase &operator=(const CutEdgeBase &) = default;
    CutEdgeBase &operator=(CutEdgeBase &&) = default;

    //a tuple where the first entry is the cell-grid index of hte external cell, and the bool is whether it is below it or not
    //  i.e given two cells cm cp bounded by f,
    //  [cm] f [cp], if [cm] is not part of the staggered grid then we store
    //  cm,1
    //  otherwise we store
    //  cp,0
    // this sign is with respect to "the boundary operator of hte cell is of positive sign", which may be backwards
    std::optional<std::tuple<int, bool>> external_boundary = {};
    bool is_mesh_edge() const;
    bool is_axial_edge() const;
    // when this face comes from a mesh,  as_edge_id returns the face index in the mesh
    int as_edge_id() const;

    // when the face comes from a grid plane mesh/axial
    const std::array<int, 2> &as_axial_id() const;
    int as_axial_axis() const;
    int as_axial_coord() const;


    template<typename Derived>
    void update_mask(const std::vector<Vertex<D>> &V, const mtao::geometry::grid::indexing::IndexerBase<D, Derived> &indexer);

    operator std::string() const;


    IDType id;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template<int D>
struct CutEdge {};

template<>
struct CutEdge<2> : public CutEdgeBase<2> {
    using Base = CutEdgeBase<2>;
    using Base::mask;
    using Base::indices;
    using Base::Base;

    CutEdge(const coord_mask<2> &m, const std::array<int, 2> &indices, const IDType &id, const Vec &N) : Base(m, indices, id), N(N) {}
    CutEdge(const coord_mask<2> &m, const std::array<int, 2> &indices, const std::array<int, 2> &id) : Base(m, indices, id), N(Vec::Unit(m.bound_axis())) {}
    CutEdge(const coord_mask<2> &m, const std::array<int, 2> &indices) : CutEdge(m, indices, std::array<int, 2>{ { m.bound_axis(), *m[m.bound_axis()] } }) {}
    CutEdge(const CoordMaskedEdge<2> &E) : CutEdge(E.mask(), E.indices) {}

    CutEdge(const CutMeshEdge<2> &E, const Vec &N) : CutEdge(E.mask(), E.indices, E.parent_eid, N) {}

    Vec N;

    // TODO:
    //void  serialize(protobuf::CutEdge& edge) const ;
};
template<>
struct CutEdge<3> : public CutEdgeBase<3> {
    using Base = CutEdgeBase<3>;
    using Base::mask;
    using Base::indices;
    using Base::Base;

    //CutEdge(const coord_mask<3>& m, const std::array<int,2>& indices, const IDType& id): Base(m,{indices}), id(id) {}
    //CutEdge(const coord_mask<3>& m, const std::array<int,2>& indices, const std::array<int,2>& id): Base(m,indices), id(id)  {}

    //CutEdge(const CutMeshEdge<3>& E): CutEdge(E.mask(),E.indices,E.parent_eid) {}
    //CutEdge(const CoordMaskedEdge<3>& E): CutEdge(E) {}

    // TODO:
    //void  serialize(protobuf::CutEdge& edge) const ;
};


}// namespace mandoline


#include "mandoline/cutface_impl.hpp"
