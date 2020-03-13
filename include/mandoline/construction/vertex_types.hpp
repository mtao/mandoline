#pragma once
#include <set>
#include <mtao/types.hpp>
#include <variant>
#include "mandoline/coord_mask.hpp"
#include "mandoline/vertex.hpp"
#include <map>


namespace mandoline::construction {
template<int D>
class Crossing;

template<int D, typename CrossingContainer>
std::map<const Vertex<D> *, int> crossing_indexer(const CrossingContainer &C);
template<int D, typename CrossingContainer>
std::map<int, const Vertex<D> *> inverse_crossing_indexer(const CrossingContainer &C);
template<int D, typename CrossingContainer>
void populate_crossing_indices(CrossingContainer &C);

template<int D, typename CrossingContainer, typename StaggeredGridType>
std::map<const Vertex<D> *, int> crossing_indexer(const CrossingContainer &C, const StaggeredGridType &g);


//Vertex that lies on the interior of an edge
template<int D>
struct EdgeIntersection : public Vertex<D> {
    //Definitions
    using VType = Vertex<D>;
    using Vec = typename VType::Vec;
    using coord_type = typename VType::coord_type;
    using VType::p;
    using VType::quot;
    using VType::coord;
    using VType::is_grid_vertex;

    //Members
    double edge_coord = 0.5;
    int edge_index = -1;

    //Constructors
    EdgeIntersection() = default;
    EdgeIntersection(const VType &v, double coord, int index = -1);
    EdgeIntersection(EdgeIntersection &&) = default;
    EdgeIntersection(const EdgeIntersection &) = default;
    EdgeIntersection &operator=(EdgeIntersection &&) = default;
    EdgeIntersection &operator=(const EdgeIntersection &) = default;

    //Member functions
    bool operator<(const EdgeIntersection &o) const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};


//Vertex that lies on the interior of an triangle
template<int D>
struct TriangleIntersection : public Vertex<D> {
    //Definitions
    using VType = Vertex<D>;
    using Vec = typename VType::Vec;
    using coord_type = typename VType::coord_type;
    using VType::p;
    using VType::quot;
    using VType::coord;
    using VType::is_grid_vertex;
    using VType::Vertex;

    //Members
    mtao::Vec3d bary_coord = Vec::Zero();
    int triangle_index = -1;

    //Member functions
    TriangleIntersection(const VType &is, const mtao::Vec3d &v, int triangle_index);

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template<int D>
struct Crossing {
    //Definitions
    using VType = Vertex<D>;
    using Vec = typename VType::Vec;
    using coord_type = typename VType::coord_type;
    using EdgeIsect = EdgeIntersection<D>;
    using TriangleIsect = TriangleIntersection<D>;
    using VPtr = VType const *;
    using EPtr = EdgeIsect const *;
    using FPtr = TriangleIsect const *;
    using VertexVariant = std::variant<VPtr, EPtr, FPtr>;

    //Members
    VertexVariant vv;
    int index;

    //Constructors
    template<typename T>
    Crossing(const T *v, int i = -1) : vv(v), index(i) {}
    template<typename T>
    Crossing(const T &v, int i = -1) : vv(&v), index(i) {}
    Crossing(const VertexVariant &v, int i = -1) : vv(v), index(i) {}
    Crossing() = default;
    Crossing(const Crossing &) = default;
    Crossing(Crossing &&) = default;
    Crossing &operator=(const Crossing &) = default;
    Crossing &operator=(Crossing &&) = default;

    //Member functions
    bool operator<(const Crossing &o) const;
    bool is_grid_vertex() const;
    bool is_edge_intersection() const;
    bool is_triangle_intersection() const;
    VPtr get_grid_vertex() const { return std::get<VPtr>(vv); }
    EPtr get_edge_intersection() const { return std::get<EPtr>(vv); }
    FPtr get_triangle_intersection() const { return std::get<FPtr>(vv); }

    //Conversions
    typename VType::MaskType mask() const;
    Vec point() const;
    operator std::string() const;
    const VType *vertex_ptr() const;
    const VType &vertex() const;
};


// treat a vertex as a point on an edge
template<int D>
struct EdgeCrossing {
    //Definitions
    using VType = Vertex<D>;
    using Vec = typename VType::Vec;
    using coord_type = typename VType::coord_type;

    //Members
    const VType *val = nullptr;
    int axis;
    int index;

    //Constructors
    EdgeCrossing(const VType *v, int i);
    EdgeCrossing(const VType &v, int i);
    EdgeCrossing(const Crossing<D> &v, int i);
    EdgeCrossing(const Crossing<D> &v);
    EdgeCrossing(const EdgeCrossing &) = default;
    EdgeCrossing(EdgeCrossing &&) = default;
    EdgeCrossing &operator=(const EdgeCrossing &) = default;
    EdgeCrossing &operator=(EdgeCrossing &&) = default;
    //Construction helper, checks for being on the interior of an edge
    static int find_unbound(const VType &v);

    //Member functions
    bool operator<(const EdgeCrossing &o) const;
    double edge_coord() const;

    //Conversions
    const VType &vertex() const;
    typename VType::MaskType mask() const;
    Vec point() const;
    operator std::string() const;
};
}// namespace mandoline::construction

#include "vertex_types_impl.hpp"
