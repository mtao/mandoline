#include <mtao/geometry/grid/staggered_grid.hpp>
#include <balsa/eigen/types.hpp>
#include "mandoline/construction/vertex_types.hpp"

namespace mandoline::construction {
struct SubGridTransformer {
  public:
    using SGrid2 = mtao::geometry::grid::StaggeredGrid<double, 2>;
    using SGrid3 = mtao::geometry::grid::StaggeredGrid<double, 3>;
    using Vec2 = balsa::eigen::Vector<double, 2>;
    using Vec3 = balsa::eigen::Vector<double, 3>;
    using Coord2 = std::array<int, 2>;
    using Coord3 = std::array<int, 3>;
    SubGridTransformer(const SGrid2 &g2, const SGrid3 &g3, int axis, int c);

    Coord2 downgrade(const Coord3 &c3) const;
    Vec2 downgrade(const Vec3 &c3) const;
    std::bitset<2> downgrade(const std::bitset<3> &c3) const;

    Coord3 upgrade(const Coord2 &c2) const;
    Vec3 upgrade(const Vec2 &c2) const;
    std::bitset<3> upgrade(const std::bitset<2> &c2) const;

    Coord3 upgrade(const Coord2 &c2, int value) const;
    Vec3 upgrade(const Vec2 &c2, double value) const;
    std::bitset<3> upgrade(const std::bitset<2> &c3, bool value) const;

    int downgrade(int idx3) const;
    int upgrade(int idx2) const;
    int upgrade(int idx2, int value) const;
    EdgeIntersection<2> downgrade(const EdgeIntersection<3> &p) const;
    EdgeIntersection<3> upgrade(const EdgeIntersection<2> &p) const;

    TriangleIntersection<2> downgrade(const TriangleIntersection<3> &p) const;
    TriangleIntersection<3> upgrade(const TriangleIntersection<2> &p) const;

    Vertex<2> downgrade(const Vertex<3> &p) const;
    Vertex<3> upgrade(const Vertex<2> &p) const;

    void set_coord(int v) { c = v; }
    void set_axis(int v) { axis = v; }

  private:
    const SGrid2 &g2;
    const SGrid3 &g3;
    int axis;
    int c;
    template<typename Two, typename Three>
    Two _downgrade(const Three &c3) const;
    template<typename Two, typename Three, typename Scalar>
    Three _upgrade(const Two &in, Scalar v) const;
};
template<typename Two, typename Three>
Two SubGridTransformer::_downgrade(const Three &c3) const {
    Two c2;
    for (int i = 0; i < 2; ++i) {
        int d = (i + axis + 1) % 3;
        c2[i] = c3[d];
    }
    return c2;
}
template<typename Two, typename Three, typename Scalar>
Three SubGridTransformer::_upgrade(const Two &in, Scalar v) const {
    Three c3;
    c3[axis] = v;
    for (int i = 0; i < 2; ++i) {
        int d = (i + axis + 1) % 3;
        c3[d] = in[i];
    }
    return c3;
}
}// namespace mandoline::construction
