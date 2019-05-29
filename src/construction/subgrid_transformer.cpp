#include "mandoline/construction/subgrid_transformer.hpp"

namespace mandoline::construction {

    SubGridTransformer::SubGridTransformer(const SGrid2& g2, const SGrid3& g3, int axis, int c): g2(g2), g3(g3), axis(axis), c(c) {}

    auto SubGridTransformer::downgrade(const Coord3& c3) const -> Coord2  {return _downgrade<Coord2,Coord3>(c3);}
    auto SubGridTransformer::downgrade(const Vec3& c3) const -> Vec2 {return _downgrade<Vec2,Vec3>(c3);}

    auto SubGridTransformer::upgrade(const Coord2& c2) const -> Coord3 {return upgrade(c2,c);}
    auto SubGridTransformer::upgrade(const Vec2& c2) const -> Vec3 {return upgrade(c2,0);}

    auto SubGridTransformer::upgrade(const Coord2& c2, int value) const -> Coord3 {return _upgrade<Coord2,Coord3>(c2,value);}
    auto SubGridTransformer::upgrade(const Vec2& c2, double value) const -> Vec3 {return _upgrade<Vec2,Vec3>(c2,value);}


    auto SubGridTransformer::downgrade(const std::bitset<3>& c3) const -> std::bitset<2> {return _downgrade<std::bitset<2>,std::bitset<3>>(c3);}

    auto SubGridTransformer::upgrade(const std::bitset<2>& c2) const -> std::bitset<3> {return upgrade(c2,true);}

    auto SubGridTransformer::upgrade(const std::bitset<2>& c2, bool value) const -> std::bitset<3> {return _upgrade<std::bitset<2>,std::bitset<3>>(c2,value);}

    int  SubGridTransformer::downgrade(int idx3) const {
        return g2.index(_downgrade<Coord2,Coord3>(g3.unindex(idx3)));
    }
    int  SubGridTransformer::upgrade(int idx2) const {
        return upgrade(idx2,c);
    }
    int  SubGridTransformer::upgrade(int idx2, int value) const {
        return g3.index(_upgrade<Coord2,Coord3>(g2.unindex(idx2),value));
    }
    EdgeIntersection<2>  SubGridTransformer::downgrade(const EdgeIntersection<3>& p) const {
        auto gv = downgrade(static_cast<const Vertex<3>&>(p));
        return {gv,p.edge_coord,p.edge_index};
    }
    EdgeIntersection<3>  SubGridTransformer::upgrade(const EdgeIntersection<2>& p) const {
        auto gv = upgrade(static_cast<const Vertex<2>&>(p));
        return {gv,p.edge_coord,p.edge_index};
    }

    TriangleIntersection<2>  SubGridTransformer::downgrade(const TriangleIntersection<3>& p) const {
        auto gv = downgrade(static_cast<const Vertex<3>&>(p));
        return {gv,p.bary_coord,p.triangle_index};
    }
    TriangleIntersection<3>  SubGridTransformer::upgrade(const TriangleIntersection<2>& p) const {
        auto gv = upgrade(static_cast<const Vertex<2>&>(p));
        return {gv,p.bary_coord,p.triangle_index};
    }

    Vertex<2>  SubGridTransformer::downgrade(const Vertex<3>& p) const {
        return {downgrade(p.coord),downgrade(p.quot),downgrade(p.clamped_indices)};

    }
    Vertex<3>  SubGridTransformer::upgrade(const Vertex<2>& p) const {
        return {upgrade(p.coord),upgrade(p.quot),upgrade(p.clamped_indices)};
    }
}
