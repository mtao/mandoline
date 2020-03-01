#pragma once
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/iterable.hpp>
#include <mtao/iterator/zip.hpp>
#include <mtao/iterator/interval.hpp>
#include <mtao/geometry/barycentric.hpp>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/types.hpp>
#include <variant>
#include <optional>
#include <map>
#include <bitset>
#include "mandoline/construction/vertex_types.hpp"
#include "mandoline/cutedge.hpp"
#include "mandoline/cutface.hpp"
#include <mtao/geometry/mesh/halfedge.hpp>
#include "mandoline/construction/face_collapser.hpp"


namespace mandoline::construction {

//Upgrades intersections one level at a time
template<int D>
std::array<std::map<coord_mask<D>, std::set<const Vertex<D> *>>, D + 1> convex_grid_condensation(std::set<const Vertex<D> *> &V);


template<int D, typename Derived_>
struct IntersectionsBase : public coord_mask<D> {
    using Derived = Derived_;
    using VType = Vertex<D>;
    using Edge = std::array<int, 2>;
    using VPtrEdge = std::array<const VType *, 2>;
    auto &derived() { return *static_cast<Derived *>(this); }
    const auto &derived() const { return *static_cast<const Derived *>(this); }
    std::set<VPtrEdge> vptr_edges() const {
        return derived().vptr_edges();
    }

    std::set<Edge> edges(const std::map<const VType *, int> &indexer) const {
        auto vedges = vptr_edges();
        std::set<Edge> ret;
        for (auto &&[a, b] : vedges) {
            Edge e{ { indexer.at(a), indexer.at(b) } };
            if (e[0] == e[1]) continue;
            std::sort(e.begin(), e.end());
            ret.emplace(e);
        }

        return ret;
    }

    template<typename Vertices>
    static coord_mask<D> get_container_mask(const Vertices &C) {
        coord_mask<D> mask;
        auto it = C.begin();
        mask = (*it)->mask();
        ++it;
        for (; it != C.end(); ++it) {
            mask &= (*it)->mask();
        }
        return mask;
    }
    template<typename Vertices>
    void set_container_mask(const Vertices &C) {
        this->coord_mask<D>::operator=(get_container_mask(C));
    }


    const coord_mask<D> &mask() const { return *this; }
};


template<int D>
struct EdgeIntersections : public IntersectionsBase<D, EdgeIntersections<D>> {
    using Edges = mtao::ColVectors<int, 2>;
    using EdgeIsect = EdgeIntersection<D>;
    using VType = Vertex<D>;
    using VPtrEdge = std::array<const VType *, 2>;
    using Vec = typename VType::Vec;
    using CoordType = typename VType::CoordType;
    using SGType = mtao::geometry::grid::StaggeredGrid<double, D>;
    using Base = IntersectionsBase<D, EdgeIntersections<D>>;
    using Base::mask;
    int edge_index = -1;
    VPtrEdge vptr_edge;
    std::vector<EdgeIsect> intersections;
    std::array<double, D> tangent;
    EdgeIntersections() = default;
    EdgeIntersections(const EdgeIntersections &) = default;
    EdgeIntersections(EdgeIntersections &&) = default;
    EdgeIntersections &operator=(const EdgeIntersections &) = default;
    EdgeIntersections &operator=(EdgeIntersections &&) = default;
    EdgeIntersections(const mtao::vector<VType> &V, const Edges &E, int index = -1) : edge_index(index) {
        auto e = E.col(edge_index);
        auto &a = V[e(0)];
        auto &b = V[e(1)];
        vptr_edge = VPtrEdge{ { &a, &b } };
        Base::set_container_mask(vptr_edge);
        assert(mask().count() <= 2);
    }
    EdgeIntersections(const VPtrEdge &e, int index = -1) : vptr_edge(e), edge_index(index) {
        Base::set_container_mask(vptr_edge);
        assert(mask().count() <= 2);
    }
    EdgeIntersections(const VType &a, const VType &b, int index = -1) : EdgeIntersections(VPtrEdge{ { &a, &b } }, index) {}
    void bake(const std::optional<SGType> &grid = {});

    bool is_cut() const { return intersections.empty(); }

    size_t edge_size() const {
        return intersections.size() + 1;
    }
    size_t vertex_size() const {
        return intersections.size() + 2;
    }
    void clear() {
        intersections.clear();
        Base::set_container_mask(vptr_edge);
    }
    std::vector<Crossing<D>> gvertices() const {
        std::vector<Crossing<D>> ret;
        ret.reserve(vertex_size());

        auto &&[a, b] = vptr_edge;
        ret.emplace_back(a);


        std::map<double, const EdgeIsect *> im;
        for (auto &&i : intersections) {
            if (i.edge_coord > 0 && i.edge_coord < 1) {
                im[i.edge_coord] = &i;
            }
        }
        for (auto &&[t, ptr] : im) {
            ret.emplace_back(ptr);
        }
        ret.emplace_back(b);


        return ret;
    }
    EdgeIsect from_coord(double t) const {
        auto &a = *vptr_edge[0];
        auto &b = *vptr_edge[1];
        EdgeIsect is{ a.lerp(b, t), t, edge_index };

        coord_mask<D>::clamp(is);
        return is;
    }
    double get_coord(const VType& v) const {
        Vec a = vptr_edge[0]->p();
        Vec b = vptr_edge[1]->p();
        Vec ba = (b-a).eval();
        return (v.p() - a).dot(ba) / ba.squaredNorm();
    }

    std::set<VPtrEdge> vptr_edges() const {
        std::set<VPtrEdge> ret;
        auto gv = gvertices();

        auto ival = mtao::iterator::interval<2>(gv);
        for (auto &&[a, b] : ival) {
            ret.emplace(VPtrEdge{ { a.vertex_ptr(), b.vertex_ptr() } });
        }
        return ret;
    }
};
template<int D>
struct TriangleIntersections : public IntersectionsBase<D, TriangleIntersections<D>> {
    using EdgeIsect = EdgeIntersection<D>;
    using EdgeIsects = EdgeIntersections<D>;
    using TriIsect = TriangleIntersection<D>;
    using SGType = typename EdgeIntersections<D>::SGType;
    using VType = Vertex<D>;
    using VPtrEdge = std::array<const VType *, 2>;
    using VPtrTri = std::array<const VType *, 3>;
    //using VPtrTri = std::vector<const VType*>;
    using Edges = mtao::ColVectors<int, 2>;
    using Faces = mtao::ColVectors<int, 3>;
    using Edge = std::array<int, 2>;
    using Face = std::array<int, 3>;
    using Vec = typename VType::Vec;
    using CoordType = typename VType::CoordType;
    using Base = IntersectionsBase<D, TriangleIntersections<D>>;
    using Base::mask;
    std::map<VPtrEdge, EdgeIntersections<D>> edge_intersections;
    VPtrTri vptr_tri;
    int triangle_index;
    std::map<const EdgeIsect *, const TriIsect *> edge_to_triangle_map;
    std::array<const EdgeIntersections<D> *, 3> edge_isects;
    std::vector<TriIsect> intersections;
    std::array<Edge, 3> bary_indices;

    TriangleIntersections(const mtao::vector<VType> &V, const Faces &F, const Edges &E, const Faces &FEA, const std::vector<EdgeIntersections<D>> &eis, int index) : triangle_index(index) {
        auto fea = FEA.col(index);
        auto f = F.col(index);

        for (auto &&[i, gv, ei] : mtao::iterator::enumerate(vptr_tri, edge_isects)) {
            gv = &V[f(i)];
            ei = &eis[fea(i)];
            //if this is the gridvertex in an edge, make hte bary_indices structure remember it
        }
        for (auto &&[i, gv] : mtao::iterator::enumerate(vptr_tri)) {
            for (auto &&[bi, ei] : mtao::iterator::zip(bary_indices, edge_isects)) {
                for (auto &&[bie, egv] : mtao::iterator::zip(bi, ei->vptr_edge)) {
                    if (egv == gv) {
                        bie = i;
                    }
                }
            }
        }
        Base::set_container_mask(vptr_tri);
        assert(mask().count() <= 1);
    }

    bool is_cut() const { return edge_intersections.empty(); }

    void clear() {
        intersections.clear();
        edge_intersections.clear();
        edge_to_triangle_map.clear();
        Base::set_container_mask(vptr_tri);
    }


    size_t edge_size() const {

        size_t size = 0;
        for (auto &&[gv, e] : edge_intersections) {
            size += e.edge_size();
        }
        return size;
    }
    size_t vertex_size() const {
        size_t size = 3 + intersections.size();
        for (auto &&[gv, e] : edge_intersections) {
            size += e.intersections.size();
        }
        return size;
    }
    std::map<const VType *, std::array<double, 3>> barycentric_coords() const {
        std::map<const VType *, std::array<double, 3>> barys;
        std::array<double, 3> bary;
        auto bary_map = mtao::eigen::stl2eigen(bary);

        //vertex indices

        // write the barycentric stuff
        for (auto &&[i, v] : mtao::iterator::enumerate(vptr_tri)) {
            bary_map = mtao::Vec3d::Unit(i);
            barys[v] = bary;
        }

        //in-face indices/bary
        for (auto &&v : intersections) {
            bary_map = v.bary_coord;
            barys[&v] = bary;
        }

        //per-edge in-edge indices/bary
        for (auto &&[eptr, ei] : mtao::iterator::zip(edge_isects, bary_indices)) {
            for (auto &&v : eptr->intersections) {
                bary_map = edge_bary(ei, v.edge_coord);
                barys[&v] = bary;
            }
        }

        return barys;
    }

    std::vector<Crossing<D>> boundary_gvertices() const {
        std::vector<Crossing<D>> ret;

        for (auto &&v : vptr_tri) {
            ret.emplace_back(v);
        }
        for (auto &&eisptr : edge_isects) {
            auto &&e = eisptr->intersections;
            ret.insert(ret.end(), e.begin(), e.end());
        }

        return ret;
    }

    //Returns the boundary loop in the reverse order
    std::vector<Crossing<D>> boundary_vptr_loop() const;
    //a single boundary edge on the outside to distinguish
    VPtrEdge boundary_vptr_edge() const;
    std::vector<Crossing<D>> gvertices() const {
        std::vector<Crossing<D>> ret = boundary_gvertices();

        for (auto &&p : intersections) {
            ret.emplace_back(p);
        }

        return ret;
    }
    std::vector<int> boundary_loop(const std::map<const VType *, int> &indexer) const {
        auto bvl = boundary_vptr_loop();
        std::vector<int> loop(bvl.size());
        std::transform(bvl.begin(), bvl.end(), loop.begin(), [&](auto &&c) {
            return indexer.at(c.vertex_ptr());
        });
        return loop;
    }
    Edge boundary_edge(const std::map<const VType *, int> &indexer) const {
        auto bve = boundary_vptr_edge();
        Edge e;
        std::transform(bve.begin(), bve.end(), e.begin(), [&](auto &&vptr) {
            return indexer.at(vptr);
        });
        return e;
    }

    //Only edge vptrs are bad and need to be transformed
    const VType *get_vptr(const Crossing<D> &c) const {
        if (c.is_edge_intersection()) {
            if (auto it = edge_to_triangle_map.find(c.get_edge_intersection()); it != edge_to_triangle_map.end()) {
                return it->second;
            } else {
                return c.vertex_ptr();
            }
        } else {
            return c.vertex_ptr();
        }
    }

    std::set<Edge> edges(const std::map<const VType *, int> &indexer) const {
        auto vedges = vptr_edges();
        std::set<Edge> ret;
        for (auto &&[a, b] : vedges) {
            Edge e{ { indexer.at(a), indexer.at(b) } };
            if (e[0] == e[1]) continue;
            std::sort(e.begin(), e.end());
            ret.emplace(e);
        }

        return ret;
    }
    std::set<Edge> nobdry_edges(const std::map<const VType *, int> &indexer) const {
        auto vedges = vptr_edges(true);
        std::set<Edge> ret;
        for (auto &&[a, b] : vedges) {
            Edge e{ { indexer.at(a), indexer.at(b) } };
            if (e[0] == e[1]) continue;
            std::sort(e.begin(), e.end());
            ret.emplace(e);
        }

        return ret;
    }
    std::set<VPtrEdge> vptr_edges(bool get_boundary = false) const {
        std::set<VPtrEdge> ret;
        for (auto &&[gv, e] : edge_intersections) {
            auto isects = e.gvertices();
            VPtrEdge vpe;
            for (int i = 0; i < isects.size() - 1; ++i) {
                vpe[0] = get_vptr(isects[i]);
                vpe[1] = get_vptr(isects[i + 1]);

                ret.emplace(vpe);
            }
        }
        if (get_boundary) {
            for (auto &&eisptr : edge_isects) {
                auto &&es = eisptr->vptr_edges();
                ret.insert(es.begin(), es.end());
            }
        }
        return ret;
    }

    mtao::Vec3d N() const {
        if constexpr (D == 3) {
            mtao::Vec3d a = vptr_tri[1]->p() - vptr_tri[0]->p();
            mtao::Vec3d b = vptr_tri[2]->p() - vptr_tri[0]->p();
            return a.cross(b).normalized();
        } else {
            mtao::Vec2d a = vptr_tri[1]->p() - vptr_tri[0]->p();
            mtao::Vec2d b = vptr_tri[2]->p() - vptr_tri[0]->p();
            return mtao::Vec3d::Unit(2);
            //return (a.x() * b.y() - a.y() * b.x()) * mtao::Vec3d::Unit(2);
        }
    }
    mtao::Vec3d get_bary(const VType &v) const {
        //try to build barycentric
        Eigen::Matrix<double, 3, 3> V;
        for (auto &&[i, ptr] : mtao::iterator::enumerate(vptr_tri)) {
            V.col(i) = ptr->p();
        }
        return mtao::geometry::barycentric_simplicial(V, v.p());
    }

    //line lerp coord is inverted from barycentric coordinates
    static mtao::Vec3d edge_bary(const std::array<int, 2> &indices, double t) {
        mtao::Vec3d v = mtao::Vec3d::Zero();
        v(indices[0]) = 1 - t;
        v(indices[1]) = t;
        return v;
    }
    mtao::Vec3d edge_bary(const VPtrEdge &indices, double t) const {
        Edge e;
        for (auto &&[i, gvptr] : mtao::iterator::enumerate(vptr_tri)) {
            for (auto &&[v, ind] : mtao::iterator::zip(e, indices)) {
                if (gvptr == ind) {
                    v = i;
                }
            }
        }
        return edge_bary(e, t);
    }
    mtao::Vec3d edge_bary(const EdgeIsects &ei, double t) const {
        Edge e;
        for (auto &&[eptr, bi] : mtao::iterator::zip(edge_isects, bary_indices)) {
            if (&ei == eptr) {
                return edge_bary(bi, t);
            }
        }
        auto p = get_bary(ei.from_coord(t));

        return p;
    }
    /*
               mtao::Vec3d edge_bary(const EdgeIsect& ei, double t) const {
               Edge e;
               for(auto&& [eptr,bi]: mtao::iterator::zip(edge_isects,bary_indices)) {
               if(&ei == eptr) {
               return edge_bary(bi,t);
               }
               }
               auto p = get_bary(ei.from_coord(t));

               }
               */
    void bake(const std::optional<SGType> &grid = {});

    mtao::Vec3d get_bary(const VPtrEdge &indices, double t) const {
        auto B = edge_barys(indices);
        return (1 - t) * B.col(0) + t * B.col(1);
    }

    std::set<std::vector<const VType *>> vptr_faces() const;
    std::set<std::vector<int>> faces(const std::map<const VType *, int> &indexer) const;

    TriIsect from_coord(const Vec &B) const {
        VType gv = B(0) * *vptr_tri[0] + B(1) * *vptr_tri[1] + B(2) * *vptr_tri[2];
        TriIsect is{ gv, triangle_index, B };
        coord_mask<D>::clamp(is);
    }
    TriIsect from_coord(const VPtrEdge &eis, const double t) const {
        return from_coord(get_bary(eis, t));
    }
    TriIsect from_coord(const EdgeIntersections<D> &eis, const double t) const {
        return from_coord(eis.vptr_edge, t);
    }
};
}// namespace mandoline::construction
#include "mandoline/construction/facet_intersections_impl.hpp"
