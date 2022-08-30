#pragma once
#include <balsa/eigen/stack.hpp>
#include <spdlog/spdlog.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <mtao/colvector_loop.hpp>
#include <balsa/eigen/stl2eigen.hpp>
#include <mtao/iterator/enumerate.hpp>

#include "mandoline/construction/cutdata.hpp"
#include "mandoline/construction/generator.hpp"
#include "mandoline/interpolated_edge.hpp"

namespace mandoline::construction {
template <int D>
std::vector<Vertex<D>> CutCellEdgeGenerator<D>::all_vertices() const {
    std::vector<Vertex<D>> ret(num_vertices());
    for (int i = 0; i < num_vertices(); ++i) {
        ret[i] = GV(i);
    }
    return ret;
}

template <int D>
bool CutCellEdgeGenerator<D>::valid() const {
    bool ret = true;
    ret &= valid_grid();
    return ret;
}
template <int D>
bool CutCellEdgeGenerator<D>::valid_grid() const {
    bool ret = true;
    auto N = vertex_shape();
    int min = *std::min_element(N.begin(), N.end());
    if (min <= 2) {
        if constexpr (D == 2) {
            spdlog::warn("Grid shape is too small ({} {}) (min must be >= 2)",
                         N[0], N[1]);
        } else if constexpr (D == 3) {
            spdlog::warn(
                "Grid shape is too small ({} {} {}) (min must be >= 2)", N[0],
                N[1], N[2]);
        } else {
            spdlog::warn("Grid shape is too small");
        }

        ret = false;
    }
    auto bbox = this->bbox();
    if ((bbox.sizes().minCoeff() <= 0)) {
        if constexpr (D == 2) {
            spdlog::warn(
                "Empty or negative shaped grid (bbox of [{} {}] => [{} {}]; "
                "grid edges of [{} {}]",
                bbox.min().x(), bbox.min().y(), bbox.max().x(), bbox.max().y(),
                bbox.sizes().x(), bbox.sizes().y());
        } else if constexpr (D == 3) {
            spdlog::warn(
                "Empty or negative shaped grid (bbox of [{} {} {}] => [{} {} "
                "{}]; grid edges of [{} {} {}]",
                bbox.min().x(), bbox.min().y(), bbox.min().z(), bbox.max().x(),
                bbox.max().y(), bbox.max().z(), bbox.sizes().x(),
                bbox.sizes().y(), bbox.sizes().z());
        } else {
            spdlog::warn("Empty or negative shaped grid {}",
                         bbox.sizes().minCoeff());
        }
    }
    return true;
}
template <int D>
auto CutCellEdgeGenerator<D>::colvecs_to_vecvector(const ColVecs &V)
    -> VecVector {
    mtao::vector<balsa::eigen::Vector<double, D>> stlp(V.cols());
    for (auto &&[i, v] : mtao::iterator::enumerate(stlp)) {
        v = V.col(i);
    }
    return stlp;
}

template <int D>
CutCellEdgeGenerator<D>::CutCellEdgeGenerator(const ColVecs &V,
                                              const StaggeredGrid &g,
                                              std::optional<double> threshold)
    : CutCellEdgeGenerator(colvecs_to_vecvector(V), g, threshold) {}
template <int D>
CutCellEdgeGenerator<D>::CutCellEdgeGenerator(const VecVector &V,
                                              const StaggeredGrid &g,
                                              std::optional<double> threshold)
    : StaggeredGrid(g),
      m_data(g.vertex_grid()),
      m_origV(V),
      m_active_grid_cell_mask(g.cell_shape()) {
    assert(g.vertex_grid().valid());

    // CutCellEdgeGenerator<D>::CutCellEdgeGenerator(const VecVector& V, const
    // StaggeredGrid& g, std::optional<double> threshold): StaggeredGrid(g),
    // m_origV(V) {
    m_data.m_V.resize(V.size());

    update_vertices(V);
}

template <int D>
void CutCellEdgeGenerator<D>::update_vertices(
    const VecVector &V, const std::optional<double> &threshold) {
    m_origV = V;
    assert(V.size() == m_data.m_V.size());
    double mythresh = -1;
    if (threshold) {
        mythresh = *threshold;
    }
    if (mythresh < 0) {
        int v = *std::max_element(vertex_shape().begin(), vertex_shape().end());
        v = 2 * std::max<int>(1, v);
        // mythresh =  v * threshold_epsilon;
        mythresh = v * 1e-10;
    }
    std::transform(V.begin(), V.end(), m_data.m_V.begin(), [&](const Vec &v) {
        auto gv = get_vertex(v);

        if (threshold) {
            // std::cout << "Thresholding!" << std::endl;
            // std::cout << std::string(gv) << " => ";
            gv.apply_thresholding(mythresh);
            /*
                           for(int i = 0; i < D; ++i) {
                           bool higher = gv.quot(i) > .5;
                           double p = (gv.coord[i] + (higher ? 1 : 0)) *
               vertex_grid().dx()(i) + vertex_grid().origin()(i);
                           if(std::abs(v(i) - p) < threshold_epsilon)
                           {
                           gv.quot(i) = 0;
                           gv.clamped_indices[i] = true;
                           if(higher) {
                           gv.coord[i]++;
                           }

                           }
                           }
                           */
            // std::cout << std::string(gv )<< std::endl;
        }

        return gv;
    });
    m_data.update_topology_masks();
}
template <int D>
void CutCellEdgeGenerator<D>::update_grid(const StaggeredGrid &sg) {
    StaggeredGrid::operator=(sg);
    m_data.update_grid(sg);
    m_active_grid_cell_mask.resize(sg.cell_shape());
}
template <int D>
void CutCellEdgeGenerator<D>::update_vertices(
    const ColVecs &V, const std::optional<double> &threshold) {
    update_vertices(colvecs_to_vecvector(V), threshold);
}
template <int D>
CutCellEdgeGenerator<D>::CutCellEdgeGenerator(const StaggeredGrid &g)
    : CutCellEdgeGenerator(VecVector{}, g.vertex_shape()) {}
/*
           template <int D>
           auto CutCellEdgeGenerator<D>::get_crossings(int eidx, int edge_index)
   const -> std::set<Intersection>{ return
   get_crossings(get_orig_line(eidx),edge_index);
           }
           */

template <int D>
void CutCellEdgeGenerator<D>::bake() {
    { m_data.bake(vertex_grid(), true, toss_external_facets); }

    {
        // auto t = mtao::logging::timer("generator bake vertices");
        auto t =
            mtao::logging::profiler("grid bake vertices", false, "profiler");
        bake_vertices();
        mtao::logging::trace() << "Number of crossings: " << m_crossings.size();
#if !defined(NDEBUG)
        int vc = 0;
        int ec = 0;
        int fc = 0;
        for (auto &&c : m_crossings) {
            using VType = Vertex<D>;
            using EdgeIsect = EdgeIntersection<D>;
            using TriIsect = TriangleIntersection<D>;
            std::visit(
                [&](auto &&v) {
                    using T = typename std::decay_t<decltype(v)>;
                    if constexpr (std::is_same_v<T, VType const *>) {
                        vc++;
                    } else if constexpr (std::is_same_v<T, EdgeIsect const *>) {
                        ec++;
                    } else if constexpr (std::is_same_v<T, TriIsect const *>) {
                        fc++;
                    }
                },
                c.vv);
        }
        // mtao::logging::trace() << "Number of crossings: " <<
        // m_crossings.size();
        spdlog::warn("Number of crossings {}: (vc:{},ec:{},fc:{})",
                     m_crossings.size(), vc, ec, fc);
#endif
    }
    {
        // auto t = mtao::logging::timer("generator bake edges");
        auto t = mtao::logging::profiler("grid active grid cells", false,
                                         "profiler");
        bake_active_grid_cell_mask();
    }
    {
        // auto t = mtao::logging::timer("generator bake edges");
        auto t = mtao::logging::profiler("grid bake edges", false, "profiler");
        bake_edges();
        mtao::logging::debug()
            << "Number of edges [cut,grid]: [" << m_cut_edges.size() << ","
            << m_grid_edges.size() << "]";
    }
    {
        // auto t = mtao::logging::timer("generator bake faces");
        auto t = mtao::logging::profiler("grid bake faces", false, "profiler");
        bake_faces();
    }
}
template <int D>
void CutCellEdgeGenerator<D>::bake_vertices() {
    const std::vector<CrossingType> &crossings = data().crossings();

    for (auto &&c : crossings) {
        // std::cout << std::string(c) << std::endl;
    }
    int maxind = -1;
    for (auto &&c : crossings) {
        maxind = std::max(c.index, maxind);
    }
    maxind++;
    maxind -= data().grid_size();
    maxind = std::max<int>(maxind, 0);
    // maxind -= data.nV();
    // std::cout << "Max ind: " << maxind << std::endl;
    m_crossings.resize(maxind);
    m_newV.resize(maxind);
    int gsize = StaggeredGrid::vertex_size();
    for (auto &&c : crossings) {
        int newidx = c.index - gsize;
        if (newidx >= 0) {
            m_crossings[newidx] = c;
            ////std::cout << std::string(c) << " => " << std::string(c2) <<
            /// std::endl;
            m_newV[newidx] = get_world_vertex(c.vertex());
        }
    }
}
template <int D>
void CutCellEdgeGenerator<D>::add_boundary_elements(const BoundaryElements &F) {
    auto t = mtao::logging::profiler("creating facets", false, "profiler");
    m_data.set_topology(F);
}
template <int D>
void CutCellEdgeGenerator<D>::set_boundary_elements(const BoundaryElements &F) {
    auto t = mtao::logging::profiler("creating facets", false, "profiler");
    m_data.reset_topology();
    m_data.set_topology(F);
}
template <int D>
void CutCellEdgeGenerator<D>::add_edges(const Edges &E) {
    // auto t = mtao::logging::timer("Adding edge");
    mtao::map<Edge, int> newEMap;
    auto Esize = [&]() { return m_origEMap.size() + newEMap.size(); };
    auto has_e = [&](Edge e) -> bool {
        // if we can't find it in this order
        if (!m_origEMap.contains(e) && !newEMap.contains(e)) {
            std::swap(e[0], e[1]);
            // or in the reverse order
            if (!m_origEMap.contains(e) && !newEMap.contains(e)) {
                return false;
            }
        }
        // then this must be a new one
        return true;
    };
    int count = 0;
    for (int i = 0; i < E.cols(); ++i) {
        Edge e;
        balsa::eigen::stl2eigen(e) = E.col(i);
        if (!has_e(e)) {
            // int oldsize = newEMap.size();
            // spdlog::info("{} <= {}", fmt::join(e, ","), Esize());
            newEMap[e] = Esize();
            // int newsize = newEMap.size();
            // if (oldsize == newsize) {
            //    for (int j = 0; j < i; ++j) {
            //        Edge e2;
            //        balsa::eigen::stl2eigen(e2) = E.col(j);
            //    }
            //}

            ////mtao::logging::debug() << i << " " << newEMap[e];
        } else {
            // std::cout << "Existing edge!" << i << ") " <<
            // E.col(i).transpose() << std::endl;
        }
    }
    Edges origE;
    int oldsize = origE.cols();
    origE.resize(2, Esize());
    for (auto &&[e, idx] : newEMap) {
        origE.col(idx) = balsa::eigen::stl2eigen(e);
    }
    std::copy(newEMap.begin(), newEMap.end(),
              std::inserter(m_origEMap, m_origEMap.end()));
    m_data.set_topology(origE);
}
template <int D>
void CutCellEdgeGenerator<D>::update_vertices_from_intersections() {
    m_newV.clear();
    m_newV.resize(crossings().size());
    for (auto &&c : crossings()) {
        m_newV[c.index - grid_vertex_size()] =
            get_world_vertex(c.grid_vertex());
    }
}
template <int D>
auto CutCellEdgeGenerator<D>::cell_edge_intersections(
    const std::set<coord_type> &cells) -> std::set<EdgeIntersectionType> {
    std::set<EdgeIntersectionType> isects;
    for (auto &&c : cells) {
        per_cell_vertex_looper(
            [&](const coord_type &vc, const std::bitset<D> &bs) {
                isects.emplace(EdgeIntersectionType{vc, Vec::Zero(), -1, 0});
            },
            c);
    }
    return isects;
}

template <int D>
auto CutCellEdgeGenerator<D>::grid_coord(const Vec &v) const
    -> std::tuple<coord_type, std::array<double, D>> {
    auto CQ = StaggeredGrid::coord(v);
    return CQ;
}

template <int D>
auto CutCellEdgeGenerator<D>::active_cells() const -> std::set<coord_type> {
    std::set<coord_type> cells;
    auto add_cells = [&](const coord_type &c, const Vec &quot) {
        std::set<coord_type> mycells;
        mycells.insert(c);

        coord_type myc = c;
        std::function<void(bool, int)> mac;
        mac = [&](bool sign, int dim) {
            myc[dim] = c[dim] - sign;
            if (dim == D - 1) {
                mycells.insert(myc);
            } else {
                if (quot(dim + 1) == 0) {
                    mac(0, dim + 1);
                    mac(1, dim + 1);
                } else {
                    mac(0, dim + 1);
                }
            }
        };

        /*
                       for(int i = 0; i < D; ++i) {
                       if(quot(i) == 0) {
                       std::set<coord_type> newcoords;
                       std::transform(mycells.begin(),mycells.end(),std::inserter(newcoords,newcoords.end()),[&](coord_type
           c) { c[i]--; return c;
                       });
                       mycells.insert(newcoords.begin(),newcoords.end());
                       }
                       }
                       */
        if (quot(0) == 0) {
            mac(0, 0);
            mac(1, 0);
        } else {
            mac(0, 0);
        }
        cells.insert(mycells.begin(), mycells.end());
    };

    for (auto &&c : data().crossings()) {
        auto &gv = c.vertex();
        add_cells(gv.coord, gv.quot);
    }
    /*
                   for(auto&& p: data().V()) {
                   add_cells(p.coord,p.quot);
                   }
                   */
    return cells;
}
template <int D>
auto CutCellEdgeGenerator<D>::crossing_indices() const
    -> std::array<mtao::map<coord_type, std::set<int>>, D> {
    std::array<mtao::map<coord_type, std::set<int>>, D> crossings;

    for (auto &&[idx, isect] : mtao::iterator::enumerate(m_crossings)) {
        for (int d = 0; d < D; ++d) {
            if (isect.quot(d) == 0) {
                crossings[d][isect.coord].insert(idx);
            }
        }
    }
    return crossings;
}

/*
           template <int D, typename Func>
           static void CutCellEdgeGenerator<D>::per_cell_vertex_looper(Func&& f,
   const coord_type& c) { for(int i = 0; i < (2 << D); ++i) { coord_type cc = c;
           std::bitset<D> bs(i);
           for(int j = 0; j < D; ++j) {
           cc[j] += bs[j]?1:0;
           }
           f(cc);

           }
           };
           */
template <int D>
mtao::map<int, int> CutCellEdgeGenerator<D>::vertex_reindexer() const {
    // TODO: this only handles vertex-grid vertex checks. should i check for
    // vertex-edge checks?
    const int grid_size = grid_vertex_size();
    int num_verts = grid_size;  // initialize to grid size
    mtao::map<int, int> vert_reindex;
    for (int i = 0; i < grid_size; ++i) {
        vert_reindex[i] = i;
    }
    auto reindexvertex = [&](coord_type c, Vec q, int oldidx) {
        // for(int i = 0; i < D; ++i) {
        //    if(q(i) > .5) {
        //        c[i]+=1;
        //        q(i) = q(i) - 1;
        //    }
        //}
        // if its a new one then its going to go at "verts.size"
        int newidx = (q.template lpNorm<Eigen::Infinity>() == 0)
                         ? vertex_index(c)
                         : num_verts;
        vert_reindex[oldidx] = newidx;
        if (newidx >= grid_size) {
            num_verts++;
        }
    };
    /*
                   size_t offset = origV().size();
                   for(auto&& [i,oV]: mtao::iterator::enumerate(data.V())) {
                   reindexvertex(oV.coord,oV.quot,i+grid_size);
                   }
                   */
    // assert(intersections().size()  == newV().size());
    for (auto &&[i, c] : mtao::iterator::enumerate(m_crossings)) {
        auto &isect = c.vertex();
        reindexvertex(isect.coord, isect.quot, c.index);
    }
    return vert_reindex;
}
template <int D>
auto CutCellEdgeGenerator<D>::vertex(int i) const -> Vec {
    size_t grid_size = grid_vertex_size();
    size_t osize = origV().size();
    size_t nsize = newV().size();

    assert(i >= 0);
    if (i < grid_size) {
        return StaggeredGrid::vertex(i);
    } else if (size_t gnsize = grid_size + osize; i < gnsize) {
        return origV(i - grid_size);
    } else if (size_t size = grid_size + osize + nsize; i < size) {
        return newV(i - gnsize);
    } else {
        assert(i <= size);
    }
    return {};
}
template <int D>
size_t CutCellEdgeGenerator<D>::num_vertices() const {
    return total_vertex_size();
}
template <int D>
auto CutCellEdgeGenerator<D>::compact_vertices() const -> ColVecs {
    auto vri = vertex_reindexer();
    int size = 0;
    for (auto [a, b] : vri) {
        size = std::max(size, b);
    }
    size++;
    size -= grid_vertex_size();
    ColVecs V(D, size);
    for (auto [a, b] : vri) {
        if (b >= grid_vertex_size()) {
            V.col(b - grid_vertex_size()) = vertex(a);
        }
    }
    return V;
}

/*
           template <int D>
        //Pupulate the mask and add RBD cut edges to the mesh
        auto  CutCellEdgeGenerator<D>::cut_vertex_metas() const ->
   CutVertexMetas { int grid_size = vertex_size(); CutVertexMetas
   cvm(*this,m_intersections); return cvm;
        }
        */
/*
           template <int D>
           template <typename ReidxFunc>
           void
   CutCellEdgeGenerator<D>::paci_reindex(std::array<crossing_store_type,D>&
   paci, ReidxFunc&& f) const { for(auto&& pac: paci) { for(auto&& [k,cs]: pac)
   { std::set<Crossing> ncs;
           std::transform(cs.begin(),cs.end(),std::inserter(ncs,ncs.end()),[&](const
   Crossing& c) { Crossing cn = c; int& index = cn.index; index= f(index);
           return cn;

           });
           cs = std::move(ncs);
           }
           }
           }
           */
template <int D>
// Pupulate the mask and add RBD cut edges to the mesh
auto CutCellEdgeGenerator<D>::get_per_axis_crossing_indices() const
    -> std::array<crossing_store_type, D> {
    auto paci = get_per_axis_crossing_indices(data().crossings());
    return paci;
}

template <int D>
// Pupulate the mask and add RBD cut edges to the mesh
auto CutCellEdgeGenerator<D>::get_per_axis_crossing_indices(
    const mtao::vector<Crossing<D>> &crossings) const
    -> std::array<crossing_store_type, D> {
    // Three types of vertices: original vertices, intersection vertices,
    // stencil vertices that are neither of the above The first two can be
    // stored as extraneous vertices, the last one will be vvirtual vertices

    // auto t = mtao::logging::timer("per axis crossing indices");
    std::array<crossing_store_type, D> per_axis;

    for (auto &&crossing : crossings) {
        auto &gv = crossing.vertex();

        if (gv.clamped_indices.count() == D - 1) {
            for (int j = 0; j < D; ++j) {
                if (!gv.clamped(j)) {
                    per_axis[j][gv.coord].emplace(crossing);
                    for (auto &&crossing : per_axis[j][gv.coord]) {
                    }
                    break;
                }
            }
        }
    }
    return per_axis;
}
template <int D>
auto CutCellEdgeGenerator<D>::stl_V() const -> VecVector {
    VecVector R = origV();
    const auto &n = newV();
    std::copy(n.begin(), n.end(), std::back_inserter(R));
    return R;
}
template <int D>
auto CutCellEdgeGenerator<D>::V() const -> ColVecs {
    // auto sV = stl_V();
    auto sV = newV();
    ColVecs R(D, sV.size());
    for (int i = 0; i < sV.size(); ++i) {
        R.col(i) = sV[i];
    }
    return R;
}
template <int D>
auto CutCellEdgeGenerator<D>::all_V() const -> ColVecs {
    return balsa::eigen::hstack(StaggeredGrid::vertices(), V());
}
template <int D>
auto CutCellEdgeGenerator<D>::all_GV() const -> ColVecs {
    ColVecs R(D, num_vertices());
    for (int i = 0; i < num_vertices(); ++i) {
        R.col(i) = grid_vertex(i).p();
    }
    return R;
}

template <int D>
auto CutCellEdgeGenerator<D>::V(int i) const -> Vec {
    if (i < grid_vertex_size()) {
        return StaggeredGrid::vertex(i);
    } else {
        return newV(i - grid_vertex_size());
    }
}
template <int D>
auto CutCellEdgeGenerator<D>::GV(int i) const -> VType {
    return grid_vertex(i);
}
template <int D>
void CutCellEdgeGenerator<D>::extra_metadata(CutCellMesh<D> &mesh) const {
    auto reindexer = vertex_reindexer();
    {
        size_t grid_size = grid_vertex_size();
        size_t osize = origV().size();
        size_t noff = osize + grid_size;
        size_t nsize = newV().size();
        size_t end = noff + nsize;

        for (int i : mtao::iterator::range(grid_size, noff)) {
            mesh.original_vertices.insert(reindexer[i]);
        }
        for (int i : mtao::iterator::range(noff, end)) {
            mesh.new_vertices.insert(reindexer[i]);
        }
    }
}

template <int D, typename PACI>
void show_crossings(const PACI &paci) {
    auto show_crossing = [&](auto &&crossing) {
        for (auto &&[ij, crossings] : crossing) {
            std::cout << "(";
            std::copy(ij.begin(), ij.end(),
                      std::ostream_iterator<int>(std::cout, ","));
            std::cout << "):";
            for (auto c : crossings) {
                std::cout << "{";
                std::copy(ij.begin(), ij.end(),
                          std::ostream_iterator<int>(std::cout, ","));
                std::cout << "}";
                std::cout << "[" << std::string(c.vertex()) << "=" << c.index
                          << "]";
            }
            std::cout << std::endl;
        }
    };
    for (int i = 0; i < D; ++i) {
        char grid = 'U' + i;
        // std::cout << grid << "\n====" << std::endl;show_crossing(paci[i]);
    }
}
template <int D, typename PACI, typename GridType>
balsa::eigen::ColVectors<double, D> paci_verts(const PACI &paci, const GridType &g) {
    int count = 0;
    size_t grid_size = g.size();
    for (auto &&cs : paci) {
        for (auto &&[coord, crossings] : cs) {
            for (auto &&c : crossings) {
                if (c.index >= grid_size) {
                    count++;
                }
            }
            count += crossings.size();
        }
    }
    balsa::eigen::ColVectors<double, D> V(D, count);
    int i = 0;
    for (auto &&cs : paci) {
        for (auto &&[coord, crossings] : cs) {
            for (auto &&c : crossings) {
                if (c.index >= grid_size) {
                    auto v = g.origin() + g.dx().asDiagonal() * c.vertex().p();
                    V.col(i++) = v;
                }
            }
        }
    }
    return V;
}
template <int D>
void CutCellEdgeGenerator<D>::clear() {
    m_data.clear();
    m_origEMap.clear();
    m_per_axis_crossings = {};
    m_newV.clear();
    for (auto &&a : m_per_axis_crossings) {
        a.clear();
    }
    m_crossings.clear();
    cut_edges = {};
    m_grid_edges.clear();
    m_cut_edges.clear();
    m_cut_faces.clear();
}

template <int D>
void CutCellEdgeGenerator<D>::reset(const mtao::vector<VType> &gvs) {
    m_data = CutData<D>(gvs);
    m_origV.clear();
    m_origEMap.clear();
    for (auto &&a : m_per_axis_crossings) {
        a.clear();
    }
    cut_edges = {};
    m_newV.clear();
    m_crossings.clear();
    m_grid_edges.clear();
    m_cut_edges.clear();
    m_cut_faces.clear();
    std::transform(gvs.begin(), gvs.end(), std::back_inserter(m_origV),
                   [&](const VType &gv) { return get_world_vertex(gv); });
}

/*
           template <int D>
           mtao::StackedReIndexer<2>
   CutCellEdgeGenerator<D>::prune_intersections() { mtao::StackedReIndexer<2>
   reindexer; mtao::map<Vertex,int> gvs; auto check = [&](const Vertex& gv, int
   idx) { if(auto it = gvs.find(isect); it != gvs.end()) { return it->second; }
   else { gvs[isect] = idx; return idx;
           }
           };
           mtao::geometry::grid::utils::multi_loop(vertex_shape(), [&](auto&&
   coord) { int vi = vertex_index(coord); check(Vertex(coord),vi);
           });
           int grid_size = vertex_size();
           std::vector<Vertex>
           for(auto&& gv: m_origGV) {
           check(
           }
           int offset = grid_size+m_origGV.size();
           for(auto&& [idx, isect]: mtao::iterator::enumerate(m_intersectiosn))
   { check(isect,idx+offset);
           }
           mtao::vector<Intersection> new_isects(reindexer.size());
           for(auto&& [idx, newidx]:
   mtao::iterator::enumerate(reindexer.inverse())) { new_isects[newidx] =
   m_intersections[idx];
           }
           m_intersections = std::move(new_isects);

           return reindexer;
           }
           */

template <int D>
void CutCellEdgeGenerator<D>::bake_active_grid_cell_mask() {
    auto s = StaggeredGrid::cell_shape();
    m_active_grid_cell_mask = GridDatab::Constant(true, s);

    mtao::map<std::array<int, 2>, int> cutedge_to_edge_map;

    m_per_axis_crossings = get_per_axis_crossing_indices();

    GridDatab vertmask = GridDatab::Constant(true, vertex_shape());
    {
        auto t = mtao::logging::timer("Adding grid vertices");
        auto add_grid = [&](int dim, const coord_type &c) {
            m_per_axis_crossings[dim][c];
            /*
                           per_boundary_cell_vertex_looper(dim,[&](auto&& c,
               auto&& bs) { Vec v = Vec::Zero();
                           }, c);
                           */
        };

        auto AC = active_cells();

        tbb::parallel_for(int(0), D, [&](int i) {
            for (auto &&c : AC) {
                if (m_active_grid_cell_mask.valid_index(c)) {
                    m_active_grid_cell_mask(c) = false;

                    per_boundary_cell_vertex_looper(
                        i,
                        [&](const coord_type &a, const std::bitset<D> &bs) {
                            add_grid(i, a);
                        },
                        c);
                }
            }
        });
    }
}
template <int D>
void CutCellEdgeGenerator<D>::bake_edges() {
    {
        auto t = mtao::logging::timer("Adding data edges");
        m_cut_edges = data().cut_edges();
    }

    // spdlog::info("Bad edges from data directly");
    for (auto &&[idx, e] : mtao::iterator::enumerate(m_cut_edges)) {
        if (e.indices[0] >= total_vertex_size() ||
            e.indices[1] >= total_vertex_size()) {
            // std::cout << idx << ": " << std::string(e.mask()) << " ";
            // std::cout << e.indices[0] << ":" << e.indices[1] << std::endl;
        }
    }
    // spdlog::info("====================directly");

    auto get_grid_edge = [&](coord_type c, int type) -> Edge {
        int fidx = StaggeredGrid::vertex_index(c);
        c[type]++;
        int nidx = StaggeredGrid::vertex_index(c);
        return Edge{fidx, nidx};
    };

    {
        std::mutex edges_mutex;
        // Create every new edge possible where at least one side is on the
        // interior of the mask
        auto add_edges = [&](auto &&s, const coord_type &grid_edge, int dim) {
            std::vector<int> vec;
            mtao::map<double, int> edge_ordered_elements;
            auto e = get_grid_edge(grid_edge, dim);
            edge_ordered_elements[0] = e[0];
            edge_ordered_elements[1] = e[1];
            coord_mask<D> mask(grid_edge);
            mask[dim].reset();
            std::transform(s.begin(), s.end(),
                           std::inserter(edge_ordered_elements,
                                         edge_ordered_elements.end()),
                           [&](const EdgeCrossing<D> &c) {
                               return std::make_pair(c.edge_coord(), c.index);
                           });

            // grab the indices from sorted intervals
            std::transform(edge_ordered_elements.begin(),
                           edge_ordered_elements.end(), std::back_inserter(vec),
                           [](auto &&pr) { return std::get<1>(pr); });

            std::scoped_lock lock(edges_mutex);
            for (int i = 0; i < vec.size() - 1; ++i) {
                auto e = std::array<int, 2>{{vec[i], vec[i + 1]}};
                if (e[0] != e[1]) {
                    m_grid_edges.emplace_back(CoordMaskedEdge<D>{mask, e});
                }
            }
        };

        tbb::parallel_for(size_t(0), m_per_axis_crossings.size(), [&](int i) {
            // for(auto&& [i,crossings]:
            // mtao::iterator::enumerate(m_per_axis_crossings)) {
            auto &&crossings = m_per_axis_crossings[i];
            auto t = mtao::logging::timer("Per axis crossings timing:");

            for (auto &&[c, crossings] : crossings) {
                if (StaggeredGrid::template grid<1>(i).valid_index(c)) {
                    add_edges(crossings, c, i);
                }
            }
        });
    }
}
template <int D>
CutCellMesh<D> CutCellEdgeGenerator<D>::generate_vertices() const {
    // auto t = mtao::logging::timer("Generating edges");
    int size = 0;
    // get the number of non-grid vertices
    for (auto &&c : crossings()) {
        int idx = c.index;
        size = std::max(size, idx + 1);
    }
    size -= grid_vertex_size();
    // if we only have grid vertices then size was a grtid vertex, so
    // subtracting made it negative
    size = std::max<int>(0, size);

    // for every vertex that isn't a grid vertex lets extract it
    std::vector<VType> V(size);
    for (auto &&c : crossings()) {
        if (!is_grid_vertex(c.index)) {
            int idx = c.index - grid_vertex_size();
            if (idx >= 0 && idx < size) {
                V[idx] = c.vertex();
            } else {
                spdlog::error(
                    "More non-grid cut-vertices than there should be ({} / {})",
                    idx, size);
            }
        }
    }
    CutCellMesh<D> ret(*this, V);
    return ret;
}
template <int D>
CutCellMesh<D> CutCellEdgeGenerator<D>::generate_edges() const {
    CutCellMesh<D> ret = generate_vertices();
    ret.m_cut_edges.reserve(m_cut_edges.size() + m_grid_edges.size());

    if constexpr (D == 2) {
        balsa::eigen::ColVectors<double, D> N(D, data().nE());
        for (int i = 0; i < N.cols(); ++i) {
            auto n = N.col(i);
            auto e = data().E(i);
            auto v0 = origV(e(0));
            auto v1 = origV(e(1));
            balsa::eigen::Vec2d r = v1 - v0;
            n.x() = -r.y();
            n.y() = r.x();
        }

        std::transform(m_cut_edges.begin(), m_cut_edges.end(),
                       std::back_inserter(ret.m_cut_edges),
                       [&](const CutMeshEdge<2> &E) {
                           return CutEdge<2>(E, N.col(E.parent_eid));
                       });
        std::copy(m_grid_edges.begin(), m_grid_edges.end(),
                  std::back_inserter(ret.m_cut_edges));

    } else {
        // make sure that hybrid cut_edges havea nonempty masks (because they
        // must have some mask). something. It could be a grid-cut-edge and TODO
        // lets make sure that it's identified properly if it is
        std::transform(m_cut_edges.begin(), m_cut_edges.end(),
                       std::back_inserter(ret.m_cut_edges),
                       [](const CutMeshEdge<D> &E) -> CutEdge<D> {
                           if (E.parent_eid == -1) {
                               auto m = E.mask();
                               for (int j = 0; j < D; ++j) {
                                   if (m[j]) {
                                       return CutEdge<D>(
                                           m, E.indices,
                                           std::array<int, 2>{{j, *m[j]}});
                                   }
                               }
                               // we can only get here if an edge without a
                               // parent eid does ont belong in a grid plane,
                               // which would be teh result of some very weird
                               // numerical artifact.
                               assert(false);
                           }
                           return E;
                       });
        std::copy(m_grid_edges.begin(), m_grid_edges.end(),
                  std::back_inserter(ret.m_cut_edges));
        // spdlog::error("Generator mesh edges: ");
        // for (auto &&[idx, e] : mtao::iterator::enumerate(m_cut_edges)) {
        //    std::cout << idx << ": " << std::string(e.mask()) << " ";
        //    std::cout << e.indices[0] << ":" << e.indices[1] << std::endl;
        //}
        // spdlog::error("Generator grid edges: ");
        // for (auto &&[idx, e] : mtao::iterator::enumerate(m_grid_edges)) {
        //    std::cout << idx << ": " << std::string(e.mask()) << " ";
        //    std::cout << e.indices[0] << ":" << e.indices[1] << std::endl;
        //}
        // spdlog::error("Generator cut edges: ");
        // for (auto &&[idx, e] : mtao::iterator::enumerate(ret.m_cut_edges)) {
        //    std::cout << idx << ": " << std::string(e.mask()) << " "
        //              << std::string(e) << std::endl;
        //}
    }

    for (auto &&[eidx, ce] : mtao::iterator::enumerate(m_cut_edges)) {
        if (ce.parent_eid >= 0) {
            balsa::eigen::Vec2d ts;
            for (auto &&[i, vi, t] : mtao::iterator::enumerate(
                     ce.indices,
                     mtao::iterator::shell(ts.data(), ts.data() + ts.size()))) {
                if (vi < grid_vertex_size()) {
                    t = data().get_edge_coord(ce.parent_eid, grid_vertex(vi));
                } else {
                    // spdlog::info("Looking into crossing {} of {}", vi,
                    //             total_vertex_size());
                    // std::cout << std::string(crossing(vi)) << std::endl;
                    t = data().get_edge_coord(ce.parent_eid, crossing(vi));
                }
            }
            ret.m_mesh_cut_edges[eidx] = InterpolatedEdge{ts, ce.parent_eid};
        }
    }
    return ret;
}
template <int D>
auto CutCellEdgeGenerator<D>::edges() const -> std::set<Edge> {
    // auto t = mtao::logging::timer("Generating edges");
    std::set<Edge> ret;
    std::transform(m_cut_edges.begin(), m_cut_edges.end(),
                   std::inserter(ret, ret.end()),
                   [](auto &&e) { return e.indices; });
    std::transform(m_grid_edges.begin(), m_grid_edges.end(),
                   std::inserter(ret, ret.end()),
                   [](auto &&e) { return e.indices; });
    return ret;
}
template <int D>
auto CutCellEdgeGenerator<D>::edge_slice(int dim, int coord) const
    -> std::set<Edge> {
    // auto t = mtao::logging::timer("Generating edges");
    std::set<Edge> ret;
    for (auto &&e : m_grid_edges) {
        if (e[dim] && *e[dim] == coord) {
            ret.emplace(e.indices);
        }
    }
    for (auto &&e : m_cut_edges) {
        if (e[dim] && *e[dim] == coord) {
            ret.emplace(e.indices);
        }
    }
    return ret;
}
template <int D>
auto CutCellEdgeGenerator<D>::axial_edges() const
    -> std::array<mtao::map<int, std::set<Edge>>, D> {
    std::array<mtao::map<int, std::set<Edge>>, D> ret;
    for (auto &&e : m_cut_edges) {
        for (int dim = 0; dim < D; ++dim) {
            if (e[dim]) {
                ret[dim][*e[dim]].emplace(e.indices);
            }
        }
    }
    for (auto &&e : m_grid_edges) {
        for (int dim = 0; dim < D; ++dim) {
            if (e[dim]) {
                ret[dim][*e[dim]].emplace(e.indices);
            }
        }
    }
    return ret;
}
template <int D>
auto CutCellEdgeGenerator<D>::vertex_ownership() const
    -> mtao::map<coord_type, std::set<int>> {
    mtao::map<coord_type, std::set<int>> crossings;

    for (auto &&[idx, isect] : mtao::iterator::enumerate(m_crossings)) {
        auto bs = isect.vertex().clamped_indices.as_bitset1();
        coord_type c = isect.coord;
        per_dual_cell_vertex_looper(
            [&, idx = idx](const coord_type &c, std::bitset<D> &b) {
                if ((b & bs) == b) {
                    crossings[c].insert(idx);
                }
            },
            isect.coord);
    }
    return crossings;
}
template <int D>
auto CutCellEdgeGenerator<D>::possible_cells(const std::vector<int> &face) const
    -> std::set<coord_type> {
    if (face.empty()) {
        return {};
    }
    std::set<coord_type> possibles = GV(face[0]).possible_cells();

    for (auto &&f : face) {
        auto s = GV(f).possible_cells();
        std::set<coord_type> i;
        std::set_intersection(possibles.begin(), possibles.end(), s.begin(),
                              s.end(), std::inserter(i, i.end()));
        possibles = std::move(i);
        if (possibles.empty()) {
            return {};
        }
    }

    return possibles;
}
template <int D>
auto CutCellEdgeGenerator<D>::possible_cells(
    const std::set<std::vector<int>> &faces) const -> std::set<coord_type> {
    if (faces.empty()) {
        return {};
    }
    std::set<coord_type> possibles = GV((*faces.begin())[0]).possible_cells();

    for (auto &&face : faces) {
        for (auto &&f : face) {
            auto s = GV(f).possible_cells();
            std::set<coord_type> i;
            std::set_intersection(possibles.begin(), possibles.end(), s.begin(),
                                  s.end(), std::inserter(i, i.end()));
            possibles = std::move(i);
            if (possibles.empty()) {
                return {};
            }
        }
    }

    return possibles;
}
template <int D>
auto CutCellEdgeGenerator<D>::possible_cells_cell(
    const std::set<int> &faces, const std::vector<CutFace<D>> &CFs) const
    -> std::set<coord_type> {
    if (faces.empty()) {
        return {};
    }
    std::set<coord_type> possibles =
        possible_cells(CFs[*faces.begin()].indices);

    for (auto &&face : faces) {
        auto s = possible_cells(CFs[face].indices);
        std::set<coord_type> i;
        std::set_intersection(possibles.begin(), possibles.end(), s.begin(),
                              s.end(), std::inserter(i, i.end()));
        possibles = std::move(i);
        if (possibles.empty()) {
            spdlog::warn("No possible cell!");
            return {};
        }
    }

    return possibles;
}
template <int D>
auto CutCellEdgeGenerator<D>::face_mask(
    const std::set<std::vector<int>> &faces) const -> coord_mask<D> {
    if (faces.empty()) {
        return {};
    }
    auto mask = face_mask(*faces.begin());
    auto it = faces.begin();
    ++it;
    for (; it != faces.end(); ++it) {
        mask &= face_mask(*it);
    }
    return mask;
}
template <int D>
auto CutCellEdgeGenerator<D>::face_mask(const std::vector<int> &faces) const
    -> coord_mask<D> {
    if (faces.empty()) {
        return {};
    }
    auto mask = get_mask(faces[0]);
    auto it = faces.begin();
    ++it;
    for (; it != faces.end(); ++it) {
        mask &= get_mask(*it);
    }
    return mask;
}
/*
                   template <int D>
                   auto CutCellEdgeGenerator<D>::possible_faces(const
   std::set<std::vector<int>>& faces) const -> std::set<coord_type> {

                   if(faces.empty()) { return {}; }
                   std::set<coord_type> possibles =
   GV((*faces.begin())[0]).possible_faces();

                   for(auto&& face: faces) {
                   for(auto&& f: face) {
                   auto s = GV(f).possible_faces();
                   std::set<coord_type> i;
                   std::set_intersection(possibles.begin(),possibles.end(),s.begin(),s.end(),std::inserter(i,i.end()));
                   possibles = std::move(i);
                   if(possibles.empty()) {
                   return {};
                   }
                   }
                   }

                   return possibles;
                   }
                   template <int D>
                   auto CutCellEdgeGenerator<D>::possible_faces(const
   std::vector<int>& face) const -> std::set<coord_type> {

                   if(face.empty()) { return {}; }
                   std::set<coord_type> possibles =
   GV(face[0]).possible_faces();

                   for(auto&& f: face) {
                   auto s = GV(f).possible_faces();
                   std::set<coord_type> i;
                   std::set_intersection(possibles.begin(),possibles.end(),s.begin(),s.end(),std::inserter(i,i.end()));
                   possibles = std::move(i);
                   if(possibles.empty()) {
                   return {};
                   }
                   }

                   return possibles;
                   }
                   */
template <int D>
bool CutCellEdgeGenerator<D>::is_in_cell(
    const std::set<std::vector<int>> &face) const {
    return !possible_cells(face).empty();
}
template <int D>
bool CutCellEdgeGenerator<D>::is_in_cell(const std::vector<int> &face) const {
    return !possible_cells(face).empty();
}

template <int D>
bool CutCellEdgeGenerator<D>::is_active_cell(const coord_type &coord) const {
    return m_active_grid_cell_mask.valid_index(coord) &&
           m_active_grid_cell_mask(coord);
}
/*
// helper for assigning boundary information to
template<int D, typename IndexContainerType>
std::optional<std::tuple<int, bool>>
  make_boundary_pair(const std::vector<Vertex<D>> &V, const
CoordMaskedGeometry<D, IndexContainerType> &boundary_facet) { auto pc =
boundary_facet.possible_cells(V); if (pc.size() != 2) { return {};
    }
    std::array<coord_type, 2> pca;
    std::copy(pc.begin(), pc.end(), pca.begin());
    //find the boundary cells axis
    int idx;
    for (idx = 0; idx < 2; ++idx) {
        if (pca[0][idx] != pca[1][idx]) {
            break;
        }
    }
    if (pca[0][idx] + 1 != pca[1][idx]) {
        std::cout << "SET WASNT LEXICOGRAPHICAL ORDER SOMEHOW?" << std::endl;
    }
    if (pca[0][idx] < 0) {
        return { -1, 1 };
    } else if (pca[1][idx] >= cell_shape()[idx]) {
        return { -1, 0 };
    } else if (I.size() == 1 && I.begin()->size() == (1 << (D - 1))) {//should
always be a grid cell! int pi = cell_index(pca[0]); int ni = cell_index(pca[1]);
        bool pa = m_active_grid_cell_mask.get(pi);
        bool na = m_active_grid_cell_mask.get(ni);
        if (na ^ pa) {//active inactive boundary
            if (pa) {
                return { pi, 1 };
            } else {
                return { ni, 0 };
            }
        }
    }
    return {};
}
*/

}  // namespace mandoline::construction
