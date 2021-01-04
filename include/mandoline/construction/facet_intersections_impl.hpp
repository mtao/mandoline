#pragma once
#include <set>

#include "mandoline/construction/facet_intersections.hpp"

namespace mandoline::construction {
template <int D>
void EdgeIntersections<D>::bake(const std::optional<SGType> &grid) {
    intersections.clear();

    using namespace mtao::eigen;
    using Veci = mtao::Vector<int, D>;
    auto &gvstart = *vptr_edge[0];
    auto &gvend = *vptr_edge[1];

    // std::cout << "Edge intersections " ;
    // std::cout << std::string(gvstart) << " => " << std::string(gvend) <<
    // std::endl;

    mtao::eigen::stl2eigen(tangent) = gvend.p() - gvstart.p();

    // normalize and move to the upper left quadrant
    Vec direction = gvend.p() - gvstart.p();
    std::bitset<D> reflection_bs;
    for (int i = 0; i < D; ++i) {
        auto &&sc = gvstart.coord[i];
        auto &&ec = gvend.coord[i];
        auto &&sq = gvstart.quot(i);
        auto &&eq = gvend.quot(i);
        reflection_bs[i] = !((ec > sc) || (ec == sc && eq > sq));
    }

    auto reflect_axis = [&](VType &gv, int i) {
        auto &&c = gv.coord;
        auto &&q = gv.quot;
        if (reflection_bs[i]) {
            if (q[i] == 0) {
                c[i] = -c[i];
            } else if (c[i] >= 0) {
                c[i] = -(c[i] + 1);
                q(i) = 1 - q(i);
            } else {
                c[i] = -(c[i] + 1);
                q(i) = 1 - q(i);
            }
        }
    };
    auto reflect = [&](VType &gv) {
        for (int i = 0; i < D; ++i) {
            reflect_axis(gv, i);
        }
    };
    std::array<VType, 2> gvsr;

    std::transform(vptr_edge.begin(), vptr_edge.end(), gvsr.begin(),
                   [&](const VType *gvptr) {
                       VType gv = *gvptr;
                       reflect(gv);
                       return gv;
                   });

    Veci begin = stl2eigen(gvsr[0].coord) +
                 gvsr[0].quot.array().ceil().template cast<int>().matrix();
    Veci end = stl2eigen(gvsr[1].coord) +
               gvsr[1].quot.array().floor().template cast<int>().matrix();

    auto Dir = direction.array().abs().eval();

    std::set<EdgeIsect> T;
    std::array<std::map<coord_mask<D>, std::set<EdgeIsect>>, D> edges;

    std::array<coord_mask<D>, 2> vertex_masks{
        {vptr_edge[0]->mask(), vptr_edge[1]->mask()}};
    for (auto &&a : vertex_masks) {
        a -= mask();
    }

    auto invalid_isect = [&](auto &&isect) {
        auto m = isect.mask() - this->mask();
        return (vertex_masks[0] & m).active() || (vertex_masks[1] & m).active();
    };

    for (int d = 0; d < D; ++d) {
        int b = begin(d);
        int e = end(d);
        double dir = Dir(d);
        if (dir == 0) continue;
        double sq = gvsr[0].quot(d);
        double offset = std::ceil(sq) - sq;
        for (int i = 0; i <= e - b; ++i) {
            double t = (offset + i) / dir;
            if (t <= 0 || t >= 1) {
                continue;
            }

            VType np = gvstart * (1 - t) + gvend * t;
            // VType np = gvstart +  ( gvend - gvstart  ) * t;
            {
                auto M = mask();
                M.clamp(np);
            }

            int abs_grid = i + b;
            np.coord[d] = abs_grid;
            np.quot(d) = 0;
            np.clamped_indices[d] = true;
            if (reflection_bs[d]) {
                if (abs_grid <= 0) {
                    np.coord[d] = -abs_grid;
                } else {
                    np.coord[d] = -abs_grid;
                }
            }

            mask().clamp(np);

            EdgeIsect sect{np, t, edge_index};

            // sect.apply_thresholding();
            {
                auto M = mask();
                for (int i = 0; i < D; ++i) {
                    if (!M[i]) {
                        if (gvstart.coord[i] == gvend.coord[i]) {
                            np.clamped_indices[i] = {};
                        }
                    }
                }
            }
            auto mask = sect.mask();
            if (grid) {
                // auto g = grid->template grid<D-1>(D==2?(1-d):d);
                auto g = grid->cell_shape();
                bool valid = true;
                for (int i = 0; i < D; ++i) {
                    int idx = np.coord[i];
                    int q = np.quot(i);
                    if (!(i >= 0 || i < g[i] || (i == g[i] && q == 0))) {
                        valid = false;
                        break;
                    }
                }
                if (invalid_isect(sect)) {
                    valid = false;
                }
                if (valid) {
                    // if(g.valid_index(np.coord)) {
                    // if(grid->template grid<D-1>(d).valid_index(np.coord)) {
                    // if(grid->template grid<1>(d).valid_index(np.coord)) {
                    edges[mask.count() - 1][mask].emplace(sect);
                }
            } else {
                edges[mask.count() - 1][mask].emplace(sect);
            }
        }
    }

    // std::map<coord_mask<D>,std::set<EdgeIsect>> edges;
    {
        for (int i = 0; i < D - 1; ++i) {
            auto &Ei = edges[i];
            auto &Eip = edges[i + 1];
            for (auto &&[c, es] : Eip) {
                auto &cc = c;
                for (int k = 0; k < D; ++k) {
                    if (c[k]) {
                        coord_mask<D> m = c;
                        m[k] = {};
                        if (auto it = Ei.find(m); it != Ei.end()) {
                            auto &o = it->second;
                            std::transform(o.begin(), o.end(),
                                           std::inserter(es, es.end()),
                                           [&](EdgeIsect isect) -> EdgeIsect {
                                               cc.clamp(isect);
                                               return isect;
                                           });
                            Ei.erase(it);
                        }
                    }
                }
            }
        }
    }

    std::map<double, EdgeIsect> mm;
    for (auto &&ee : edges) {
        for (auto &&[m, es] : ee) {
            double lc = 0;
            const EdgeIsect *eptr = nullptr;
            for (auto &&e : es) {
                lc += e.edge_coord;
                eptr = &e;
            }
            lc /= es.size();
            auto &v = mm[lc];
            if (v.edge_index == -1) {
                v = *eptr;
                v.edge_coord = lc;
            } else {
                v.clamped_indices |= eptr->clamped_indices;
                v.repair();
            }
        }
    }

    intersections.clear();
    double mythresh = 1;
    if (grid) {
        int v = *std::max_element(grid->vertex_shape().begin(),
                                  grid->vertex_shape().end());
        v = 2 * std::max<int>(1, v);
        // mythresh =  v * threshold_epsilon;
        mythresh = v * 1e-10;
    }
    std::transform(mm.begin(), mm.end(), std::back_inserter(intersections),
                   [&](auto &&pr) {
                       EdgeIsect e = pr.second;
                       this->mask().clamp(e);
                       e.VType::apply_thresholding(mythresh);
                       return e;
                   });
    std::map<coord_mask<D>, EdgeIsect> cem;
    for (auto &&v : intersections) {
        cem[v.mask()] = v;
    }
    intersections.clear();
    std::transform(cem.begin(), cem.end(), std::back_inserter(intersections),
                   [](auto &&pr) { return std::get<1>(pr); });

    std::sort(
        intersections.begin(), intersections.end(),
        [](auto &&ei, auto &&ei2) { return ei.edge_coord < ei2.edge_coord; });

    // std::cout << intersections.size() << " found" << std::endl;
    // for(auto&& i: intersections)
    //{
    //    std::cout << std::string(i) << " ";
    //}
    // std::cout << std::endl;
}
template <int D>
void TriangleIntersections<D>::bake(const std::optional<SGType> &grid) {
    edge_intersections.clear();
    // per axis, per plane, set of intersections
    std::array<std::map<int, std::set<std::tuple<const VType *, int>>>, D> bins;
    std::map<const VType *, std::array<double, 3>> coords;
    auto bin_gv = [&](const VType &gv, int v) {
        for (int i = 0; i < D; ++i) {
            if (gv.clamped(i)) {
                auto &bin = bins[i];
                bin[gv.coord[i]].insert({&gv, v});
            }
        }
    };
    for (auto &&[bary_index, eisptr] :
         mtao::iterator::zip(bary_indices, edge_isects)) {
        for (auto &&ei : eisptr->intersections) {
            std::array<double, 3> coord;
            mtao::eigen::stl2eigen(coord) =
                edge_bary(bary_index, ei.edge_coord);
            coords[&ei] = coord;
            bin_gv(ei, eisptr->edge_index);
        }
    }
    std::set<const VType *> vertex_ptrs;
    // std::cout << "Face intersections" ;
    for (auto &&[i, gv] : mtao::iterator::enumerate(vptr_tri)) {
        //    std::cout << std::string(*gv) << ",";
        bin_gv(*gv, -i - 1);  // negative to avoid clashing with edges
        std::array<double, 3> coord;
        mtao::eigen::stl2eigen(coord) = mtao::Vec3d::Unit(i);
        coords[gv] = coord;
        vertex_ptrs.insert(gv);
    }
    // std::cout << std::endl;
    int count = 0;
    for (auto &&[dim, bin] : mtao::iterator::enumerate(bins)) {
        // std::cout << "Dim " << dim << std::endl;
        for (auto &&[coord, gvs] : bin) {
            // std::cout << "Coord " << coord <<  ": ";
            int vertex_count = 0;
            for (auto &&[v, i] : gvs) {
                // std::cout << std::string(*v) << ", ";
                if (i < 0) {
                    vertex_count++;
                }
            }
            // std::cout << std::endl;
            if (vertex_count == 2) {
                continue;
            }
            if (gvs.size() == 2) {
                count += gvs.size();
                VPtrEdge gvpe;
                std::transform(gvs.begin(), gvs.end(), gvpe.begin(),
                               [](auto &&pr) { return std::get<0>(pr); });
                // make sure we're not dealing with two vertices
                if (vertex_ptrs.find(gvpe[0]) == vertex_ptrs.end() ||
                    vertex_ptrs.find(gvpe[1]) == vertex_ptrs.end()) {
                    EdgeIntersections<D> eis(gvpe);
                    eis.bake(grid);

                    edge_intersections.emplace(gvpe, std::move(eis));
                }
            }
        }
        // std::cout << std::endl;
    }

    intersections.clear();
    intersections.reserve(vertex_size());
    std::map<const EdgeIsect *, int> edge_triangle_index_map;
    std::set<const VType *> vertices;
    std::set<const VType *> triangle_vertices;
    // std::array<std::map<coord_mask<D>,std::set<std::tuple<EdgeIsect*,std::array<double,3>>>>,3>
    // edge_isects;
    std::map<const VType *, std::array<double, 3>> barys;
    vertices.insert(vptr_tri.begin(), vptr_tri.end());
    for (auto &&eisptr : edge_isects) {
        auto &&isects = eisptr->intersections;
        for (auto &&i : isects) {
            vertices.insert(&i);
        }
    }
    for (auto &&[ge, eis] : edge_intersections) {
        auto &&isects = eis.intersections;
        for (auto &&i : isects) {
            triangle_vertices.insert(&i);
        }
    }

    for (auto &&[ge, eis] : edge_intersections) {
        auto a = mtao::eigen::stl2eigen(coords[ge[0]]);
        auto b = mtao::eigen::stl2eigen(coords[ge[1]]);
        for (auto &&ei : eis.intersections) {
            double t = ei.edge_coord;
            std::array<double, 3> varr;
            auto v = mtao::eigen::stl2eigen(varr);
            v = (1 - t) * a + t * b;
            for (int i = 0; i < 3; ++i) {
                auto &&a = v(i);
                auto &&b = v((i + 1) % 3);
                auto &&c = v((i + 2) % 3);
                if (a < 0) {
                    a = 0;
                    double d = b + c;
                    b = b / d;
                    c = c / d;
                } else if (a > 1) {
                    v = mtao::Vec3d::Unit(i);
                    break;
                }
            }
            barys[&ei] = varr;
        }
    }
    std::set<VType *> to_delete;
    auto condensed_triangle_verts = convex_grid_condensation(triangle_vertices);
    auto condensed_verts = convex_grid_condensation(vertices);

    for (auto &&verts_per_mask : condensed_triangle_verts) {
        for (auto &&[m, vs] : verts_per_mask) {
            mtao::Vec3d b = mtao::Vec3d::Zero();
            const VType *eiptr;
            for (auto &&ei : vs) {
                edge_triangle_index_map[static_cast<const EdgeIsect *>(ei)] =
                    intersections.size();
                eiptr = ei;
                b += mtao::eigen::stl2eigen(barys[ei]);
            }
            b /= vs.size();

            intersections.emplace_back(*eiptr, b, triangle_index);
        }
    }
    for (auto &&[e, i] : edge_triangle_index_map) {
        edge_to_triangle_map[e] = &intersections[i];
    }
}

template <int D>
std::array<std::map<coord_mask<D>, std::set<const Vertex<D> *>>, D + 1>
convex_grid_condensation(std::set<const Vertex<D> *> &V) {
    std::array<std::map<coord_mask<D>, std::set<const Vertex<D> *>>, D + 1> ret;

    for (auto &&v : V) {
        auto m = v->mask();
        ret[m.count()][m].insert(v);
    }
    for (int i = 0; i < D; ++i) {
        auto &Vi = ret[i];
        auto &Vip = ret[i + 1];
        for (auto &&[c, vs] : Vip) {
            for (int k = 0; k < D; ++k) {
                if (c[k]) {
                    coord_mask<D> m = c;
                    m[k] = {};
                    if (auto it = Vi.find(m); it != Vi.end()) {
                        auto &o = it->second;
                        std::transform(o.begin(), o.end(),
                                       std::inserter(vs, vs.end()),
                                       [&](auto a) {
                                           //                                c.clamp(*const_cast<Vertex<D>*>(a));
                                           return a;
                                       });
                        Vi.erase(it);
                    }
                }
            }
        }
    }
    return ret;
}

template <int D>
auto TriangleIntersections<D>::boundary_vptr_edge() const -> VPtrEdge {
    std::array<std::tuple<const EdgeIntersections<D> *, bool>, 3>
        eisptr_with_sign;
    auto &&bi = bary_indices[0];
    auto &&eisptr = edge_isects[0];
    auto [a, b] = bi;
    bool reverse = ((a + 1) % 3 == b);
    auto gv = eisptr->gvertices();
    if (reverse) {
        return VPtrEdge{{gv[1].vertex_ptr(), gv[0].vertex_ptr()}};
    } else {
        return VPtrEdge{{gv[0].vertex_ptr(), gv[1].vertex_ptr()}};
    }
}

template <int D>
std::vector<Crossing<D>> TriangleIntersections<D>::boundary_vptr_loop() const {
    std::vector<Crossing<D>> ret;
    std::array<std::tuple<const EdgeIntersections<D> *, bool>, 3>
        eisptr_with_sign;
    for (auto &&[bi, eisptr] : mtao::iterator::zip(bary_indices, edge_isects)) {
        auto [a, b] = bi;
        bool reverse = ((a + 1) % 3 == b);
        // bool reverse = ((a+1)%3 != b);
        eisptr_with_sign[2 - (reverse ? b : a)] =
            std::make_tuple(eisptr, reverse);
        // eisptr_with_sign[reverse?b:a] = std::make_tuple(eisptr,reverse);
    }
    for (auto &&[eisptr, reverse] : eisptr_with_sign) {
        auto gv = eisptr->gvertices();
        if (reverse) {
            // std::cout << std::string(gv.back()) << " => " <<
            // std::string(gv.front()) << std::endl;
            // ret.insert(ret.end(),gv.rbegin()+1,gv.rend());
            // ret.insert(ret.end(),gv.begin()+1,gv.end());
            ret.insert(ret.end(), gv.rbegin() + 1, gv.rend());
        } else {
            // std::cout << std::string(gv.front()) << " => " <<
            // std::string(gv.back()) << std::endl;
            // ret.insert(ret.end(),gv.begin()+1,gv.end());
            // ret.insert(ret.end(),gv.rbegin()+1,gv.rend());
            ret.insert(ret.end(), gv.begin() + 1, gv.end());
        }
    }

    return ret;
}

template <int D>
auto TriangleIntersections<D>::vptr_faces(
    const std::map<const VType *, int> *other_vptr_indexer) const
    -> std::set<std::vector<const VType *>> {
    std::set<std::vector<const VType *>> ret;
    /*
                               ret.emplace(std::vector<const
       VType*>(vptr_tri.begin(),vptr_tri.end())); return ret;
                               */
    auto V = gvertices();

    if (other_vptr_indexer != nullptr) {
        // spdlog::info(
        //    "Got a nonzero indexer so i might hvae to prune out my own
        //    items");
        const auto &o = *other_vptr_indexer;
        V.erase(std::remove_if(V.begin(), V.end(),
                               [&](const auto &item) {
                                   return o.find(item.vertex_ptr()) == o.end();
                               }),
                V.end());
        // spdlog::info("succeeded in cleaning myself");
    }
    populate_crossing_indices<D>(V);

    auto VM = crossing_indexer<D>(V);

    // have to use the forward indexer to produce the inverse here because there
    // are missing elements
    auto IVM = inverse_crossing_indexer<D>(V);
    auto BC = barycentric_coords();
    mtao::ColVecs2d B(2, V.size());
    for (auto &&c : V) {
        B.col(c.index) =
            mtao::eigen::stl2eigen(BC.at(c.vertex_ptr())).template topRows<2>();
    }
    // for(auto&& v: V) { std::cout << std::string(v) << " ";}
    // std::cout << std::endl;
    // std::cout << B << std::endl;
    // TODO: If
    auto Es = this->nobdry_edges(VM, false);
    // for(auto&& e: Es) {
    //    std::cout << e[0] << ":" << e[1] << " ";
    //}
    // std::cout << std::endl;
    // std::set<Edge> Es(Evec.begin(),Evec.end());
    FaceCollapser fc(Es);
    Edge be = boundary_edge(VM);
    fc.set_edge_for_removal(be);

    fc.bake(B);

    // for(auto&& [e,fc]: fc.edge_to_face()) {
    //    std::cout << e[0] << ":" << e[1] << " => " << (
    //    std::get<1>(fc)?'-':'+')<< std::get<0>(fc)  << std::endl;
    //}

    // for each face (in this case we dont have holes)
    for (auto &&[cid, fs] : fc.faces_no_holes()) {
        std::vector<const VType *> v(fs.size());
        // facecollapser does things backwards, so we gotta rbegin
        std::transform(fs.rbegin(), fs.rend(), v.rbegin(),
                       [&](int idx) { return IVM.at(idx); });
        ret.emplace(std::move(v));
    }

    return ret;
}

template <int D>
std::set<std::vector<int>> TriangleIntersections<D>::faces(
    const std::map<const VType *, int> &indexer,
    bool expect_all_vertices_in_indexer) const {
    auto ptrf = vptr_faces();
    // auto ptrf = vptr_faces(expect_all_vertices_in_indexer ? nullptr :
    // &indexer);
    std::erase_if(ptrf, [&](auto &&face_loop) {
        // for (auto &&face_loop : face_loops) {
        for (auto &&v : face_loop) {
            if (indexer.find(v) == indexer.end()) {
                return true;
            }
        }
        // }
        return false;
    });
    std::set<std::vector<int>> ret;
    std::transform(ptrf.begin(), ptrf.end(), std::inserter(ret, ret.end()),
                   [&](auto &&ptrvec) {
                       std::vector<int> vec(ptrvec.size());
                       std::transform(ptrvec.begin(), ptrvec.end(), vec.begin(),
                                      [&](const VType *v) {
                                          int val = indexer.at(v);
                                          return val;
                                      });
                       return vec;
                   });
    return ret;
}
}  // namespace mandoline::construction
