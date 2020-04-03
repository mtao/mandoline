#include "plcurve2.hpp"
#include <iostream>
#include <mtao/eigen/stl2eigen.hpp>
#include <algorithm>
#include <iterator>
#include <mtao/logging/logger.hpp>
#include <mandoline/tools/plcurve_io.hpp>
#include <mandoline/tools/edges_to_plcurves.hpp>
#include <mtao/geometry/mesh/read_obj.hpp>

PLCurve2::PLCurve2(const std::string &filename) {
    load(filename);
}

auto PLCurve2::lines() const -> mtao::vector<LineType> {
    mtao::vector<LineType> lines;
    for (auto &&[pts, closed] : curves) {
        if (!pts.empty()) {
            for (int i = 0; i < pts.size() - !closed; ++i) {
                lines.emplace_back(pts[i], pts[(i + 1) % pts.size()]);
            }
        }
    }
    return lines;
}
void PLCurve2::add_point(const Vec &p) {
    {
        auto d = mtao::logging::debug();
        d << "Adding point " << p.transpose();
        if (curves.empty()) {
            curves.emplace_back();
        }
        auto &pts = current_curve();
        pts.emplace_back(p);
        d << ". Now have " << pts.size() << " points.";
    }
    for (auto &&[c, b] : curves) {
        for (auto &&c : c) { std::cout << c.transpose() << ","; }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
mtao::ColVectors<unsigned int, 2> PLCurve2::edges() const {
    int tot_size = 0;
    for (auto &&[pts, closed] : curves) {
        if (pts.size() > 0) {
            int size = pts.size() - !closed;
            tot_size += size;
        }
    }
    mtao::ColVectors<unsigned int, 2> E(2, tot_size);
    int off = 0;
    int voff = 0;
    for (auto &&[pts, closed] : curves) {
        if (pts.size() > 0) {
            for (int i = 0; i < pts.size() - !closed; ++i) {

                E.col(off++) = mtao::Vector<unsigned int, 2>(i + voff, (i + 1) % pts.size() + voff);
            }
            voff += pts.size();
        }
    }
    return E;
}
auto PLCurve2::stl_points() const -> mtao::vector<Vec> {
    int tot_size = 0;
    for (auto &&[pts, closed] : curves) {
        int size = pts.size();
        tot_size += size;
    }
    mtao::vector<Vec> p(tot_size);
    int off = 0;
    for (auto &&[pts, closed] : curves) {
        int size = pts.size();
        if (size > 0) {
            std::copy(pts.begin(), pts.end(), p.begin() + off);
        }
        off += size;
    }

    return p;
}
auto PLCurve2::points() const -> ColVecs {
    int tot_size = 0;
    for (auto &&[pts, closed] : curves) {
        int size = pts.size();
        tot_size += size;
    }
    ColVecs P(2, tot_size);
    int off = 0;
    for (auto &&[pts, closed] : curves) {
        int size = pts.size();
        if (size > 0) {
            P.block(0, off, 2, pts.size()) = mtao::eigen::stl2eigen(pts);
        }
        off += size;
    }
    return P;
}
void PLCurve2::load(const std::string &filename) {
    mtao::ColVecs2d V;
    mtao::ColVecs2i E;
    std::tie(V,E) = mtao::geometry::mesh::read_obj2D(filename);
    if (E.minCoeff() < 0) {
        E.array() += 1;
    }
    auto c = mandoline::tools::edge_to_plcurves(V, E, true);
    std::cout << "Read " << c.size() << " plcurves" << std::endl;
    curves.clear();
    curves.resize(c.size());
    std::transform(c.begin(), c.end(), curves.begin(), [&](auto &&pr) -> std::tuple<mtao::vector<Vec>, bool> {
        auto &&[cv, closedness] = pr;
        mtao::vector<Vec> cs(cv.size());
        std::transform(cv.begin(), cv.end(), cs.begin(), [&](int idx) -> Vec {
            return V.col(idx);
        });
        return { cs, closedness };
    });
}
void PLCurve2::save(const std::string &filename) const {
    mandoline::tools::write_plcurve(points(), mtao::ColVecs2i(edges().cast<int>()), filename);
}
