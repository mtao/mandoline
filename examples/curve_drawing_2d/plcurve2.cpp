#include "plcurve2.hpp"
#include <iostream>
#include <mtao/logging/logger.hpp>

auto PLCurve2::lines() const -> mtao::vector<LineType> {
    mtao::vector<LineType> lines;
    for(int i = 0; i < pts.size() - !closed; ++i) {
        lines.emplace_back(pts[i],pts[(i+1)%pts.size()]);
    }
    return lines;
}
void PLCurve2::add_point(const Vec& p) {
    auto d = mtao::logging::debug();
    d << "Adding point " << p.transpose();
    pts.emplace_back(p);
    d << ". Now have " << pts.size() << " points.";
}
mtao::ColVectors<unsigned int,2> PLCurve2::edges() const {
    int size = pts.size() - !closed;
    mtao::ColVectors<unsigned int, 2> E(2,size);
    for(int i = 0; i < size; ++i) {

        E.col(i) = mtao::Vector<unsigned int,2>(i, (i+1)%pts.size());
    }
    return E;
}
auto PLCurve2::stl_points() const -> mtao::vector<Vec> {
    return pts;
}
auto PLCurve2::points() const -> ColVecs {
    ColVecs P(2,pts.size());
    for(int i = 0; i < pts.size(); ++i) {
        P.col(i) = pts[i];
    }
    return P;
}
