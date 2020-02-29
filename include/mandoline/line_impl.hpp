#pragma once
#include "mandoline/line.hpp"
namespace mandoline {

template<int D>
Line<D>::Line(const Vec &a, const Vec &b) {
    start() = a;
    end() = b;
}

template<int D, typename Derived>
Line<D> operator+(const Eigen::MatrixBase<Derived> &T, const Line<D> &l) {
    //make sure it's a vector
    static_assert(Derived::RowsAtCompileTime == D);
    static_assert(Derived::ColsAtCompileTime == 1);
    return Line<D>(l.start_end().colwise() + T);
}
template<int D, typename Derived>
Line<D> operator+(const Line<D> &l, const Eigen::MatrixBase<Derived> &T) {
    //make sure it's a vector
    static_assert(Derived::RowsAtCompileTime == D);
    static_assert(Derived::ColsAtCompileTime == 1);
    return Line<D>(l.start_end().colwise() + T);
}

template<int D, typename Derived>
Line<D> operator-(const Eigen::MatrixBase<Derived> &T, const Line<D> &l) {
    //make sure it's a vector
    static_assert(Derived::RowsAtCompileTime == D);
    static_assert(Derived::ColsAtCompileTime == 1);
    return Line<D>(T - l.start_end().colwise());
}
template<int D, typename Derived>
Line<D> operator-(const Line<D> &l, const Eigen::MatrixBase<Derived> &T) {
    //make sure it's a vector
    static_assert(Derived::RowsAtCompileTime == D);
    static_assert(Derived::ColsAtCompileTime == 1);
    return Line<D>(l.start_end().colwise() - T);
}
template<int D, typename Derived>
Line<D> operator*(const Eigen::EigenBase<Derived> &T, const Line<D> &l) {
    return Line<D>(T * l.start_end());
}
template<int D>
Line<D> operator*(double s, const Line<D> &l) {
    return Line<D>(s * l.start_end());
}
template<int D>
Line<D> operator*(const Line<D> &l, double s) {
    return Line<D>(s * l.start_end());
}
template<int D>
Line<D> operator/(const Line<D> &l, double s) {
    return Line<D>(l.start_end() / s);
}

template<int D>
auto Line<D>::interp(double t) const -> Vec {
    return start() + t * direction();
}

template<int D>
template<typename Derived>
auto Line<D>::project(const Eigen::MatrixBase<Derived> &p) const -> Vec {
    static_assert(Derived::RowsAtCompileTime == D);
    static_assert(Derived::ColsAtCompileTime == 1);

    auto pa = p - start();
    Vec d = direction();
    double len = d.norm();
    d.normalize();
    Vec nd = std::clamp<double>(pa.dot(d), 0, len) * d;
    return nd + start();
}
template<int D>
template<typename Derived>
double Line<D>::distance(const Eigen::MatrixBase<Derived> &p) const {
    return (project(p) - p).norm();
}
template<int D>
template<typename Derived>
Line<D> Line<D>::cwiseProduct(const Eigen::MatrixBase<Derived> &p) const {

    return Vec2(p.asDiagonal() * start_end());
}
template<int D>
template<typename Derived>
Line<D> Line<D>::cwiseQuotient(const Eigen::MatrixBase<Derived> &p) const {
    Vec pinv = 1.0 / p.array();
    return Vec2(pinv.asDiagonal() * start_end());
}
}// namespace mandoline
