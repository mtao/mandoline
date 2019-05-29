#pragma once
#include <mtao/types.h>

namespace mandoline {
    template <int D>
        class Line {
            public:
                using Vec = mtao::Vector<double,D>;
                using Vec2 = mtao::Matrix<double,D,2>;
                Line(const Vec& a, const Vec& b);
                Line(const Vec2& se): m_start_end(se) {}

                auto&& start_end() { return m_start_end; }
                auto&& start_end() const { return m_start_end; }
                auto start() { return start_end().col(0); }
                auto end() { return start_end().col(1); }
                auto start() const { return start_end().col(0); }
                auto end() const { return start_end().col(1); }
                auto direction() const { return end() - start(); }
                double length() const { return direction().norm(); }
                template <typename Derived>
                    double distance(const Eigen::MatrixBase<Derived>& p) const;
                template <typename Derived>
                    Vec project(const Eigen::MatrixBase<Derived>& p) const;
                template <typename Derived>
                    Line cwiseProduct(const Eigen::MatrixBase<Derived>& p) const;
                template <typename Derived>
                    Line cwiseQuotient(const Eigen::MatrixBase<Derived>& p) const;

                Vec interp(double t) const;
            private:
                Vec2 m_start_end;
            public:
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };


    template <int D, typename Derived>
        Line<D> operator+(const Eigen::MatrixBase<Derived>& T, const Line<D>& l);
    template <int D, typename Derived>
        Line<D> operator+(const Line<D>& l, const Eigen::MatrixBase<Derived>& T);

    template <int D, typename Derived>
        Line<D> operator-(const Eigen::MatrixBase<Derived>& T, const Line<D>& l);
    template <int D, typename Derived>
        Line<D> operator-(const Line<D>& l, const Eigen::MatrixBase<Derived>& T);
    template <int D, typename Derived>
        Line<D> operator*(const Eigen::EigenBase<Derived>& T, const Line<D>& l);

    template <int D, typename Derived>
        Line<D> operator*(double s, const Line<D>& l);
    template <int D, typename Derived>
        Line<D> operator*(const Line<D>& l, double s);
    template <int D, typename Derived>
        Line<D> operator/(const Line<D>& l, double s);

}

#include "mandoline/line_impl.hpp"
