#pragma once

#include <mtao/types.h>

template <typename T>
mtao::ColVectors<T,3> enright_velocities(const mtao::ColVectors<T,3>& P) {
    mtao::ColVectors<T,3> V(P.rows(),P.cols());
    for(int i = 0; i < P.cols(); ++i) {
        auto v = V.col(i);
        auto p = P.col(i);
        const auto& x = p.x();
        const auto& y = p.y();
        const auto& z = p.z();
        T s2x = std::sin(2*M_PI*x);
        T s2y = std::sin(2*M_PI*y);
        T s2z = std::sin(2*M_PI*z);
        v(0) = 2 * s2y * s2z * std::pow<T>(std::sin(M_PI*x),2);
        v(1) = -   s2x * s2z * std::pow<T>(std::sin(M_PI*y),2);
        v(2) = -   s2x * s2y * std::pow<T>(std::sin(M_PI*z),2);
    }
    return V;
}
template <typename T>
mtao::ColVectors<T,3> enright_velocities(const mtao::ColVectors<T,3>& P, T time) {
    if(time < .5) {
        return enright_velocities(P);
    } else {
        return -enright_velocities(P);
    }
}
