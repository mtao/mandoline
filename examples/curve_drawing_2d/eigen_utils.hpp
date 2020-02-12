#pragma once

#include <glm/gtc/type_ptr.hpp>

#include <mtao/types.hpp>



template <typename T> inline auto glm2eigen(const glm::tvec2<T>& v) {return Eigen::Map<const mtao::Vector<T,2>>(glm::value_ptr(v));}
template <typename T> inline auto glm2eigen(const glm::tvec3<T>& v) {return Eigen::Map<const mtao::Vector<T,3>>(glm::value_ptr(v));}
template <typename T> inline auto glm2eigen(const glm::tvec4<T>& v) {return Eigen::Map<const mtao::Vector<T,4>>(glm::value_ptr(v));}

/*
template <typename Object>
auto eigenMap( Object& obj) {
   using T =  
}
*/
#ifdef ARRAY2_H
template <typename T,typename Storage>
inline auto eigenVecMap( const Array2<T,Storage>& obj) {
    return Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1>>(obj.a.data,obj.size());
}
template <typename T, typename Storage>
inline auto eigenVecMap(Array2<T, Storage>& obj) {
    return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>>(obj.a.data,obj.size());
}
#endif
#ifdef ARRAY1_H
template <typename T>
inline auto eigenVecMap( const Array1<T>& obj) {
    return Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1>>(obj.data,obj.size());
}
template <typename T>
inline auto eigenVecMap(Array1<T>& obj) {
    return Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>>(obj.data,obj.size());
}

#endif//ARRAY2_H

template <typename Derived>
auto pinv(const Eigen::MatrixBase<Derived>& M) {
    auto N = M.eval();
    for(int i = 0; i < N.rows() ;++i) {
    for(int j = 0; j < N.cols() ;++j) {
        auto&& n = N(i,j);
        if(std::isfinite(n) && n != 0) {
            n = 1.0/n;
        }
    }
    }
    return N;
}

