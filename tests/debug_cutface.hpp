#pragma once
#include <mandoline/cutface.hpp>
#include <iostream>

inline bool check_convex_face_normal(const mtao::ColVecs3d& V, const mandoline::CutFace<3>& F) {
    const auto& N = F.N;
    const auto& inds = *F.indices.begin();
    bool good = true;
    //std::cout << "[[[" << N.transpose() << "]]]";
    //for(int i = 0; i < inds.size(); ++i) {
    //    std::cout << "(" << inds[i] << "): " << V.col(inds[i]).transpose() << "====";
    //}
    //std::cout << std::endl;
    for(int i = 0; i < inds.size(); ++i) {
        int j = (i+1)%inds.size();
        auto o = V.col(inds[i]);
        auto od = V.col(inds[j]);
        mtao::Vec3d d = od - o;
        //std::cout << "(" << d.dot(N) << ")";
        mtao::Vec3d pn = N.cross(d);
        for(int k = 2; k < inds.size(); ++k) {
            double val = (V.col(inds[(i+k)%inds.size()]) - o).dot(pn);
            //std::cout << val << " ";
            good &= val >= 0;
        }
        //std::cout << " || ";
    }
    //std::cout << std::endl;
    return good;
}

