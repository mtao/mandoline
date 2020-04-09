#pragma once
#include "mandoline/cutface.hpp"

namespace mandoline {

template<int D>
bool CutFaceBase<D>::is_mesh_face() const {
    return std::holds_alternative<int>(id);
}
template<int D>
bool CutFaceBase<D>::is_axial_face() const {
    return std::holds_alternative<std::array<int, 2>>(id);
}
template<int D>
int CutFaceBase<D>::as_face_id() const {
    return std::get<int>(id);
}
template<int D>
const std::array<int, 2> &CutFaceBase<D>::as_axial_id() const {
    return std::get<std::array<int, 2>>(id);
}
template<int D>
int CutFaceBase<D>::as_axial_axis() const {
    return std::get<std::array<int, 2>>(id)[0];
}
template<int D>
int CutFaceBase<D>::as_axial_coord() const {
    return std::get<std::array<int, 2>>(id)[1];
}
template<int D>
CutFaceBase<D>::operator std::string() const {
    std::stringstream ss;
    ss << "{";

    std::visit([&](auto &&f) {
        using T = std::decay_t<decltype(f)>;
        if constexpr (std::is_same_v<T, int>) {
            ss << "(" << f << ")";
        } else {
            auto [a, b] = f;
            ss << "(" << a << ":" << b << ")";
        }
    },
               id);
    for (auto &&f : indices) {
        ss << "[";
        std::copy(f.begin(), f.end(), std::ostream_iterator<int>(ss, ","));
        ss << "],";
    }

    if (external_boundary) {
        auto [b, s] = *external_boundary;
        if (s) {
            std::cout << "{+" << b << "}";
        } else {
            std::cout << "{-" << b << "}";
        }
    }
    ss << "}";

    return ss.str();
}
template<int D>
template<typename Derived>
void CutFaceBase<D>::update_mask(const std::vector<Vertex<D>> &V, const mtao::geometry::grid::indexing::IndexerBase<D, Derived> &indexer) {
    int offset = indexer.size();
    auto get_mask = [&](int idx) -> coord_mask<D> {
        if (idx >= offset) {
            return V[idx - offset].mask();
        } else {
            return indexer.unindex(idx);
        }
    };
    coord_mask<D> &me = *this;
    me = get_mask(indices.begin()->front());
    for (auto &&ind : indices) {
        for (auto &&i : ind) {
            auto mask = get_mask(i);
            me &= mask;
        }
    }
}

template<typename Derived>
mtao::Vec2d CutFace<2>::brep_centroid(const Eigen::MatrixBase<Derived> &V) const {
    mtao::Vec2d ret = mtao::Vec2d::Zero();
    double tot_vol = 0;
    for (auto &&i : indices) {
        mtao::Vec2d cent = mtao::geometry::curve_centroid(V,i);
        double vol = mtao::geometry::curve_volume(V,i);
        ret += cent * vol;
        tot_vol += vol;
    }

    ret /= tot_vol;
    return ret;
}

template<typename Derived>
mtao::Vec3d CutFace<3>::brep_centroid(const Eigen::MatrixBase<Derived> &V, bool use_triangulation) const {
    mtao::Vec3d ret = mtao::Vec3d::Zero();
    double tot_vol = 0;
    int count = 0;
    for (auto &&i : indices) {
        count += i.size();
        for (auto &&i : i) {
            ret += V.col(i);
        }
        /*
                   auto vol = brep_volume(V,use_triangulation);
                   tot_vol += vol;
                   mtao::Vec3d cent = mtao::geometry::curve_centroid(V,i);
                   ret += vol * cent;
                   */
    }
    /*
               if(tot_vol > 1e-10) {
               ret /= tot_vol;
               }
               */
    ret /= count;
    return ret;
}


template<typename Derived>
double CutFace<3>::brep_volume(const Eigen::MatrixBase<Derived> &V, bool use_triangulation) const {

    if (use_triangulation) {
        return mtao::geometry::brep_volume(V, *triangulation);
    }
    double vol = 0;
    auto loop_volume = [&](const std::vector<int> &loop) {
        double vol = 0;
        mtao::Matrix<double, 3, 4> S;
        S.col(3).setZero();
        S.col(0) = V.col(loop[0]);
        S.col(2) = V.col(loop[1]);

        for (auto it = loop.begin() + 2; it != loop.end(); ++it) {
            S.col(1) = S.col(2);
            S.col(2) = V.col(*it);
            vol += mtao::geometry::volume_signed(S);
        }
        return vol;
        //return mtao::geometry::volume_signed(S);
    };
    for (auto &&loop : indices) {
        vol += loop_volume(loop);
    }
    return vol;
}

template<typename Derived, typename PDerived>
double CutFace<2>::winding_number(const Eigen::MatrixBase<Derived> &V, const Eigen::MatrixBase<PDerived>& p) const {

    double wn = 0;
        for (auto &&f : indices) {
            wn += mtao::geometry::winding_number(V, f, p);
        }
        return wn;
}
template<typename Derived, typename PDerived>
bool CutFace<2>::is_inside(const Eigen::MatrixBase<Derived> &V, const Eigen::MatrixBase<PDerived>& p) const {

    return winding_number(V,p)  > .5;
}
}// namespace mandoline
