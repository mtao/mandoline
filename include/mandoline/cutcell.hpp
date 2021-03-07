#pragma once

#include <mtao/eigen/stack.h>
#include <mtao/geometry/volume.h>

#include <Eigen/Sparse>
#include <mtao/geometry/mesh/compactify.hpp>
#include <vector>

#include "cutmesh.pb.h"
#include "mandoline/cutface.hpp"
namespace mandoline {
struct CutCell : public std::map<int, bool> {
    using coord_type = std::array<int, 3>;
    int index = -1;
    int region = -1;
    coord_type grid_cell;

    operator std::string() const;
    void serialize(protobuf::CutCell &) const;
    static CutCell from_proto(const protobuf::CutCell &);

    std::vector<Eigen::Triplet<double>> boundary_triplets() const;
    mtao::ColVecs3i triangulated(const std::vector<CutFace<3>> &Fs) const;
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i>
    triangulated_with_additional_vertices(const std::vector<CutFace<3>> &Fs,
                                          int vertex_offset) const;

    /*
           coord_type grid_cell(const std::vector<CutFace<3>>& F) const;
           */

    template <typename Derived>
    double volume(const Eigen::MatrixBase<Derived> &V,
                  const mtao::vector<CutFace<3>> &Fs) const;
    double volume(const mtao::VecXd &face_brep_vols) const;
    template <typename Derived>
    mtao::Vec3d centroid(const Eigen::MatrixBase<Derived> &V,
                         const mtao::vector<CutFace<3>> &Fs) const;
    mtao::Vec3d moment(const mtao::ColVecs3d &face_brep_cents) const;
    template <typename Derived>
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> get_mesh(
        const Eigen::MatrixBase<Derived> &V,
        const mtao::vector<CutFace<3>> &Fs) const;

    template <typename Derived, typename VecType>
    bool contains(const Eigen::MatrixBase<Derived> &V,
                  const mtao::vector<CutFace<3>> &Fs,
                  const Eigen::MatrixBase<VecType> &v) const;

    template <typename Derived, typename VecType>
    double solid_angle(const Eigen::MatrixBase<Derived> &V,
                       const mtao::vector<CutFace<3>> &Fs,
                       const Eigen::MatrixBase<VecType> &v) const;

    template <typename Derived>
    double volume(const Eigen::MatrixBase<Derived> &V,
                  const mtao::vector<CutFace<3>> &Fs,
                  const std::set<int> &duplicated_faces) const;
    double volume(const mtao::VecXd &face_brep_vols,
                  const std::set<int> &duplicated_faces) const;
    template <typename Derived>
    mtao::Vec3d centroid(const Eigen::MatrixBase<Derived> &V,
                         const mtao::vector<CutFace<3>> &Fs,
                         const std::set<int> &duplicated_faces) const;
    mtao::Vec3d moment(const mtao::ColVecs3d &face_brep_cents,
                       const std::set<int> &duplicated_faces) const;

    template <typename Derived, typename VecType>
    bool contains(const Eigen::MatrixBase<Derived> &V,
                  const mtao::vector<CutFace<3>> &Fs,
                  const Eigen::MatrixBase<VecType> &v,
                  const std::set<int> &duplicated_faces) const;

    template <typename Derived, typename VecType>
    double solid_angle(const Eigen::MatrixBase<Derived> &V,
                       const mtao::vector<CutFace<3>> &Fs,
                       const Eigen::MatrixBase<VecType> &v,
                       const std::set<int> &duplicated_faces) const;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <typename Derived>
double CutCell::volume(const Eigen::MatrixBase<Derived> &V,
                       const mtao::vector<CutFace<3>> &Fs) const {
    double vol = 0;

    for (auto &&[f, b] : *this) {
        auto &&F = Fs[f];
        double sign = b ? 1 : -1;
        vol += sign * Fs[f].brep_volume(V);
    }
    return vol;
}
template <typename Derived>
mtao::Vec3d CutCell::centroid(const Eigen::MatrixBase<Derived> &V,
                              const mtao::vector<CutFace<3>> &Fs) const {
    double vol = 0;
    mtao::Vec3d mom = mtao::Vec3d::Zero();

    for (auto &&[f, b] : *this) {
        auto &&F = Fs[f];
        double sign = b ? 1 : -1;
        vol += sign * Fs[f].brep_volume(V);
        mom += sign * Fs[f].brep_centroid(V);
    }
    mom /= vol;
    return mom;
}
template <typename Derived>
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutCell::get_mesh(
    const Eigen::MatrixBase<Derived> &V,
    const mtao::vector<CutFace<3>> &Fs) const {
    std::vector<mtao::ColVecs3i> FF;
    for (auto &&[fidx, s] : *this) {
        auto &&f = Fs[fidx];
        assert(bool(f.triangulation));
        FF.push_back(*f.triangulation);
    }
    auto F = mtao::eigen::hstack_iter(FF.begin(), FF.end());
    return mtao::geometry::mesh::compactify(V, F);
}

template <typename Derived, typename VecType>
bool CutCell::contains(const Eigen::MatrixBase<Derived> &V,
                       const mtao::vector<CutFace<3>> &Fs,
                       const Eigen::MatrixBase<VecType> &v) const {
    return solid_angle(V, Fs, v) > .5;
}
template <typename Derived, typename VecType>
double CutCell::solid_angle(const Eigen::MatrixBase<Derived> &V,
                            const mtao::vector<CutFace<3>> &Fs,
                            const Eigen::MatrixBase<VecType> &v) const {
    double sa = 0;
    for (auto &&[fid, sgn] : *this) {
        sa += (sgn ? -1 : 1) * Fs[fid].solid_angle(V, v);
    }
    return sa;
}

template <typename Derived>
double CutCell::volume(const Eigen::MatrixBase<Derived> &V,
                       const mtao::vector<CutFace<3>> &Fs,
                       const std::set<int> &duplicated_faces) const {
    double vol = 0;

    for (auto &&[f, b] : *this) {
        if (duplicated_faces.contains(f)) {
            continue;
        }
        auto &&F = Fs[f];
        double sign = b ? 1 : -1;
        vol += sign * Fs[f].brep_volume(V);
    }
    return vol;
}
template <typename Derived>
mtao::Vec3d CutCell::centroid(const Eigen::MatrixBase<Derived> &V,
                              const mtao::vector<CutFace<3>> &Fs,
                              const std::set<int> &duplicated_faces) const {
    double vol = 0;
    mtao::Vec3d mom = mtao::Vec3d::Zero();

    for (auto &&[f, b] : *this) {
        if (duplicated_faces.contains(f)) {
            continue;
        }
        auto &&F = Fs[f];
        double sign = b ? 1 : -1;
        vol += sign * Fs[f].brep_volume(V);
        mom += sign * Fs[f].brep_centroid(V);
    }
    mom /= vol;
    return mom;
}

template <typename Derived, typename VecType>
bool CutCell::contains(const Eigen::MatrixBase<Derived> &V,
                       const mtao::vector<CutFace<3>> &Fs,
                       const Eigen::MatrixBase<VecType> &v,
                       const std::set<int> &duplicated_faces) const {
    return solid_angle(V, Fs, v, duplicated_faces) > .5;
}
template <typename Derived, typename VecType>
double CutCell::solid_angle(const Eigen::MatrixBase<Derived> &V,
                            const mtao::vector<CutFace<3>> &Fs,
                            const Eigen::MatrixBase<VecType> &v,
                            const std::set<int> &duplicated_faces) const {
    double sa = 0;
    for (auto &&[fid, sgn] : *this) {
        if (duplicated_faces.contains(fid)) {
            continue;
        }
        sa += (sgn ? -1 : 1) * Fs[fid].solid_angle(V, v);
    }
    return sa;
}
}  // namespace mandoline
