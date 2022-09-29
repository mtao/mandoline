#include "mandoline/operators/centroids3.hpp"

#include <balsa/eigen/shape_checks.hpp>

namespace mandoline::operators {
namespace {
template <balsa::eigen::concepts::ColVecs2Compatible ColVecs>
balsa::eigen::Vec3d brep_centroid(const CutFace<3> &cutface, const ColVecs &V,
                          bool use_triangulation = false) {
    balsa::eigen::Vec3d ret = balsa::eigen::Vec3d::Zero();
    double tot_vol = 0;
    int count = 0;
    for (auto &&i : cutface.indices) {
        count += i.size();
        for (auto &&i : i) {
            ret += V.col(i);
        }
        /*
                   auto vol = brep_volume(V,use_triangulation);
                   tot_vol += vol;
                   balsa::eigen::Vec3d cent = mtao::geometry::curve_centroid(V,i);
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
}  // namespace

balsa::eigen::ColVecs3d cell_centroids(const CutCellMesh<3> &ccm) {
    /*
           VecX V(cell_size());
           V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();
           */
    balsa::eigen::ColVecs3d face_brep_cents(3, ccm.face_size());
    auto Vs = ccm.vertices();
    for (auto &&[i, f] : mtao::iterator::enumerate(ccm.faces())) {
        face_brep_cents.col(i) = f.brep_centroid(Vs);
        // face_brep_cents.col(i) = f.brep_volume(Vs) * f.brep_centroid(Vs);
    }
    balsa::eigen::ColVecs3d V(3, ccm.cell_size());
    V.setZero();
    auto vols = ccm.cell_volumes();
    for (auto &&[k, c] : mtao::iterator::enumerate(ccm.cells())) {
        V.col(k) = c.moment(face_brep_cents);
        // V.col(k) = c.moment(face_brep_cents) / vols(k);
    }

    ccm.exterior_grid().cell_centroids(V);
    return V;
}
balsa::eigen::ColVecs3d face_centroids(const CutCellMesh<3> &ccm) {
    /*
           VecX V(cell_size());
           V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();
           */
    balsa::eigen::ColVecs3d ret(3, ccm.cut_face_size());
    auto Vs = ccm.vertices();

    auto subVs = ccm.compute_subVs();
    for (auto &&[i, f] : mtao::iterator::enumerate(ccm.faces())) {
        auto v = ret.col(i);
        if (f.is_mesh_face()) {
            v = mtao::geometry::curve_centroid(Vs, *f.indices.begin());
        } else {
            const int axis = f.as_axial_axis();
            const auto &V = subVs[axis];
            v.setZero();
            double vol = 0;
            for (auto &&ind : f.indices) {
                double myvol = mtao::geometry::curve_volume(V, ind);
                v += myvol * mtao::geometry::curve_centroid(Vs, ind);
                vol += myvol;
            }
            v /= vol;
        }
    }
    auto r = ccm.exterior_grid().boundary_centroids();
    return balsa::eigen::hstack(ret, r);
}

}  // namespace mandoline::operators
