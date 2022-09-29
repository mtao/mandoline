#include "mandoline/operators/volume3.hpp"

#include "mandoline/operators/interpolation3.hpp"

namespace mandoline::operators {

balsa::eigen::VecXd cell_volumes(const CutCellMesh<3> &ccm) {
    balsa::eigen::VecXd V(ccm.cells().size());
    V.setZero();
    auto Vs = ccm.vertices();
    balsa::eigen::VecXd face_brep_vols(ccm.faces().size());
    for (auto &&[i, f] : mtao::iterator::enumerate(ccm.faces())) {
        face_brep_vols(i) = f.brep_volume(Vs);
    }

    for (auto &&[k, c] : mtao::iterator::enumerate(ccm.cells())) {
        V(k) = c.volume(face_brep_vols, ccm.folded_faces());
        // V(k) = c.volume(Vs,faces);
        // V(k) = mtao::geometry::brep_volume(Vs,c.triangulated(faces));
        // std::cout << (V(k) / std::abs(c.volume(Vs,faces))) << std::endl;
    }

#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
    return balsa::eigen::vstack(V, ccm.exterior_grid().cell_volumes());
#else
    return balsa::eigen::vstack(
        V, balsa::eigen::VecXd::Constant(ccm.dx().prod(),
                                 ccm.exterior_grid().num_cells()));
#endif
}
balsa::eigen::VecXd face_volumes(const CutCellMesh<3> &ccm, bool from_triangulation) {
    balsa::eigen::VecXd FV(ccm.faces().size());
    if (from_triangulation) {
        for (auto &&[i, face] : mtao::iterator::enumerate(ccm.faces())) {
            auto V = ccm.vertices();
            if (face.triangulation) {
                auto &&T = *face.triangulation;
                FV(i) = mtao::geometry::volumes(V, T).sum();
            }
        }
        FV = balsa::eigen::vstack(FV, ccm.exterior_grid().face_volumes());
    } else {
        // Use barycentric for tri-mesh cut-faces, planar areas for axial
        // cut-faces, and let the adaptive grid do it's thing
        auto trimesh_vols = mtao::geometry::volumes(ccm.origV(), ccm.origF());

        if (trimesh_vols.size() > 0) {
            FV = face_barycentric_volume_matrix(ccm) * trimesh_vols;
            auto subVs = ccm.compute_subVs();
            for (auto &&[i, f] : mtao::iterator::enumerate(ccm.faces())) {
                if (f.is_axial_face()) {
                    auto [dim, coord] = f.as_axial_id();
                    auto &vol = FV(i) = 0;
                    auto &V = subVs[dim];
                    for (auto &&indices : f.indices) {
                        vol += mtao::geometry::curve_volume(V, indices);
                    }
                }
                if (!std::isfinite(FV(i))) {
                    FV(i) = 0;
                }
            }
        } else {
            FV.resize(ccm.exterior_grid().num_faces());
        }

#if defined(MANDOLINE_USE_ADAPTIVE_GRID)
        FV.tail(ccm.exterior_grid().num_faces()) =
            ccm.exterior_grid().face_volumes();
#else
        auto FVT = FV.tail(ccm.exterior_grid().num_faces());
        balsa::eigen::Vec3d ddx = ccm.dx().prod() / ccm.dx();
        for (auto &&[f, axis] : mtao::iterator::zip(
                 FVT, ccm.exterior_grid().boundary_facet_axes())) {
            f = ddx(axis);
        }
#endif
    }

    if (FV.minCoeff() < 0) {
        mtao::logging::warn()
            << "Negative face area! warn mtao because this shouldn't happen";
        FV = FV.cwiseAbs();
    }
    return FV;
}

balsa::eigen::VecXd dual_edge_lengths(const CutCellMesh<3> &ccm) {
    balsa::eigen::VecXd DL = balsa::eigen::VecXd::Zero(ccm.face_size());

    auto &dx = ccm.Base::dx();
    auto g = ccm.exterior_grid().cell_ownership_grid();
    for (auto &&c : ccm.cells()) {
        auto &gc = c.grid_cell;
        for (auto &&[fidx, s] : c) {
            auto &f = ccm.faces()[fidx];
            if (f.external_boundary) {
                auto [cid, s] = *f.external_boundary;
                if (cid >= 0) {
                    auto &&c = ccm.exterior_grid().cell(g.get(cid));
                    balsa::eigen::Vec3d C = c.center();
                    C.array() -= .5;
                    C -= balsa::eigen::stl2eigen(gc).cast<double>();

                    DL(fidx) = (dx.asDiagonal() * C).norm();
                }
            }
        }
    }
    for (auto &&[fidx, f] : mtao::iterator::enumerate(ccm.faces())) {
        if (f.is_axial_face() && !f.external_boundary) {
            int ba = f.as_axial_axis();
            DL(fidx) = dx(ba);
        }
    }
    auto ADL = ccm.exterior_grid().dual_edge_lengths();
    DL.tail(ADL.rows()) = ADL;
    return DL;
}
balsa::eigen::VecXd dual_hodge2(const CutCellMesh<3> &ccm) {
    auto PV = face_volumes(ccm);
    auto DV = dual_edge_lengths(ccm);
    balsa::eigen::VecXd CV = (DV.array() > 1e-5).select(PV.cwiseQuotient(DV), 0);
    for (int i = 0; i < CV.size(); ++i) {
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
balsa::eigen::VecXd primal_hodge2(const CutCellMesh<3> &ccm) {
    auto PV = face_volumes(ccm);
    auto DV = dual_edge_lengths(ccm);
    balsa::eigen::VecXd CV = (PV.array() > 1e-5).select(DV.cwiseQuotient(PV), 0);
    for (int i = 0; i < CV.size(); ++i) {
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
balsa::eigen::VecXd dual_hodge3(const CutCellMesh<3> &ccm) {
    auto CV = cell_volumes(ccm);
    for (int i = 0; i < CV.size(); ++i) {
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
balsa::eigen::VecXd primal_hodge3(const CutCellMesh<3> &ccm) {
    auto CV = cell_volumes(ccm);
    for (int i = 0; i < CV.size(); ++i) {
        CV(i) = (std::abs(CV(i)) < 1e-5) ? 0 : (1. / CV(i));
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
}  // namespace mandoline::operators
