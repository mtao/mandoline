#include "mandoline/operators/volume2.hpp"
#include "mandoline/operators/interpolation2.hpp"
#include "mandoline/operators/volume_impl.hpp"

#include <iostream>

namespace mandoline::operators {

mtao::VecXd face_volumes(const CutCellMesh<2> &ccm, bool from_triangulation) {

    mtao::VecXd FV(ccm.cut_faces().size());
    if (from_triangulation) {
        for (auto &&[i, face] : mtao::iterator::enumerate(ccm.cut_faces())) {
            auto V = ccm.vertices();
            if (face.triangulation) {
                auto &&T = *face.triangulation;
                FV(i) = mtao::geometry::volumes(V, T).sum();
            }
        }
        FV = mtao::eigen::vstack(FV, face_volumes(ccm.exterior_grid));
    } else {
        FV = face_volumes(ccm.exterior_grid);
    }

    if (FV.minCoeff() < 0) {
        mtao::logging::warn() << "Negative face area! warn mtao because this shouldn't happen";
        FV = FV.cwiseAbs();
    }
    return FV;
}


mtao::VecXd edge_lengths(const CutCellMesh<2> &ccm) {

    mtao::VecXd EV;
    //Use barycentric for tri-mesh cut-faces, planar areas for axial cut-faces, and let the adaptive grid do it's thing
    auto mesh_edge_vols = mtao::geometry::volumes(ccm.origV(), ccm.origE());

    if (mesh_edge_vols.size() > 0) {
        // handles all of the mesh cut-edges, gotta do the others differently
        EV = edge_barycentric_volume_matrix(ccm) * mesh_edge_vols;
        EV.resize(ccm.num_edges());
        for (auto &&[i, e] : mtao::iterator::enumerate(ccm.cut_edges())) {
            if (e.is_axial_edge()) {
                auto [dim, coord] = e.as_axial_id();
                dim = 1 - dim;// reverse the dim to be the unbound axis
                auto a = ccm.masked_vertex(e.indices[0]);
                auto b = ccm.masked_vertex(e.indices[1]);
                int coff = b.coord[dim] - a.coord[dim];

                EV(i) = ccm.dx()(dim) * ((b.quot - a.quot)(dim) + coff);
            }
            if (!std::isfinite(EV(i))) {
                EV(i) = 0;
            }
        }
    } else {
        EV.resize(ccm.exterior_grid.boundary_facet_size());
    }
    EV.tail(ccm.exterior_grid.boundary_facet_size()) = ccm.exterior_grid.boundary_facet_volumes();
    return EV;
}

mtao::VecXd dual_edge_lengths(const CutCellMesh<2> &ccm) {

    mtao::VecXd DL(ccm.cut_edges().size());
    DL.setZero();
    if (ccm.cut_edges().size() > 0) {
        for (auto &&[i, edge] : mtao::iterator::enumerate(ccm.cut_edges())) {
            if (edge.is_axial_edge()) {
                int axis = edge.bound_axis();
                DL(i) = ccm.dx()(axis);
            }
        }
        return mtao::eigen::vstack(DL, boundary_facet_volumes(ccm.exterior_grid));
    } else {
        return face_volumes(ccm.exterior_grid);
    }
}
mtao::VecXd dual_hodge1(const CutCellMesh<2> &ccm) {
    auto PV = edge_lengths(ccm);
    auto DV = dual_edge_lengths(ccm);
    mtao::VecXd CV = (DV.array() > 1e-5).select(PV.cwiseQuotient(DV), 0);
    for (int i = 0; i < CV.size(); ++i) {
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
mtao::VecXd primal_hodge1(const CutCellMesh<2> &ccm) {
    auto PV = edge_lengths(ccm);
    auto DV = dual_edge_lengths(ccm);
    mtao::VecXd CV = (PV.array() > 1e-5).select(DV.cwiseQuotient(PV), 0);
    for (int i = 0; i < CV.size(); ++i) {
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
mtao::VecXd dual_hodge2(const CutCellMesh<2> &ccm) {
    auto CV = face_volumes(ccm);
    for (int i = 0; i < CV.size(); ++i) {
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
mtao::VecXd primal_hodge2(const CutCellMesh<2> &ccm) {
    auto CV = face_volumes(ccm);
    for (int i = 0; i < CV.size(); ++i) {
        CV(i) = (std::abs(CV(i)) < 1e-5) ? 0 : (1. / CV(i));
        if (!std::isfinite(CV(i))) {
            CV(i) = 0;
        }
    }
    return CV;
}
}// namespace mandoline::operators
