#include "mandoline/operators/volume.hpp"
#include "mandoline/operators/interpolation.hpp"


namespace mandoline::operators {

    mtao::VecXd cell_volumes(const CutCellMesh<3>& ccm) {
        mtao::VecXd V(ccm.cells().size());
        V.setZero();
        auto Vs = ccm.vertices();
        mtao::VecXd face_brep_vols(ccm.faces().size());
        for(auto&& [i,f]: mtao::iterator::enumerate(ccm.faces())) {
            face_brep_vols(i) = f.brep_volume(Vs);
        }

        for(auto&& [k,c]: mtao::iterator::enumerate(ccm.cells())) {
            V(k) = c.volume(face_brep_vols);
            //V(k) = c.volume(Vs,faces);
            //V(k) = mtao::geometry::brep_volume(Vs,c.triangulated(faces));
            //std::cout << (V(k) / std::abs(c.volume(Vs,faces))) << std::endl;
        }


        return mtao::eigen::vstack(V,ccm.adaptive_grid().cell_volumes());
    }
    mtao::VecXd face_volumes(const CutCellMesh<3>& ccm, bool from_triangulation) {

        mtao::VecXd FV(ccm.faces().size());
        if(from_triangulation) {
            for(auto&& [i,face]: mtao::iterator::enumerate(ccm.faces())) {
                auto V = ccm.vertices();
                if(face.triangulation) {
                    auto&& T = *face.triangulation;
                    FV(i) = mtao::geometry::volumes(V,T).sum();
                }
            }
            FV = mtao::eigen::vstack(FV,ccm.adaptive_grid().face_volumes());
        } else {
            //Use barycentric for tri-mesh cut-faces, planar areas for axial cut-faces, and let the adaptive grid do it's thing
            auto trimesh_vols = mtao::geometry::volumes(ccm.origV(),ccm.origF());

            if(trimesh_vols.size() > 0) {
                FV = face_barycentric_volume_matrix(ccm) * trimesh_vols;
                auto subVs = ccm.compute_subVs();
                for(auto&& [i,f]: mtao::iterator::enumerate(ccm.faces())) {
                    if(f.is_axial_face()) {
                        auto [dim,coord] = f.as_axial_id();
                        auto& vol = FV(i) = 0;
                        auto& V = subVs[dim];
                        for(auto&& indices: f.indices) {
                            vol += mtao::geometry::curve_volume(V,indices);
                        }

                    }
                    if(!std::isfinite(FV(i))) {
                        FV(i) = 0;
                    }
                }
            } else {
                FV.resize(ccm.adaptive_grid().num_faces());
            }
            FV.tail(ccm.adaptive_grid().num_faces()) = ccm.adaptive_grid().face_volumes();
        }

        if(FV.minCoeff() < 0) {
            mtao::logging::warn() << "Negative face area! warn mtao because this shouldn't happen";
            FV = FV.cwiseAbs();
            
        }
        return FV;
    }



    mtao::VecXd dual_edge_lengths(const CutCellMesh<3>& ccm) {
        mtao::VecXd DL = mtao::VecXd::Zero(ccm.face_size());

        auto& dx = ccm.Base::dx();
        auto g = ccm.adaptive_grid().cell_ownership_grid();
        for(auto&& c: ccm.cells()) {
            auto& gc = c.grid_cell;
            for(auto&& [fidx,s]: c) {
                auto& f = ccm.faces()[fidx];
                if(f.external_boundary) {
                    auto [cid,s] = *f.external_boundary;
                    if(cid >= 0) {
                        auto&& c = ccm.adaptive_grid().cell(g.get(cid));
                        mtao::Vec3d C = c.center();
                        C.array() -= .5;
                        C -= mtao::eigen::stl2eigen(gc).cast<double>();

                        DL(fidx) = (dx.asDiagonal() * C).norm();
                    }
                }
            }
        }
        for(auto&& [fidx,f]: mtao::iterator::enumerate(ccm.faces())) {
            if(f.is_axial_face() && !f.external_boundary) {
                int ba = f.as_axial_axis();
                DL(fidx) = dx(ba);
            }
        }
        auto ADL = ccm.adaptive_grid().dual_edge_lengths();
        DL.tail(ADL.rows()) = ADL;
        return DL;
    }
    mtao::VecXd dual_hodge2(const CutCellMesh<3>& ccm) {
        auto PV = face_volumes(ccm);
        auto DV = dual_edge_lengths(ccm);
        mtao::VecXd  CV = (DV.array() > 1e-5).select(PV.cwiseQuotient(DV),0);
        for(int i = 0; i < CV.size(); ++i) {
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }
    mtao::VecXd primal_hodge2(const CutCellMesh<3>& ccm) {
        auto PV = face_volumes(ccm);
        auto DV = dual_edge_lengths(ccm);
        mtao::VecXd CV = (PV.array() > 1e-5).select(DV.cwiseQuotient(PV),0);
        for(int i = 0; i < CV.size(); ++i) {
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }
    mtao::VecXd dual_hodge3(const CutCellMesh<3>& ccm) {
        auto CV = cell_volumes(ccm);
        for(int i = 0; i < CV.size(); ++i) {
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }
    mtao::VecXd primal_hodge3(const CutCellMesh<3>& ccm) {
        auto CV = cell_volumes(ccm);
        for(int i = 0; i < CV.size(); ++i) {
            CV(i) = (std::abs(CV(i)) < 1e-5) ? 0 : (1. / CV(i));
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }
}
