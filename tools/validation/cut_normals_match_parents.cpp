#include "mandoline/mesh3.hpp"
#include <mtao/geometry/mesh/face_normals.hpp>


namespace {
    int get_convex_hull_wedge(const mandoline::CutFace<3>& f, const mtao::ColVecs2d& v) {
        if(f.is_axial_face()) {
            for(auto&& ind: f.indices) {
            }
        }
    }
}

bool cut_tangents_match_parents(const mandoline::CutCellMesh<2>& ccm) {
    bool good = true;
    mtao::ColVecs3d V = ccm.vertices();
    mtao::ColVecs3d ON = mtao::geometry::mesh::face_normals(ccm.origV(),ccm.origF());
    /*
    std::array<mtao::ColVecs2d,3> VV;
    for(int i = 0; i < 3; ++i) {
        auto& v = VV[i];
        v.resize(2,V.cols());
        v.row(0) = V.row((i+1)%3);
        v.row(1) = V.row((i+2)%3);
    }
    */


    for(const mandoline::CutFace<3>& f: ccm.cut_faces()) {
        const mtao::Vec3d& N = f.N;
        if(f.is_axial_face()) {
            int axis = f.as_axial_axis();
            if((N - mtao::Vec3d::Unit(axis)).norm() > 1e-5) {
                std::cout << "Bad axial("<< axis<< ") normal: " << N.transpose() << ": " << std::string(f) << std::endl;
                good = false;
            }
            /*
            for(auto&& ind: f.indices) {
                auto& inds = *F.indices.begin();
                // convexity check
                bool is_convex = true;
                bool is_first = true;
                bool sign = false;
                for(int j = 0; j < inds.size(); ++j) {
                    auto x = VV.col(j);
                    for(int k = 1; k < inds.size()-1; ++k) {
                        auto a = VV.col(inds[(j+k)%inds.size()]);
                        auto b = VV.col(inds[(j+k+1)%inds.size()]);

                        mtao::Vec2d c = b - a;
                        mtao::Vec2d d = x - a;
                        double cross = c.x() * d.y() - c.y() * d.x();
                        if(is_first) {
                            sign = cross;
                        } else {
                            is_convex = sign == cross;
                            if(!is_convex) {
                                break;
                            }
                        }
                    }
                }
                if(is_convex) {
                    auto a = V.col(inds[0]);
                    auto b = V.col(inds[1]);
                    auto c = V.col(inds[2]);
                    auto d = b-a;
                    auto e = c-b;
                    auto f = d.cross(e);
                    if(f.cross(N).norm() > 1e-5) {
                        mtao::logging::error() << "Faces should be parallel with the normal" << f.transpose();
                    } else if(f.dot(N) < 0) {
                        mtao::logging::error() << "face facing teh wrong way" << f.transpose() << " : " << idx;
                    }
                }
            }
        }
                */


        } else {// if(f.is_mesh_face())
            int parent_face = f.as_face_id();
            auto PN = ON.col(parent_face);
            if((N -PN).norm() > 1e-5) {
                std::cout << "Bad mesh("<< parent_face<< ") normal: " << N.transpose() << "!=" << PN.transpose() << ": " << std::string(f) << std::endl;
                good = false;
            }

        }
    }
}

