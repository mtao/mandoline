#include "mandoline/cutface.hpp"
#include <mtao/geometry/mesh/triangle_fan.hpp>
#include <mtao/geometry/mesh/triangle/triangle_wrapper.h>
#include <mtao/geometry/mesh/earclipping.hpp>

namespace mandoline {
    template <>
        void CutFace<3>::cache_triangulation(const std::array<mtao::ColVecs2d,3>& V, bool add_verts) {
            auto [VV,F] = triangulate(V,add_verts);
            cache_triangulation(VV,F);
        }
    template <>
        void CutFace<3>::cache_triangulation(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F) {
            triangulated_vertices = V;
            cache_triangulation(F);
        }
    template <>
        void CutFace<3>::cache_triangulation(const mtao::ColVecs3i& F) {
            triangulation = F; 
        }
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate(const std::array<mtao::ColVecs2d,3>& V) const {
            return std::get<1>(triangulate(V,false));
        }
    template <>
        std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutFace<3>::triangulate(const std::array<mtao::ColVecs2d,3>& V, bool add_vertices) const {
            if(is_mesh_face()) {
                return {{},triangulate_fan()};
            } else {
                int id = as_axial_id()[0];
                if(indices.size() == 1) {
                    return {{},triangulate_earclipping(V[id])};
                } else {
                    return triangulate_triangle(V[id],add_vertices);
                }
            }

        }
    template <>
        std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutFace<3>::triangulate(const mtao::ColVecs2d& V, bool add_vertices) const {

            if(is_mesh_face()) {
                return {{},triangulate_fan()};
            } else {
                if(indices.size() == 1) {
                    return {{},triangulate_earclipping(V)};
                } else {
                    return triangulate_triangle(V,add_vertices);
                }
            }
        }
    template <>
        mtao::ColVecs3i CutMeshFace<3>::triangulate() const {
            return mtao::geometry::mesh::triangle_fan(indices);
        }

    template <>
        mtao::ColVecs3i CutFace<3>::triangulate_fan() const {
            assert(indices.size() == 1);
            return mtao::geometry::mesh::triangle_fan(*indices.begin());
        }
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate_earclipping(const mtao::ColVecs2d& V) const {
            assert(indices.size() == 1);
            return mtao::geometry::mesh::earclipping(V,indices);
        }

    template <>
        std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutFace<3>::triangulate_triangle(const mtao::ColVecs2d& V, bool do_add_vertices) const {

            //collect the number of edges
            int size = 0;
            for(auto&& v: indices) {
                size += v.size();
            }


            //check the total number of vertices
            std::set<int> used_vertices;
            for(auto&& v: indices) {
                for(auto&& c: v) {
                    used_vertices.insert(c);
                }
            }
            if(used_vertices.size() == 0) {
                return {};
            }

            //compactify the vertices
            mtao::map<int,int> reindexer;
            mtao::map<int,int> unreindexer;
            mtao::ColVecs2d newV(2,used_vertices.size());
            for(auto&& [i,v]: mtao::iterator::enumerate(used_vertices)) {
                reindexer[v] = i;
                newV.col(i) = V.col(v);
                unreindexer[i] = v;
            }

            //compactify the edges
            mtao::ColVecs2i E(2,size);
            mtao::vector<mtao::Vec2d> holes;
            size = 0;
            for(auto&& curve: indices) {
                for(int i = 0; i < curve.size(); ++i) {
                    auto e = E.col(size++);
                    e(0) = reindexer[curve[i]];
                    e(1) = reindexer[curve[(i+1)%curve.size()]];
                }
            }
            mtao::geometry::mesh::triangle::Mesh m(newV,E);

            //set some arbitrary attributes?
            m.fill_attributes();
            for(int i = 0; i < m.EA.cols(); ++i) {m.EA(i) = i;}
            for(int i = 0; i < m.VA.cols(); ++i) {m.VA(i) = i;}
            bool points_added = false;
            mtao::ColVecs2d newV2;
            mtao::ColVecs3i newF;

            if(do_add_vertices) {
                //static const std::string str ="zPa.01qepcDQ";
                static const std::string str ="pcePzQYY";
                std::cerr << "I bet im about to crash" << std::endl;
                auto nm = mtao::geometry::mesh::triangle::triangle_wrapper(m,std::string_view(str));
                std::cerr << "I lost a bet" << std::endl;
                if(nm.V.cols() > newV.cols()) {
                    newV2 = nm.V;
                    points_added = true;
                }

                newF = nm.F;
            } else {
                static const std::string str ="pcePzQYY";
                auto nm = mtao::geometry::mesh::triangle::triangle_wrapper(m,std::string_view(str));
                newF = nm.F;
            }
            for(int i = 0; i < newF.size(); ++i) {
                auto&& f = newF(i);
                if(unreindexer.find(f) == unreindexer.end()) {
                    points_added = true;
                    break;
                }
            }
            if(!points_added) {
                for(int i = 0; i < newF.size(); ++i) {
                    auto&& f = newF(i);
                    f = unreindexer[f];
                }
            }


            std::set<int> interior;
            for(int i = 0; i < newF.cols(); ++i) {
                mtao::Vec2d B = mtao::Vec2d::Zero();
                auto f = newF.col(i);
                for(int  j = 0; j < 3; ++j) {
                    if(points_added) {
                        B+=newV2.col(f(j));
                    } else {
                        B+=V.col(f(j));
                    }
                }
                B /= 3;
                double wn = 0;
                for(auto&& c: indices) {
                    double mywn = mtao::geometry::winding_number(V,c,B);
                    wn += mywn;
                }
                //std::cout << wn << " ";
                if(std::abs(wn) > .5) {
                    interior.insert(i);
                }

            }
            mtao::ColVecs3i FF(3,interior.size());
            for(auto&& [i,b]: mtao::iterator::enumerate(interior)) {
                FF.col(i) = newF.col(b);
            }
            //std::cout << std::endl;
            if(points_added && !do_add_vertices) {
                std::cout << "points were added when they shouldnt have!" << std::endl;
                return {};
            }  else {
                mtao::ColVecs3d newV3;
                newV3.resize(3,newV2.cols());
                int axis = as_axial_axis();
                int coord = as_axial_coord();
                newV3.row((axis+1)%3) = newV2.row(0);
                newV3.row((axis+2)%3) = newV2.row(1);
                newV3.row(axis).setConstant(double(coord));

                std::cout << "Triangulation with vertex added: " << newV3.cols() << " / " << FF.maxCoeff() <<std::endl;
                return {newV3,FF};
            }
        }

}
