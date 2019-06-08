#include "mandoline/cutface.hpp"
#include <mtao/geometry/mesh/triangle_fan.hpp>
#include <mtao/geometry/mesh/triangle/triangle_wrapper.h>
#include <mtao/geometry/mesh/earclipping.hpp>

namespace mandoline {
    template <>
        void CutFace<3>::cache_triangulation(const std::array<mtao::ColVecs2d,3>& V) {
            cache_triangulation(triangulate(V));
        }
    template <>
        void CutFace<3>::cache_triangulation(const mtao::ColVecs3i& F) {
            triangulation = F; 
        }
    template <>
        mtao::ColVecs3i CutFace<3>::triangulate(const std::array<mtao::ColVecs2d,3>& V) const {

            if(is_mesh_face()) {
                return triangulate_fan();
            } else {
                int id = as_axial_id()[0];
                if(indices.size() == 1) {
                    return triangulate_earclipping(V[id]);
                } else {
                    return triangulate_triangle(V[id]);
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
        mtao::ColVecs3i CutFace<3>::triangulate_triangle(const mtao::ColVecs2d& V) const {

            int size = 0;
            for(auto&& v: indices) {
                size += v.size();
            }
            std::set<int> used_vertices;
            for(auto&& v: indices) {
                for(auto&& c: v) {
                    used_vertices.insert(c);
                }
            }
            if(used_vertices.size() == 0) {
                return {};
            }
            mtao::map<int,int> reindexer;
            mtao::map<int,int> unreindexer;
            mtao::ColVecs2d newV(2,used_vertices.size());
            for(auto&& [i,v]: mtao::iterator::enumerate(used_vertices)) {
                reindexer[v] = i;
                newV.col(i) = V.col(v);
                unreindexer[i] = v;
            }
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

            m.fill_attributes();
            for(int i = 0; i < m.EA.cols(); ++i) {m.EA(i) = i;}
            for(int i = 0; i < m.VA.cols(); ++i) {m.VA(i) = i;}
            bool points_added = false;
            static const std::string str ="pcePzQYY";
            auto F = mtao::geometry::mesh::triangle::triangle_wrapper(m,std::string_view(str)).F;
            for(int i = 0; i < F.size(); ++i) {
                auto&& f = F(i);
                if(unreindexer.find(f) == unreindexer.end()) {
                    points_added = true;
                    break;
                } else {
                    f = unreindexer[f];
                }
            }
            if(!points_added) {
                std::set<int> interior;
                for(int i = 0; i < F.cols(); ++i) {
                    mtao::Vec2d B = mtao::Vec2d::Zero();
                    auto f = F.col(i);
                    for(int  j = 0; j < 3; ++j) {
                        B+=V.col(f(j));
                    }
                    B /= 3;
                    double wn = 0;
                    for(auto&& c: indices) {
                        double mywn = mtao::geometry::winding_number(V,c,B);
                        wn += mywn;
                    }
                    //std::cout << wn << " ";
                    if(std::abs(wn) > 1) {
                        interior.insert(i);
                    }

                }
                mtao::ColVecs3i FF(3,interior.size());
                for(auto&& [i,b]: mtao::iterator::enumerate(interior)) {
                    FF.col(i) = F.col(b);
                }
                //std::cout << std::endl;
                return FF;
            } else {
                std::cout << "points were added!" << std::endl;
                return {};
            }
        }

}
