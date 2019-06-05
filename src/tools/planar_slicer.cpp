#include "mandoline/tools/planar_slicer.hpp"
#include "mandoline/construction/cutdata.hpp"
#include "mandoline/cutface.hpp"
#include <mtao/geometry/normal_basis.h>
#include <mtao/geometry/bounding_box.hpp>
#include <mtao/eigen/stack.h>
#include <mtao/geometry/mesh/compactify.hpp>

using namespace mandoline::construction;
using namespace mandoline;
namespace mandoline::tools {

    SliceGenerator::SliceGenerator(const mtao::Vec3d& origin, const mtao::Vec3d& direction) {
        transform.linear().setIdentity();
        transform.linear() = mtao::geometry:: normal_basis(direction);
        transform.translation() = -transform.linear().inverse() * origin;

    }
    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> SliceGenerator::slice(const mtao::ColVecs3d& V, const mtao::ColVecs3i& F) const {

        mtao::ColVecs3d VV = transform * V;
        auto bbox = mtao::geometry::bounding_box(VV);
        double zmax = std::max(bbox.max()(2),-bbox.min()(2));
        bbox.min()(2) = -zmax;
        bbox.max()(2) = zmax;

        mtao::Vec3d shape = bbox.sizes() ;
        shape = (shape.array() > 1e-10).select(shape,1);


        std::vector<Vertex<3>> vertices(V.cols());
        mtao::geometry::grid::StaggeredGrid<double,3> sg = mtao::geometry::grid::StaggeredGrid<double,3>::from_bbox(bbox,std::array<int,3>{{1,1,2}});
        for(int i = 0; i < V.cols(); ++i) {
            mtao::Vec3d v = VV.col(i);
            vertices[i] = Vertex<3>::from_vertex(sg.vertex_grid().local_coord(v));
        }

        CutData data(sg.vertex_grid(),vertices,F);
        data.bake(sg.vertex_grid());
        auto& faces = data.cut_faces();
        auto CV = data.cut_vertices();
        for(int i = 0; i < CV.cols(); ++i) {
            CV.col(i) = sg.vertex_grid().world_coord(CV.col(i));
        }
        auto newV = mtao::eigen::hstack(sg.vertices(),CV);


        std::set<int> plane_verts;
        for(auto&& c: data.crossings()) {
            auto m = c.mask();
            if(m[2] && *m[2] == 1) {
                plane_verts.insert(c.index);
            }
        }
        auto edges = data.edges();
        auto is_boundary = [&](int idx) -> bool {
            if(idx < sg.vertex_size()) {
                return sg.unindex(idx)[2] == 1;
            } else {
                return plane_verts.find(idx) != plane_verts.end();
            }
        };
        std::vector<mtao::ColVecs3i> FF;
        {
            std::set<std::array<int,2>> E;
            for(int i = 0; i < edges.cols(); ++i) {
                auto e = edges.col(i);
                int a = e(0);
                int b = e(1);
                if(is_boundary(a) && is_boundary(b)) {
                    E.insert(std::array<int,2>{{a,b}});
                }
            }




            auto VV = newV.topRows<2>();
            auto ehem = mtao::geometry::mesh::EmbeddedHalfEdgeMesh<double,2>::from_edges(VV,mtao::eigen::stl2eigen(E));
            ehem.make_topology();
            ehem.tie_nonsimple_cells();
            auto A = ehem.signed_areas();
            auto mesh_cells = ehem.cells_multi_component_map();
            for(auto&& [cid,inds]: mesh_cells) {
                if(A.find(cid) != A.end() && A.at(cid) > 0) {
                    CutFace<3> cf(coord_mask<3>(2,1),inds,{{2,1}});
                    if(inds.size() == 1) {
                        FF.push_back(cf.triangulate_earclipping(VV));
                    } else {
                        FF.push_back(cf.triangulate_triangle(VV));
                    }
                }
            }
        }





        double mid = (bbox.min()+bbox.max()/2)(2);
        for(auto&& F: faces) {
            bool above = true;
            for(auto&& i: F.indices) {
                if(newV.col(i)(2) < 0) {
                    above = false;
                    break;
                }
            }
            if(above) {
                FF.push_back(F.triangulate());
            }
        }
        if(FF.size() == 0) {
            return {};
        }
        auto T = transform.inverse();
        for(int i = 0; i < newV.cols(); ++i) {
            newV.col(i) = T * newV.col(i);
        }
        auto newF = mtao::eigen::hstack_iter(FF.begin(),FF.end());
        return  mtao::geometry::mesh::compactify(newV,newF);

    }
}
