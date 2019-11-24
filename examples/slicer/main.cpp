#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/bounding_box.hpp"
#include "mtao/geometry/mesh/write_obj.hpp"
#include <mtao/types.hpp>
#include <mtao/cmdline_parser.hpp>
#include "mtao/opengl/Window.h"
#include <iostream>
#include "imgui.h"
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/stack.h>
#include <memory>
#include <algorithm>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <mtao/opengl/drawables.h>
#include <mtao/opengl/objects/mesh.h>
#include <mtao/opengl/objects/bbox.h>
#include <Corrade/Utility/Arguments.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <mandoline/tools/planar_slicer.hpp>
#include <Magnum/Primitives/Cube.h>
#include <mtao/geometry/mesh/vertex_normals.hpp>


#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:
        int current_ccm_cell = 0;
        int current_ccm_face = 0;
        std::optional<std::string> filename;

        mandoline::CutCellMesh<3> ccm;
        mtao::ColVecs3d V;
        mtao::ColVecs4d C;
        mtao::ColVecs3i F;
        mandoline::tools::SliceGenerator sg;



        mtao::Vec3f origin = mtao::Vec3f::Zero(), direction = mtao::Vec3f::Unit(2);



        void update_slice() {
            std::tie(V,F) = sg.slice(origin.cast<double>(),direction.cast<double>());
            if(F.cols() > 0) {

                slice_mesh.setTriangleBuffer(V.cast<float>().eval(), F.cast<unsigned int>());
                Eigen::SparseMatrix<double> BM = sg.barycentric_map();
                mtao::ColVecs4d C2 = C * BM.transpose();
                slice_mesh.setColorBuffer(C2.cast<float>().eval());
            }


            Magnum::Math::Matrix4<float> trans;
            //auto T = mandoline::tools::SliceGenerator::get_transform(origin.cast<double>(),direction.cast<double>());
            for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                //trans[i][j] = T(i,j);
            }
            }
            //slice_object.setTransformation(trans);
        }





        MeshViewer(const Arguments& args): Window3(args) {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            std::tie(V,F) = mtao::geometry::mesh::read_objD(filename);
            auto bb = mtao::geometry::bounding_box(V);
            mtao::Vec3d mean = (bb.min() + bb.max())/2;
            V.colwise() -= mean;
            input_mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
            sg = mandoline::tools::SliceGenerator(V,F);
            auto E = mtao::geometry::mesh::boundary_facets(F);

            input_mesh.setEdgeBuffer(E.cast<unsigned int>());

            C.resize(4,V.cols());
            C.topRows<3>() = mtao::geometry::mesh::vertex_normals(V,F).array()/2+.5;
            C.row(3).setConstant(1);
            input_mesh.setColorBuffer(C.cast<float>().eval());


            


            //input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};
            input_mesh.setParent(&root());


            //bbox_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{slice_object,flat_shader, drawables()};
            slice_vcolor = new mtao::opengl::Drawable<Magnum::Shaders::VertexColor3D>{slice_mesh,vcolor_shader, drawables()};
            //slice_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{slice_mesh,phong_shader, drawables()};
            slice_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{slice_mesh,wireframe_shader, wireframe_drawables};
            slice_mesh.setParent(&root());


            //slice_object.setParent(&root());
            update_slice();
        }
        void gui() override {
            if(input_phong) {input_phong->gui("Input Phong");}
            if(slice_phong) {slice_phong->gui("Cell Phong");}
            //if(slice_vcolor) {slice_vcolor->gui("Cell Vertex Color");}
            if(slice_wireframe) {slice_wireframe->gui("Cell Wireframe");}

            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderFloat3("origin", origin.data(),-1,1)
                    ||
                    ImGui::SliderFloat3("direction", direction.data(),-1,1)
                    ;
                if(dirty) {
                    update_slice();
                }
            }
            if(ImGui::Button("Save")) {
                update_slice();
                mtao::geometry::mesh::write_objD(V,F,"output.obj");
            }


        }
        void  draw() override {
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::PolygonOffsetFill);
            Magnum::GL::Renderer::setPolygonOffset(2.f,1.f);
            Window3::draw();

            
            Magnum::GL::Renderer::setPolygonOffset(0,0);
            camera().draw(wireframe_drawables);

        }
    private:
        Magnum::SceneGraph::DrawableGroup3D wireframe_drawables;
        Magnum::Shaders::Phong phong_shader;
        Magnum::Shaders::Flat3D flat_shader;
        Magnum::Shaders::MeshVisualizer wireframe_shader{supportsGeometryShader()?Magnum::Shaders::MeshVisualizer::Flag::Wireframe:Magnum::Shaders::MeshVisualizer::Flag{}};
        Magnum::Shaders::VertexColor3D vcolor_shader;

        //mtao::opengl::objects::BoundingBox<3> slice_object;
        mtao::opengl::objects::Mesh<3> input_mesh;
        mtao::opengl::objects::Mesh<3> slice_mesh;
        mtao::opengl::objects::Mesh<3> cube_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* bbox_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* slice_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* slice_wireframe = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::VertexColor3D>* slice_vcolor = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)

