#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/bounding_box.hpp"
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
#include <Magnum/Shaders/MeshVisualizer.h>
#include <mtao/opengl/drawables.h>
#include <mtao/opengl/objects/mesh.h>
#include <Corrade/Utility/Arguments.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include "slice_generator.h"

#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:
        int current_ccm_cell = 0;
        int current_ccm_face = 0;
        std::optional<std::string> filename;

        mandoline::CutCellMesh<3> ccm;
        mtao::ColVecs3d V;
        mtao::ColVecs3i F;



        mtao::Vec3f origin, direction = mtao::Vec3f::Unit(2);



        void update_slice() {
            SliceGenerator sg(origin.cast<double>(),direction.cast<double>());
            auto [VV,FF] = sg.slice(V,F);
            if(FF.cols() > 0) {

                slice_mesh.setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
            }
        }





        MeshViewer(const Arguments& args): Window3(args), wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe} {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            std::tie(V,F) = mtao::geometry::mesh::read_objD(filename);
            auto bb = mtao::geometry::bounding_box(V);
            mtao::Vec3d mean = (bb.min() + bb.max())/2;
            V.colwise() -= mean;
            input_mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
            auto E = mtao::geometry::mesh::boundary_facets(F);

            input_mesh.setEdgeBuffer(E.cast<unsigned int>());




            input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};
            input_mesh.setParent(&root());


            update_slice();
            slice_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{slice_mesh,phong_shader, drawables()};
            slice_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{slice_mesh,wireframe_shader, wireframe_drawables};
            slice_mesh.setParent(&root());


        }
        void gui() override {
            if(input_phong) {input_phong->gui("Input Phong");}
            if(slice_phong) {slice_phong->gui("Cell Phong");}
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
        Magnum::Shaders::MeshVisualizer wireframe_shader;
        mtao::opengl::objects::Mesh<3> input_mesh;
        mtao::opengl::objects::Mesh<3> slice_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* slice_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* slice_wireframe = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)

