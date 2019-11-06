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

#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:
        int current_ccm_cell = 0;
        int current_ccm_face = 0;
        std::optional<std::string> filename;

        mandoline::CutCellMesh<3> ccm;


        void update_cell() {
            if(ccm.cells().empty()) {
                return;
            }
            current_ccm_cell = std::clamp<int>(current_ccm_cell,0,ccm.cell_size()-1);

            auto [V,F] = ccm.compact_triangulated_cell(current_ccm_cell);


            cell_mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
        }









        MeshViewer(const Arguments& args): Window3(args) {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            ccm = mandoline::CutCellMesh<3>::from_proto(filename);
            ccm.triangulate_faces();
            auto bbox = ccm.bbox();


            input_mesh.setTriangleBuffer((ccm.origV()).cast<float>(), ccm.origF().cast<unsigned int>());


            input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};
            input_mesh.setParent(&root());


            update_cell();
            cell_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{cell_mesh,phong_shader, drawables()};
            cell_mesh.setParent(&root());


            //point_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{point_mesh,phong_shader, drawables()};
            //point_mesh.setParent(&root());

        }
        void gui() override {
            if(input_phong) {input_phong->gui("Input Phong");}
            if(cell_phong) {cell_phong->gui("Cell Phong");}
            if(point_phong) {cell_phong->gui("Point Phong");}
            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderInt("Cell position",&current_ccm_cell,0,ccm.cell_size()-1);
                if(dirty) {
                    update_cell();
                }
            }


        }
        void  draw() override {
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::PolygonOffsetFill);
            Magnum::GL::Renderer::setPolygonOffset(2.f,1.f);
            Magnum::GL::Renderer::setPointSize(10.);
            Window3::draw();

            
            Magnum::GL::Renderer::setPolygonOffset(0,0);

        }
    private:
        Magnum::Shaders::Phong phong_shader;
        mtao::opengl::objects::Mesh<3> input_mesh;
        mtao::opengl::objects::Mesh<3> cell_mesh;
        mtao::opengl::objects::Mesh<3> point_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* cell_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* point_phong = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)
