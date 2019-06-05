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
#include <mandoline/tools/exploded_mesh.hpp>

#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:

        mandoline::CutCellMesh<3> ccm;
        mtao::ColVecs3d V;
        mtao::ColVecs3i F;

        float scale = 1.0;

        mandoline::tools::MeshExploder exploder;

        mtao::Vec3f origin, direction = mtao::Vec3f::Unit(2);



        void update_exploded() {
            for(auto&& d: exploded_wireframes) {
                wireframe_drawables.remove(*d);
            }

            auto [VV,FF] = exploder.mesh(scale);
            if(FF.cols() > 0) {

                exploded_mesh.setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
            }
        }





        MeshViewer(const Arguments& args): Window3(args), wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe} {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            input_mesh.setTriangleBuffer(ccm.origV().cast<float>(),ccm.origF().cast<unsigned int>());
            auto E = mtao::geometry::mesh::boundary_facets(F);

            input_mesh.setEdgeBuffer(E.cast<unsigned int>());


            exploder = mandoline::tools::MeshExploder(ccm);




            input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};
            input_mesh.setParent(&root());


            update_exploded();
            exploded_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{exploded_mesh,phong_shader, drawables()};
            exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_mesh,wireframe_shader, wireframe_drawables};
            exploded_mesh.setParent(&root());


        }
        void gui() override {
            if(input_phong) {input_phong->gui("Input Phong");}
            if(exploded_phong) {exploded_phong->gui("Cell Phong");}
            if(exploded_wireframe) {exploded_wireframe->gui("Cell Wireframe");}

            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderFloat3("origin", origin.data(),-1,1)
                    ||
                    ImGui::SliderFloat3("direction", direction.data(),-1,1)
                    ;
                if(dirty) {
                    update_exploded();
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
        mtao::opengl::objects::Mesh<3> exploded_mesh;
        std::vector<mtao::opengl::objects::Mesh<3>> exploded_meshes;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* exploded_phong = nullptr;
        std::vector<mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>*> exploded_wireframes;

};




MAGNUM_APPLICATION_MAIN(MeshViewer)

