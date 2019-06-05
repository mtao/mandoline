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
#include <Magnum/EigenIntegration/Integration.h>

#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:

        mandoline::CutCellMesh<3> ccm;
        mtao::ColVecs3d V;
        mtao::ColVecs3i F;

        float scale = 1.0;
        mtao::ColVectors<double,4> colors;

        mandoline::tools::MeshExploder exploder;
        std::set<int> nonzero_regions;

        mtao::Vec3f origin, direction = mtao::Vec3f::Unit(2);


        int size() const { return colors.cols(); }

        void update_exploded() {

            //auto [VV,FF] = exploder.mesh(scale, nonzero_regions);
            //if(VV.cols() > 0 && FF.cols() > 0) {
            //    exploded_mesh.setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
            //}
            for(int i = 0; i < size(); ++i) {
                auto [VV,FF] = exploder.mesh(i,scale);
                if(VV.cols() > 0 && FF.cols() > 0) {

                    exploded_meshes[i].setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
                    auto& C = exploded_wireframes.back()->data().color;
                    auto& c = colors.col(i).cast<float>().eval();
                    C = Magnum::Math::Vector4<float>(c.x(),c.y(),c.z(),c(3));

                    //= Magnum::Math::Vector4<float>{colors.col(i).cast<float>().eval()};
                }

            }
        }





        MeshViewer(const Arguments& args): Window3(args), wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe} {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            ccm = mandoline::CutCellMesh<3>::from_proto(filename);
            ccm.triangulate_faces();
            input_mesh.setTriangleBuffer(ccm.origV().cast<float>(),ccm.origF().cast<unsigned int>());
            auto E = mtao::geometry::mesh::boundary_facets(F);

            input_mesh.setEdgeBuffer(E.cast<unsigned int>());


            exploder = mandoline::tools::MeshExploder(ccm);

            auto R = ccm.regions();
            colors.resize(4,R.size());
            colors.setZero();
            colors.row(3).setConstant(1);

            exploded_meshes = decltype(exploded_meshes)(size());
            for(auto&& [i,r]: mtao::iterator::enumerate(R)) {
                exploded_wireframes.push_back(new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_meshes[i],wireframe_shader, wireframe_drawables});
                std::cout << i << ":" << r << " ";
                colors.col(i)(r) = 1;
                exploded_meshes[i].setParent(&root());
            }
            std::cout << std::endl;


            std::copy(R.begin(),R.end(), std::inserter(nonzero_regions,nonzero_regions.end()));

            if(nonzero_regions.find(0) != nonzero_regions.end()) {
                nonzero_regions.erase(0);
            }


            input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};
            input_mesh.setParent(&root());


            update_exploded();
            //exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{exploded_mesh,phong_shader, drawables()};
            //exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_mesh,wireframe_shader, wireframe_drawables};
            //exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_mesh,wireframe_shader, wireframe_drawables};
            //exploded_mesh.setParent(&root());


        }
        void gui() override {
            if(input_phong) {input_phong->gui("Input Phong");}
            //if(exploded_wireframe) {exploded_wireframe->gui("Exploded Wireframe");}
            //if(exploded_wireframe) {exploded_wireframe->gui("Cell Wireframe");}

            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderFloat("scale", &scale,1,10);
                if(dirty) {
                    update_exploded();
                }
            }


        }
        void  draw() override {
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::PolygonOffsetFill);
            Magnum::GL::Renderer::setPolygonOffset(2.f,1.f);

            camera().draw(wireframe_drawables);
            
            Magnum::GL::Renderer::setPolygonOffset(0,0);
            Window3::draw();

        }
    private:
        Magnum::SceneGraph::DrawableGroup3D wireframe_drawables;
        Magnum::Shaders::Phong phong_shader;
        Magnum::Shaders::MeshVisualizer wireframe_shader;
        mtao::opengl::objects::Mesh<3> input_mesh;
        mtao::opengl::objects::Mesh<3> exploded_mesh;
        std::vector<mtao::opengl::objects::Mesh<3>> exploded_meshes;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* exploded_wireframe = nullptr;
        std::vector<mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>*> exploded_wireframes;

};




MAGNUM_APPLICATION_MAIN(MeshViewer)

