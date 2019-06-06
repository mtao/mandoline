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
#include <mandoline/tools/planar_slicer.hpp>

#include <thread>
#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:

        mandoline::CutCellMesh<3> ccm;
        mtao::ColVecs3d V;
        mtao::ColVecs3i F;
        mtao::ColVecs3d VV;
        mtao::ColVecs3i FF;

        float scale = 1.1;
        mtao::ColVectors<double,4> colors;

        mandoline::tools::SliceGenerator slicer;
        mandoline::tools::MeshExploder exploder;
        std::vector<int> regions;
        std::set<int> nonzero_regions;
        std::map<int,Magnum::Color4> region_colors;

        mtao::Vec3f origin = mtao::Vec3f::Zero(), direction = mtao::Vec3f::Unit(2);


        int size() const { return colors.cols(); }

        void update_exploded() {

            std::tie(VV,FF) = exploder.mesh(scale, nonzero_regions);
            if(VV.cols() > 0 && FF.cols() > 0) {
                auto [VVV,FFF] = slicer.slice(origin.cast<double>(),direction.cast<double>());
                exploded_mesh.setTriangleBuffer(VVV.cast<float>(), FFF.cast<unsigned int>());
            }
            std::mutex mut;
            int i = 0;
//#pragma omp parallel for
            for(i = 0; i < size(); ++i) {
                auto [VV,FF] = exploder.mesh(i,scale);
                if(VV.cols() > 0 && FF.cols() > 0) {

                    //auto [VVV,FFF] = mandoline::tools::slice(VV,FF,origin.cast<double>(),direction.cast<double>());
                    auto&& [VVV,FFF] = std::tie(VV,FF);
                    if(FFF.cols() > 0) {

                        std::lock_guard<std::mutex> lg(mut);
                        exploded_wireframes[i]->set_visibility(true);
                        exploded_meshes[i].setTriangleBuffer(VVV.cast<float>(), FFF.cast<unsigned int>());
                    } else {
                        exploded_wireframes[i]->set_visibility(false);
                    }
                    //exploded_meshes[i].setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
                    auto& C = exploded_wireframes[i]->data().color;
                    auto& c = colors.col(i).cast<float>().eval();
                    C = Magnum::Math::Vector4<float>(c.x(),c.y(),c.z(),c(3));

                    //= Magnum::Math::Vector4<float>{colors.col(i).cast<float>().eval()};
                }

            }
        }



        void set_region_colors() {
            for(auto&& [DC,r]: mtao::iterator::zip(exploded_wireframes, regions)) {
                if(auto it = region_colors.find(r); it != region_colors.end()) {
                    DC->data().color = it->second;
                } else {
                    DC->data().color = Magnum::Math::Vector4<float>(0,0,0,0);
                }

            }
            auto C = exploder.colors(colors, nonzero_regions).cast<float>().eval();
            exploded_mesh.setColorBuffer(C);

            
        }


        MeshViewer(const Arguments& args): Window3(args), wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe} {
            mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
            mtao::logging::make_logger("profiler").set_level(mtao::logging::Level::Off);
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            ccm = mandoline::CutCellMesh<3>::from_proto(filename);
            ccm.triangulate_faces();
            input_mesh.setTriangleBuffer(ccm.origV().cast<float>(),ccm.origF().cast<unsigned int>());
            auto E = mtao::geometry::mesh::boundary_facets(F);

            input_mesh.setEdgeBuffer(E.cast<unsigned int>());




            auto R = ccm.regions();
            regions = R;
            colors.resize(4,R.size());
            colors.setZero();
            colors.row(3).setConstant(1);

            exploded_meshes = decltype(exploded_meshes)(size());
            for(auto&& [i,r]: mtao::iterator::enumerate(R)) {
                exploded_wireframes.push_back(new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_meshes[i],wireframe_shader, wireframe_drawables});
                colors.col(i)(r) = 1;
                exploded_meshes[i].setParent(&root());
            }


            std::copy(R.begin(),R.end(), std::inserter(nonzero_regions,nonzero_regions.end()));

            if(nonzero_regions.find(0) != nonzero_regions.end()) {
                //nonzero_regions.erase(0);
            }

            exploder = mandoline::tools::MeshExploder(ccm);
            std::tie(VV,FF) = exploder.mesh(scale, nonzero_regions);
            slicer = mandoline::tools::SliceGenerator(VV,FF);

            input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};
            input_mesh.setParent(&root());


            update_exploded();
            //exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{exploded_mesh,phong_shader, drawables()};
            //exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_mesh,wireframe_shader, wireframe_drawables};
            exploded_wireframe = new mtao::opengl::Drawable<Magnum::Shaders::VertexColor3D>{exploded_mesh,vcolor_shader, drawables()};
            exploded_mesh.setParent(&root());

            set_region_colors();

        }
        void gui() override {
            if(input_phong) {input_phong->gui("Input Phong");}
            //if(exploded_wireframe) {exploded_wireframe->gui("Exploded Wireframe");}
            //if(exploded_wireframe) {exploded_wireframe->gui("Cell Wireframe");}

            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderFloat3("origin", origin.data(),-1,1)
                    ||
                    ImGui::SliderFloat3("direction", direction.data(),-1,1)
                    ;
                if(dirty) {
                    auto [VVV,FFF] = slicer.slice(origin.cast<double>(),direction.cast<double>());
                    exploded_mesh.setTriangleBuffer(VVV.cast<float>(), FFF.cast<unsigned int>());
                    set_region_colors();
                }
            }
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

            //camera().draw(wireframe_drawables);
            
            Magnum::GL::Renderer::setPolygonOffset(0,0);
            Window3::draw();

        }
    private:
        Magnum::SceneGraph::DrawableGroup3D wireframe_drawables;
        Magnum::Shaders::Phong phong_shader;
        Magnum::Shaders::VertexColor3D vcolor_shader;
        Magnum::Shaders::MeshVisualizer wireframe_shader;
        mtao::opengl::objects::Mesh<3> input_mesh;
        mtao::opengl::objects::Mesh<3> exploded_mesh;
        std::vector<mtao::opengl::objects::Mesh<3>> exploded_meshes;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::VertexColor3D>* exploded_wireframe = nullptr;
        std::vector<mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>*> exploded_wireframes;

};




MAGNUM_APPLICATION_MAIN(MeshViewer)

