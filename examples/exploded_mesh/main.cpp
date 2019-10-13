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
#include <mandoline/tools/planar_slicer.hpp>

#include <thread>
#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:

        mandoline::CutCellMesh<3> ccm;
        mandoline::tools::MeshExploder exploder;


        bool show_wireframes = false;

        float scale = 1.1;
        bool do_slice = false;
        mtao::ColVectors<double,4> colors;

        std::map<int,mandoline::tools::SliceGenerator> slicers;
        std::set<int> regions;


        mtao::Vec3f origin = mtao::Vec3f::Zero(), direction = mtao::Vec3f::Unit(2);


        int size() const { return colors.cols(); }

        void update_exploded() {

            for(auto&& [i,r]: mtao::iterator::enumerate(regions)) {
                auto [VV,FF] = exploder.mesh(scale, std::set<int>{r});
                auto& slicer = slicers[r];
                slicer.set_vertices(VV);
                std::cout << "Region: " << r << " : " << VV.cols() << std::endl;

                if(VV.cols() > 0 && FF.cols() > 0) {
                    if(do_slice) {
                        std::tie(VV,FF) = slicer.slice(origin.cast<double>(),direction.cast<double>());
                    }
                    exploded_meshes[i].setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
                }
            }
            set_region_colors();
        }



        void set_region_colors() {
            for(auto&& [i,r]: mtao::iterator::enumerate(regions)) {
                auto C = exploder.colors(colors, std::set<int>{r}).cast<float>().eval();
                if(do_slice) {
                auto& slicer = slicers[r];

                Eigen::SparseMatrix<float> BM = slicer.barycentric_map().cast<float>();
                mtao::ColVecs4f C2 = C * BM.transpose();
                exploded_meshes[i].setColorBuffer(C2.cast<float>().eval());
                } else {
                exploded_meshes[i].setColorBuffer(C.cast<float>().eval());
                }
            }
        }


        MeshViewer(const Arguments& args): Window3(args) {
            mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
            mtao::logging::make_logger("profiler").set_level(mtao::logging::Level::Off);
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");

            ccm = mandoline::CutCellMesh<3>::from_proto(filename);
            ccm.triangulate_faces();
            input_mesh.setTriangleBuffer(ccm.origV().cast<float>(),ccm.origF().cast<unsigned int>());





            auto R = ccm.regions();
            colors.resize(4,R.size());
            colors.setZero();
            colors.topRows<3>().setRandom();
            colors.row(3).setConstant(1);


            regions.clear();

            std::copy(R.begin(),R.end(), std::inserter(regions,regions.end()));

            for(auto&& r: regions) {
                std::cout << "Region size: " << r << ": " << std::count(R.begin(),R.end(),r) << std::endl;
            }

            exploder = mandoline::tools::MeshExploder(ccm);

            input_mesh.setParent(&root());



            exploded_meshes = decltype(exploded_meshes)(regions.size());
            for(auto&& [i,r]: mtao::iterator::enumerate(regions)) {
                exploded_wireframes.push_back(new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{exploded_meshes[i],wireframe_shader, wireframe_drawables});
                exploded_vcolors.push_back(new mtao::opengl::Drawable<Magnum::Shaders::VertexColor3D>{exploded_meshes[i],vcolor_shader, drawables()});
                exploded_meshes[i].setParent(&root());
                auto [VV,FF] = exploder.mesh(scale, std::set<int>{r});
                slicers[r] = mandoline::tools::SliceGenerator(VV,FF);
            }
            //input_phong = new mtao::opengl::Drawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};

            show_wireframes = true;
            direction(2) = -1;
            origin(2) = -.1;

            update_exploded();
            set_region_colors();




        }
        bool show_multi = false;
        int index = 0;
        void gui() override {

            if(ImGui::Checkbox("Slice", &do_slice)) {
                update_exploded();
            }
            if(ImGui::Checkbox("Show Multi", &show_multi)) {
            }
            if(show_multi) {
                for(auto&& [i,c,w]: mtao::iterator::enumerate(exploded_vcolors, exploded_wireframes)) {
                    c->set_visibility(true);
                    w->set_visibility(true);
                }
            } else {
                ImGui::InputInt("Region id", &index);
                for(auto&& [i,c,w]: mtao::iterator::enumerate(exploded_vcolors, exploded_wireframes)) {
                    if(i == index) {
                        c->set_visibility(true);
                        w->set_visibility(true);
                    } else {
                        c->set_visibility(false);
                        w->set_visibility(false);
                    }
                }
            }

            if(ImGui::Checkbox("Show Wireframes", &show_wireframes)) {
                update_exploded();

            }

            if(ImGui::Button("Reset colors")) {
                colors.topRows<3>().setRandom();
                set_region_colors();
            }
            if(ImGui::Button("Region colors")) {
                colors.topRows<3>().setRandom();
                auto R = ccm.regions();
                double ME = *std::max_element(R.begin(),R.end());
                for(int i = 0; i < R.size(); ++i) {
                    auto c = colors.col(i).topRows<3>();
                    c.setZero();
                    auto r = R[i];
                    double s =  r / ME;
                    r = r % 7;
                    switch(r) {
                        //case 0: c = mtao::Vec3d(0,0,0); break;
                        case 0: c = mtao::Vec3d(s,0,0); break;
                        case 1: c = mtao::Vec3d(0,s,0); break;
                        case 2: c = mtao::Vec3d(0,0,s); break;
                        case 3: c = mtao::Vec3d(s,s,0); break;
                        case 4: c = mtao::Vec3d(s,0,s); break;
                        case 5: c = mtao::Vec3d(0,s,s); break;
                        case 6: c = mtao::Vec3d(s,s,s); break;

                    }
                }
                set_region_colors();
            }
            if(input_phong) {input_phong->gui("Input Phong");}
            //if(exploded_wireframe) {exploded_wireframe->gui("Exploded Wireframe");}
            //if(exploded_wireframe) {exploded_wireframe->gui("Cell Wireframe");}

            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderFloat3("origin", origin.data(),-1,1);
                dirty |= ImGui::SliderFloat3("direction", direction.data(),-1,1);
                if(dirty) {
                    update_exploded();
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

            Window3::draw();
            
            Magnum::GL::Renderer::setPolygonOffset(0,0);
            if(show_wireframes) {
                camera().draw(wireframe_drawables);
            }

        }
    private:
        Magnum::SceneGraph::DrawableGroup3D wireframe_drawables;
        Magnum::Shaders::Phong phong_shader;
        Magnum::Shaders::VertexColor3D vcolor_shader;
        Magnum::Shaders::MeshVisualizer wireframe_shader{supportsGeometryShader()?Magnum::Shaders::MeshVisualizer::Flag::Wireframe:Magnum::Shaders::MeshVisualizer::Flag{}};
        mtao::opengl::objects::Mesh<3> input_mesh;
        std::vector<mtao::opengl::objects::Mesh<3>> exploded_meshes;
        mtao::opengl::Drawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        std::vector<mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>*> exploded_wireframes;
        std::vector<mtao::opengl::Drawable<Magnum::Shaders::VertexColor3D>*> exploded_vcolors;

};




MAGNUM_APPLICATION_MAIN(MeshViewer)

