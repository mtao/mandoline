#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/mesh/write_ply.hpp"
#include "mtao/geometry/bounding_box.hpp"
#include <mtao/types.hpp>
#include <igl/colormap.h>
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
#include "setup.h"
#include <igl/parula.h>
#include <igl/colormap.h>

#include <thread>
#include "mandoline/mesh3.hpp"
using namespace mtao::logging;




class MeshViewer: public mtao::opengl::Window3 {
    public:

        mandoline::CutCellMesh<3> ccm;
        mandoline::tools::MeshExploder exploder;
        mtao::Vec3d wind_direction = mtao::Vec3d::Unit(1);

        std::optional<int> current_picked_face = {};
        mtao::VecXd face_pick_flux;

        float strip_rate = .05;
        float strip_size = .05;

        int multi_index = 0;
        bool show_multi = false;

        float scale = 1.1;
        bool do_slice = true;
        mtao::ColVectors<double,4> colors;
        mtao::ColVectors<double,3> cachedV;
        mtao::ColVectors<double,3> cachedC;
        mtao::ColVectors<int,3> cachedF;


        std::map<int,mandoline::tools::SliceGenerator> slicers;
        std::set<int> regions;


        mtao::Vec3f origin = mtao::Vec3f::Zero(), direction = mtao::Vec3f::Unit(2);


        int size() const { return colors.cols(); }

        void update_exploded() {

            for(auto&& [i,r]: mtao::iterator::enumerate(regions)) {
                auto [VV,FF] = exploder.mesh(scale, std::set<int>{r});
                auto& slicer = slicers[r];
                slicer.set_vertices(VV);
                if(!do_slice) {
                    if(i == multi_index) {
                        cachedV = VV;
                        cachedF = FF;
                    }
                    exploded_meshes[i].setTriangleBuffer(VV.cast<float>(), FF.cast<unsigned int>());
                }
            }
        }
        void update_slice() {

            for(auto&& [i,r]: mtao::iterator::enumerate(regions)) {
                auto& slicer = slicers[r];
                if(do_slice) {
                    auto [VV,FF] = slicer.slice(origin.cast<double>(),direction.cast<double>());
                if(show_multi) {
                    if(i == multi_index) {
                        cachedV = VV;
                        cachedF = FF;
                    }
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
                if(show_multi) {
                    if(i == multi_index) {
                        cachedC = C2.topRows<3>().cast<double>();
                    }
                }
                } else {
                exploded_meshes[i].setColorBuffer(C.cast<float>().eval());
                if(show_multi) {
                    if(i == multi_index) {
                        cachedC = C.topRows<3>().cast<double>();
                    }
                }
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

            face_pick_flux = mtao::VecXd::Zero(ccm.face_size());

            auto bb = ccm.bbox();
            mtao::Vec3d C = (bb.min() + bb.max())/2;

            root().translate(Magnum::Math::Vector3<float>(-C.x(),-C.y(),-C.z()));



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
                exploded_vcolors.push_back(new mtao::opengl::MeshDrawable<Magnum::Shaders::VertexColor3D>{exploded_meshes[i],vcolor_shader, drawables()});
                exploded_meshes[i].setParent(&root());
                auto [VV,FF] = exploder.mesh(scale, std::set<int>{r});
                slicers[r] = mandoline::tools::SliceGenerator(VV,FF);
            }
            input_phong = new mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{input_mesh,phong_shader, drawables()};

            direction(2) = -1;
            origin(2) = -.1;

            update_exploded();

            update_slice();



        }
        void gui() override {

            if(ImGui::Checkbox("Slice", &do_slice)) {
                update_slice();
            }
            bool multi_dirty = false;
            if(ImGui::Checkbox("Show Multi", &show_multi)) {
                for(auto&& [i,c]: mtao::iterator::enumerate(exploded_vcolors) ) {
                    c->set_visibility(true);
                    multi_dirty = true;
                }
            }
            if(show_multi) {
                if(ImGui::InputInt("Region id", &multi_index) || multi_dirty) {
                    for(auto&& [i,c]: mtao::iterator::enumerate(exploded_vcolors )) {
                        if(i == multi_index) {
                            c->set_visibility(true);
                        } else {
                            c->set_visibility(false);
                        }
                    }
                }
            }

            if(ImGui::Button("Random colors")) {
                colors.topRows<3>().setRandom();
                set_region_colors();
            }

            if(ImGui::Button("Region colors")) {
                colors.topRows<3>().setRandom();
                auto R = ccm.regions();
                for(int i = 0; i < R.size(); ++i) {
                    auto c = colors.col(i).topRows<3>();
                    c.setZero();
                    auto r = R[i];
                    switch(r) {
                        case 0: c = mtao::Vec3d(0,0,0); break;
                        case 1: c = mtao::Vec3d(1,0,0); break;
                        case 2: c = mtao::Vec3d(0,1,0); break;
                        case 3: c = mtao::Vec3d(0,0,1); break;
                        case 4: c = mtao::Vec3d(1,1,0); break;
                        case 5: c = mtao::Vec3d(1,0,1); break;
                        case 6: c = mtao::Vec3d(0,1,1); break;
                        case 7: c = mtao::Vec3d(1,1,1); break;

                    }
                }
                set_region_colors();
            }
            if(ImGui::Button("Boundary colors")) {
                mtao::VecXd C = divergence(ccm,wind_direction);
                //C.array() -= C.minCoeff();
                C /= C.cwiseAbs().maxCoeff();
                colors.row(0) = C.transpose();
                colors.row(1) = -C.transpose().array();
                colors.row(2).array() = 0;
                colors.row(3).array() = 1;
                set_region_colors();
            }

            {
                mtao::Vec3f dir = wind_direction.cast<float>();
                ImGui::SliderFloat3("wind direction", dir.data(),-1,1);
                wind_direction = dir.cast<double>();
            }

            if(ImGui::Button("Poisson colors")) {
                mtao::VecXd C = pressure(ccm,wind_direction);
                C.array() /= C.maxCoeff();

                face_pick_flux = C;
                //C.array() -= (C.minCoeff() + C.maxCoeff()) / 2;
                //C /= C.cwiseAbs().maxCoeff();
                //colors.row(0) = C.transpose();
                //colors.row(1) = -C.transpose().array();
                //colors.row(2).array() = 0;
                Eigen::MatrixXd COL;
                C = (C.array() + 1) / 2;
                /*
                C.array() -= .5;
                C *= 3;
                C.array() += .5;
                */
                igl::parula(C,false,COL);
                COL.colwise() = C;
                colors.topRows<3>() = COL.transpose();
                colors.row(3).array() = 1;
                set_region_colors();
            }
            /*
            bool a = ImGui::InputFloat("Strip size", &strip_size);
            bool b = ImGui::InputFloat("Strip rate", &strip_rate);
            if(a | b) {
                mtao::VecXd d = face_pick_flux;
                d = (strip_rate * d.array()).atan();
                d.array() /= d.maxCoeff();
                d.array() -= d.minCoeff();
                d /= d.maxCoeff();
                std::cout << d.minCoeff() << ":::" << d.maxCoeff() << std::endl;
                d = 1-(d.array()/strip_size*M_PI).array().sin().pow(4).abs().eval();
                d = d.array() * face_pick_flux.array();
                d = (d.array() + 1) / 2;
                //d.array() *= face_pick_flux.array();
                Eigen::MatrixXd COL;
                igl::parula(d,false,COL);
                COL.setConstant(1);
                COL.array().colwise() *= d.array();
                igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,d,false,COL);
                colors.topRows<3>() = COL.transpose();
                colors.row(3).array() = 1;
                set_region_colors();

            }
            */
            if(input_phong) {input_phong->gui("Input Phong");}

            auto&& io = ImGui::GetIO();


            {
                bool dirty = ImGui::SliderFloat3("origin", origin.data(),-1,1);
                dirty |= ImGui::SliderFloat3("direction", direction.data(),-1,1);
                if(dirty) {
                    update_slice();
                }
            }
            {
                bool dirty = ImGui::SliderFloat("scale", &scale,1.01,10);
                if(dirty) {
                    update_exploded();
                    if(do_slice) {
                        update_slice();
                    }
                }
            }

            if(ImGui::Button("Save")) {
                update_slice();
                set_region_colors();
                std::cout << cachedV.cols() << " " << cachedC.cols() << std::endl;
                mtao::geometry::mesh::write_plyD(cachedV,cachedC,cachedF,"output.ply");
            }

        }
        void  draw() override {
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::PolygonOffsetFill);
            Magnum::GL::Renderer::setPolygonOffset(2.f,1.f);

            Window3::draw();
            
            Magnum::GL::Renderer::setPolygonOffset(0,0);

        }
        void mousePressEvent(MouseEvent& event) override{
            Window3::mousePressEvent(event);
            if(!ImGui::GetIO().WantCaptureMouse) {
                if(event.button() == MouseEvent::Button::Left) { 
                    auto T = input_mesh.absoluteTransformationMatrix();
                    //std::cout << std::string(T) << std::endl;
                    //update();

                }
            }
        }
    private:
        Magnum::Shaders::Phong phong_shader;
        Magnum::Shaders::VertexColor3D vcolor_shader;
        mtao::opengl::objects::Mesh<3> input_mesh;
        std::vector<mtao::opengl::objects::Mesh<3>> exploded_meshes;
        mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* input_phong = nullptr;
        std::vector<mtao::opengl::MeshDrawable<Magnum::Shaders::VertexColor3D>*> exploded_vcolors;

};




MAGNUM_APPLICATION_MAIN(MeshViewer)

