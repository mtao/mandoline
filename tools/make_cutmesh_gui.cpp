#include "mtao/opengl/Window.h"
#include <iostream>
#include <Eigen/Geometry>
#include "imgui.h"
#include <memory>
#include <algorithm>
#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/bounding_box.hpp"
#include "mtao/opengl/drawables.h"
#include <mtao/types.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/ArrayView.h>

#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/GL/Renderer.h>
#include <Corrade/Utility/Arguments.h>
#include <mtao/opengl/objects/grid.h>
#include <mtao/geometry/grid/grid.h>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mandoline/construction/generator.hpp>
#include <mandoline/mesh3.hpp>

#include <glm/gtc/matrix_transform.hpp> 

using namespace mtao::opengl;


class MeshViewer: public mtao::opengl::Window3 {
    public:
        enum class Mode: int { Smoothing, LSReinitialization };
        Mode mode = Mode::LSReinitialization;

        char output_filename[128] = "output.cutmesh";
        std::string_view output_view;
        float permeability = 100.0;
        float timestep = 1000.0;
        bool animate = false;
        using Vec = mtao::VectorX<GLfloat>;
        using Vec3 = mtao::Vec3f;
        Vec data;
        Vec data_original;
        Vec dx;
        Vec dy;
        Vec signs;
        Eigen::AlignedBox<float,3> bbox, orig_bbox;
        Eigen::SparseMatrix<float> L;
        bool use_cube = false;
        bool adaptive=true;
        int adaptive_level=0;
        mtao::ColVecs3d V;
        mtao::ColVecs3i F;
        std::vector<std::array<int,2>> edges;
        std::map<std::array<int,2>,std::vector<std::array<int,2>>> mapped_edges;
        std::optional<std::array<int,2>> edge_choice = {};

        std::array<int,3> N{{5,5,5}};
        int& NI=N[0];
        int& NJ=N[1];
        int& NK=N[2];
        std::optional<mandoline::CutCellMesh<3>> ccm;


        MeshViewer(const Arguments& args): Window3(args), _wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe}, output_view(output_filename,mtao::types::container_size(output_filename)) {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");
            std::tie(V,F) = mtao::geometry::mesh::read_objD(filename);
                std::cout << "V/E/F " << V.cols() << "/" << mtao::geometry::mesh::boundary_facets(F).cols() << "/" << F.cols() << std::endl;
            mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
            orig_bbox = bbox = mtao::geometry::bounding_box(V.cast<float>().eval());
            mtao::Vec3f trans_e= -((bbox.min() + bbox.max()) / 2).cast<float>();
            Magnum::Math::Vector3<float> trans(trans_e.x(),trans_e.y(),trans_e.z());
            mesh.translate(trans);
            edge_mesh.translate(trans);
            grid.translate(trans);
            float s = 1./bbox.sizes().maxCoeff();
            mesh.scale(Magnum::Math::Vector3<float>(s));
            edge_mesh.scale(Magnum::Math::Vector3<float>(s));
            grid.scale(Magnum::Math::Vector3<float>(s));
            mv_drawable = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{mesh,_wireframe_shader, drawables()};


            edge_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{grid,_flat_shader, drawables()};
            edge_drawable->activate_triangles({});
            edge_drawable->activate_edges();

            edge_boundary_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{edge_mesh,_flat_shader, edge_drawables};
            edge_boundary_drawable->activate_triangles({});
            edge_boundary_drawable->activate_edges();

            edge_boundary_drawable->data().color = Magnum::Color4(1,0,0,1);


            edge_mesh.setParent(&root());
            edge_boundary_drawable->set_visibility(false);
            grid.setParent(&root());

            mesh.setParent(&root());
            update();
        }
        void update() {
            //mtao::geometry::grid::Grid3f g(std::array<int,3>{{NI,NJ,NK}});
            auto sg = mtao::geometry::grid::StaggeredGrid3f::from_bbox
                (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
            std::cout << "cell shape: " << mtao::eigen::stl2eigen(sg.cell_shape()).transpose() << std::endl;
            std::cout << sg.cell_grid().size() << std::endl;
            std::cout << "vertex shape: " << mtao::eigen::stl2eigen(sg.vertex_shape()).transpose() << std::endl;
            std::cout << sg.vertex_grid().size() << std::endl;
            //auto g = mtao::geometry::grid::Grid3f::from_bbox
            //    (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
            auto g = sg.vertex_grid();

            std::cout << "vertex shape: " << mtao::eigen::stl2eigen(g.shape()).transpose() << std::endl;
            grid.set(g);


        }
        void update_edges() {
            if(edge_choice) {
                if(auto it = mapped_edges.find(*edge_choice); it != mapped_edges.end()) {
                    auto&& edges = it->second;
                    if(edges.size() > 0) {
                        auto E = mtao::eigen::stl2eigen(edges);

                        edge_mesh.setEdgeBuffer(E.cast<unsigned int>().eval());

                        edge_boundary_drawable->set_visibility(true);
                        return;
                    }

                }
            }
            if(edges.size() > 0) {
                auto E = mtao::eigen::stl2eigen(edges);

                edge_mesh.setEdgeBuffer(E.cast<unsigned int>().eval());

                edge_boundary_drawable->set_visibility(true);
                return;
            }
            edge_boundary_drawable->set_visibility(false);
        }
        void gui() override {
            if(ImGui::InputInt3("N", &NI))  {
                update();
            }
            if(ImGui::InputFloat3("min", bbox.min().data()))  {
                bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(),bbox.max());
                update();
            }
            if(ImGui::InputFloat3("max", bbox.max().data()))  {
                bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(),bbox.max());
                update();
            }
            if(ImGui::Checkbox("Cubes", &use_cube)) {
                update();
            }
            if(ImGui::Checkbox("Adaptive", &adaptive)) {
                update();
            }
            if(ImGui::InputInt("Adaptive level", &adaptive_level))  {
                adaptive = true;

                update();
            }

            if(mv_drawable) {
                mv_drawable->gui();
            }
            if(edge_drawable) {
                edge_drawable->gui();
            }
            if(ImGui::Button("Reset BBox")) {
                bbox = orig_bbox;
                update();
            }
            if(ImGui::Button("Make CCM")) {
                Eigen::AlignedBox<double,3> bbox(this->bbox.min().cast<double>(),this->bbox.max().cast<double>());
                auto sg = mtao::geometry::grid::StaggeredGrid3d::from_bbox
                    (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
                auto ccg = mandoline::construction::CutCellGenerator<3>(V,sg, {});
                ccg.add_boundary_elements(F);
                ccg.adaptive = adaptive;
                if(adaptive) {
                    ccg.adaptive_level = adaptive_level;
                }
                ccg.bake();
                std::vector<mtao::ColVecs3i> Fs;
                ccm = ccg.generate();
                mtao::ColVecs3f V = ccm->vertices().cast<float>();
                for(auto&& f: ccm->faces()) {
                    if(f.is_mesh_face()) {
                        Fs.push_back(f.triangulate_fan());
                    }
                }
                auto F = mtao::eigen::hstack_iter(Fs.begin(),Fs.end()).cast<unsigned int>().eval();
                mesh.setTriangleBuffer(V,F);
                edge_mesh.setVertexBuffer(V);


                using E = std::array<int,2>;
                edges.clear();
                mapped_edges.clear();
                for(int i = 0; i < ccm->cut_edge_size(); ++i) {
                    auto e = ccm->cut_edge(i);
                    auto mask = ccg.grid_vertex(e(0)).mask() & ccg.grid_vertex(e(1)).mask();
                    if(mask.active()) {
                        edges.emplace_back(E{{e(0),e(1)}});
                        for(int i = 0; i < 3; ++i) {
                            if(mask[i]) {
                                mapped_edges[E{{i,*mask[i]}}].emplace_back(edges.back());
                            }
                        }
                    }
                }
                update_edges();

            }
            {
                bool active = bool(edge_choice);
                if(ImGui::Checkbox("Edge choice", &active)) {
                    std::cout << "Checked!" << std::endl;
                    if(active) {
                        edge_choice = std::array<int,2>{{0,0}};
                    } else {
                        edge_choice.reset();
                        update_edges();
                    }
                }
                if(active) {
                    bool changed= false;
                    changed |= ImGui::InputInt("Axis", edge_choice->begin());
                    changed |= ImGui::InputInt("Plane", edge_choice->begin()+1);
                    if(changed) {
                        auto&& [a,b] = *edge_choice;
                        a =std::clamp<int>(a,0,2);
                        b = std::clamp<int>(b,0,ccm->vertex_shape()[a] - 1);
                        update_edges();
                    }
                }
            }
            ImGui::InputText("Filename",output_filename,128);
            if(ImGui::Button("Save")) {
                if(ccm) {
                    ccm->write(std::string(output_view));

                }
            }
        }
        void draw() override {
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Window3::draw();
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::DepthTest);
            camera().draw(edge_drawables);
            Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::DepthTest);
        }
    private:
        Magnum::Shaders::MeshVisualizer _wireframe_shader;
        Magnum::Shaders::Flat3D _flat_shader;
        Magnum::SceneGraph::DrawableGroup3D edge_drawables;
        mtao::opengl::objects::Mesh<3> mesh;
        mtao::opengl::objects::Grid<3> grid;
        mtao::opengl::objects::Mesh<3> edge_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* mv_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* edge_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* edge_boundary_drawable = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)
