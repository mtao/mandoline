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
#include <igl/read_triangle_mesh.h>

#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/GL/Renderer.h>
#include <Corrade/Utility/Arguments.h>
#include <mtao/opengl/objects/grid.h>
#include <mtao/geometry/grid/grid.h>
#include <mtao/geometry/grid/staggered_grid.hpp>
#include <mtao/geometry/prune_vertices.hpp>
#include <mandoline/construction/construct.hpp>
#include <mandoline/construction/remesh_self_intersections.hpp>
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
        float bbox_offset = .1;
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
        std::optional<mandoline::construction::DeformingGeometryConstructor> constructor;
        std::optional<mandoline::CutCellMesh<3>> ccm;

            auto staggered_grid() const { 
                return mtao::geometry::grid::StaggeredGrid3f::from_bbox
                (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
            }

        MeshViewer(const Arguments& args): Window3(args), output_view(output_filename,mtao::types::container_size(output_filename)) {
            Corrade::Utility::Arguments myargs;
            myargs.addArgument("filename").parse(args.argc,args.argv);
            std::string filename = myargs.value("filename");
            Eigen::MatrixXd VV;
            Eigen::MatrixXi FF;
            igl::read_triangle_mesh(filename,VV,FF);

            //std::tie(V,F) = mtao::geometry::mesh::read_objD(filename);
            V = VV.transpose();
            F = FF.transpose();
                std::cout << "V/E/F " << V.cols() << "/" << mtao::geometry::mesh::boundary_facets(F).cols() << "/" << F.cols() << std::endl;
            mesh.setTriangleBuffer(V.cast<float>(),F.cast<unsigned int>());
            orig_bbox = bbox = mtao::geometry::bounding_box(V.cast<float>().eval());
            mtao::Vec3f trans_e= -((bbox.min() + bbox.max()) / 2).cast<float>();
            Magnum::Math::Vector3<float> trans(trans_e.x(),trans_e.y(),trans_e.z());

            //std::tie(V,F) = mtao::geometry::prune(V,F,0);
            //std::tie(V,F) = mandoline::construction::remesh_self_intersections(V,F);

            constructor.emplace(V,F,staggered_grid());



            mesh.translate(trans);
            edge_mesh.translate(trans);
            vertex_mesh.translate(trans);
            grid.translate(trans);
            float s = 1./bbox.sizes().maxCoeff();
            mesh.scale(Magnum::Math::Vector3<float>(s));
            edge_mesh.scale(Magnum::Math::Vector3<float>(s));
            vertex_mesh.scale(Magnum::Math::Vector3<float>(s));
            grid.scale(Magnum::Math::Vector3<float>(s));
            mv_drawable = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{mesh,_wireframe_shader, drawables()};


            edge_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{grid,_flat_shader, drawables()};
            edge_drawable->activate_triangles({});
            edge_drawable->activate_edges();

            edge_boundary_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{edge_mesh,_flat_shader, edge_drawables};
            point_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{vertex_mesh,_flat_shader, edge_drawables};
            edge_boundary_drawable->activate_triangles({});
            edge_boundary_drawable->activate_edges();

            edge_boundary_drawable->data().color = Magnum::Color4(1,0,0,1);


            edge_mesh.setParent(&root());
            vertex_mesh.setVertexBuffer(mtao::Vec3f(1,1,1));
            point_drawable->deactivate();
            point_drawable->activate_points();
            vertex_mesh.setParent(&root());
            edge_boundary_drawable->set_visibility(false);
            grid.setParent(&root());

            mesh.setParent(&root());
            update();
        }
        void update() {
            auto sg = staggered_grid();
            //mtao::geometry::grid::Grid3f g(std::array<int,3>{{NI,NJ,NK}});
            constructor->update_grid(sg);
            constructor->update_vertices(V);
            auto g = sg.vertex_grid();

            grid.set(g);


        }
        void update_edges() {
            mtao::ColVecs3f VV = ccm->vertices().cast<float>();
            if(edge_choice) {
                if(auto it = mapped_edges.find(*edge_choice); it != mapped_edges.end()) {
                    auto&& edges = it->second;
                    if(edges.size() > 0) {
                        auto E = mtao::eigen::stl2eigen(edges);

                        edge_mesh.setEdgeBuffer(E.cast<unsigned int>().eval());
                std::set<int> Vidx;
                std::copy(E.data(),E.data()+E.size(),std::inserter(Vidx,Vidx.end()));
                mtao::ColVecs3f V(3,Vidx.size());
                for(auto&& [i,j]: mtao::iterator::enumerate(Vidx)) {
                    V.col(i) = VV.col(j).cast<float>();
                }
                vertex_mesh.setVertexBuffer(V);

                        edge_boundary_drawable->set_visibility(true);
                        return;
                    }

                }
            }
            if(edges.size() > 0) {
                auto E = mtao::eigen::stl2eigen(edges);
                vertex_mesh.setVertexBuffer(VV.cast<float>());

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
            if(ImGui::InputInt("Adaptive level", &adaptive_level))  {

                constructor->set_adaptivity(adaptive_level);
                update();
            }

            if(mv_drawable) {
                mv_drawable->gui();
            }
            if(edge_drawable) {
                edge_drawable->gui();
            }
            if(ImGui::InputFloat("Boundingbox Offset", &bbox_offset))  {
                mtao::Vec3f C = (orig_bbox.min() + orig_bbox.max()) / 2;
                mtao::Vec3f s = orig_bbox.sizes() / 2;
                bbox.min() = C - (1 + bbox_offset) * s;
                bbox.max() = C + (1 + bbox_offset) * s;
                update();

            }
            if(ImGui::Button("Reset BBox")) {
                bbox = orig_bbox;
                update();
            }
            if(ImGui::Button("Make CCM")) {
                constructor->bake();
                ccm = constructor->emit();
                std::vector<mtao::ColVecs3i> Fs;
                mtao::ColVecs3f V = ccm->vertices().cast<float>();
                for(auto&& f: ccm->faces()) {
                    if(f.is_mesh_face()) {
                        Fs.push_back(f.triangulate_fan());
                    }
                }
                if(Fs.size() > 0) {
                    auto F = mtao::eigen::hstack_iter(Fs.begin(),Fs.end()).cast<unsigned int>().eval();
                    mesh.setTriangleBuffer(V,F);
                    edge_mesh.setVertexBuffer(V);


                    using E = std::array<int,2>;
                    edges.clear();
                    mapped_edges.clear();
                    for(int i = 0; i < ccm->cut_edge_size(); ++i) {
                        auto e = ccm->cut_edge(i);
                        auto va = ccm->masked_vertex(e(0));
                        auto vb = ccm->masked_vertex(e(1));
                        auto mask = va.mask() & vb.mask();
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
            Magnum::GL::Renderer::setPointSize(10.);
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Window3::draw();
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::DepthTest);
            camera().draw(edge_drawables);
            Magnum::GL::Renderer::enable(Magnum::GL::Renderer::Feature::DepthTest);
        }
    private:
        Magnum::Shaders::MeshVisualizer _wireframe_shader{supportsGeometryShader()?Magnum::Shaders::MeshVisualizer::Flag::Wireframe:Magnum::Shaders::MeshVisualizer::Flag{}};
        Magnum::Shaders::Flat3D _flat_shader;
        Magnum::SceneGraph::DrawableGroup3D edge_drawables;
        mtao::opengl::objects::Mesh<3> mesh;
        mtao::opengl::objects::Grid<3> grid;
        mtao::opengl::objects::Mesh<3> edge_mesh;
        mtao::opengl::objects::Mesh<3> vertex_mesh;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* mv_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* edge_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* edge_boundary_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* point_drawable = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)
