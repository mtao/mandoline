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

#include <glm/gtc/matrix_transform.hpp> 

using namespace mtao::opengl;


class MeshViewer: public mtao::opengl::Window3 {
    public:
        enum class Mode: int { Smoothing, LSReinitialization };
        Mode mode = Mode::LSReinitialization;

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

        std::array<int,3> N{{20,20,20}};
        int& NI=N[0];
        int& NJ=N[1];
        int& NK=N[2];


        MeshViewer(const Arguments& args): Window3(args), _wireframe_shader{Magnum::Shaders::MeshVisualizer::Flag::Wireframe} {
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
            grid.translate(trans);
            float s = 1./bbox.sizes().maxCoeff();
            mesh.scale(Magnum::Math::Vector3<float>(s));
            grid.scale(Magnum::Math::Vector3<float>(s));
            mv_drawable = new mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>{mesh,_wireframe_shader, drawables()};


            edge_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{grid,_flat_shader, drawables()};
            edge_drawable->activate_triangles({});
            edge_drawable->activate_edges();
            grid.setParent(&root());

            mesh.setParent(&root());
            update();
        }
        void update() {
            //mtao::geometry::grid::Grid3f g(std::array<int,3>{{NI,NJ,NK}});
            auto sg = mtao::geometry::grid::StaggeredGrid3f::from_bbox
                (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
            //auto g = mtao::geometry::grid::Grid3f::from_bbox
            //    (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
            auto g = sg.vertex_grid();

            grid.set(g);


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
            if(ImGui::Button("Make CCG")) {
                Eigen::AlignedBox<double,3> bbox(this->bbox.min().cast<double>(),this->bbox.max().cast<double>());
                auto sg = mtao::geometry::grid::StaggeredGrid3d::from_bbox
                    (bbox, std::array<int,3>{{NI,NJ,NK}}, use_cube);
                auto ccg = mandoline::construction::CutCellGenerator<3>(V,sg, {});
                ccg.add_boundary_elements(F);
                ccg.bake();
                ccg.adaptive = adaptive;
                if(adaptive) {
                    ccg.adaptive_level = adaptive_level;
                }
                std::vector<mtao::ColVecs3i> Fs;
                auto ccm = ccg.generate();
                mtao::ColVecs3f V = ccm.vertices().cast<float>();
                for(auto&& f: ccm.faces()) {
                    if(f.is_mesh_face()) {
                        Fs.push_back(f.triangulate_fan());
                    }
                }
                auto F = mtao::eigen::hstack_iter(Fs.begin(),Fs.end()).cast<unsigned int>().eval();
                mesh.setTriangleBuffer(V,F);

            }
        }
        void draw() override {
            Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
            Window3::draw();
        }
    private:
        Magnum::Shaders::MeshVisualizer _wireframe_shader;
        Magnum::Shaders::Flat3D _flat_shader;
        mtao::opengl::objects::Mesh<3> mesh;
        mtao::opengl::objects::Grid<3> grid;
        mtao::opengl::Drawable<Magnum::Shaders::MeshVisualizer>* mv_drawable = nullptr;
        mtao::opengl::Drawable<Magnum::Shaders::Flat3D>* edge_drawable = nullptr;


};




MAGNUM_APPLICATION_MAIN(MeshViewer)
