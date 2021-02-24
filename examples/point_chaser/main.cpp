#include <Corrade/Utility/Arguments.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Phong.h>
#include <mandoline/operators/nearest_facet.hpp>
#include <mtao/eigen/stack.h>
#include <mtao/opengl/drawables.h>
#include <mtao/opengl/objects/bbox.h>
#include <mtao/opengl/objects/mesh.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <mtao/cmdline_parser.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/mesh/shapes/sphere.hpp>
#include <mtao/types.hpp>

#include "imgui.h"
#include "mandoline/mesh3.hpp"
#include "mtao/opengl/Window.h"
using namespace mtao::logging;

class MeshViewer : public mtao::opengl::Window3 {
   public:
    bool single_particle_mode = false;
    bool nearest_mode = true;
    bool triangle_cut_face_mode = true;
    int current_ccm_cell = 0;
    int current_ccm_face = 0;
    std::optional<std::string> filename;
    mtao::Vec3d p;
    mtao::ColVecs4d cell_colors;
    mtao::ColVecs3d P;

    mandoline::CutCellMesh<3> ccm;

    void update_cell() {
        if (ccm.cells().empty()) {
            return;
        }
        current_ccm_cell =
            std::clamp<int>(current_ccm_cell, 0, ccm.cell_size() - 1);

        auto [V, F] = ccm.compact_triangulated_cell(current_ccm_cell);

        cell_mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
    }
    void update_points() {
        auto bb = ccm.bbox();

        if (single_particle_mode) {
            p = bb.sample();
            if (nearest_mode) {
                current_ccm_cell = ccm.get_nearest_cell_index(p);
            } else {
                current_ccm_cell = ccm.get_cell_index(p);
            }
            std::cout << p.transpose() << ")) " << current_ccm_cell
                      << std::endl;

            point_mesh.setVertexBuffer(p.cast<float>().eval());
            // point_mesh.setVertexBuffer(mtao::ColVecs3f::Random(3,
            // 100).eval());
            mtao::Vec4f C;
            if (current_ccm_cell >= 0) {
                C = cell_colors.col(current_ccm_cell).cast<float>();

            } else {
                C.setZero();
            }
            point_mesh.setColorBuffer(C);
        } else {
            mtao::ColVecs4f C(4, P.cols());
            if (nearest_mode || triangle_cut_face_mode) {
                mtao::Vec3d s = bb.sizes();
                bb.min() -= .5 * s;
                bb.max() += .5 * s;
            }
            for (int j = 0; j < P.cols(); ++j) {
                auto p = P.col(j);
                p = bb.sample();
            }
            mtao::VecXi I;
            if (nearest_mode) {
                I = ccm.get_nearest_cell_indices(P);
            } else if (triangle_cut_face_mode) {
                I = mandoline::operators::nearest_mesh_cut_faces(ccm,P);
            } else {
                I = ccm.get_cell_indices(P);
            }
            for (int j = 0; j < P.cols(); ++j) {
                int cell = I(j);

                // point_mesh.setVertexBuffer(mtao::ColVecs3f::Random(3,
                // 100).eval());
                auto c = C.col(j);
                if (cell >= 0) {
                    c = cell_colors.col(cell).cast<float>();

                } else {
                    c.setZero();
                }
            }
            point_mesh.setVertexBuffer(P.cast<float>().eval());
            point_mesh.setColorBuffer(C);
        }
    }

    MeshViewer(const Arguments& args) : Window3(args) {
        Corrade::Utility::Arguments myargs;
        myargs.addArgument("filename").parse(args.argc, args.argv);
        std::string filename = myargs.value("filename");

        ccm = mandoline::CutCellMesh<3>::from_proto(filename);
        ccm.triangulate_faces();
        auto bb = ccm.bbox();
        bbox.set_bbox(bb.cast<float>());
        bbox_drawable =
            mtao::opengl::make_drawable(bbox, _flat_shader, drawables());
        bbox.setParent(&root());
        bbox_drawable->deactivate();
        bbox_drawable->activate_edges();
        P.resize(3, 10000);

        cell_colors.resize(4, std::max<int>(ccm.num_cells(),ccm.num_faces()));
        cell_colors.setRandom();
        cell_colors.array() += 1;
        cell_colors /= 2;
        cell_colors.row(3).setOnes();

        input_mesh.setTriangleBuffer((ccm.origV()).cast<float>(),
                                     ccm.origF().cast<unsigned int>());

        input_phong = new mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{
            input_mesh, phong_shader, drawables()};
        input_mesh.setParent(&root());

        update_cell();
        cell_phong = new mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{
            cell_mesh, phong_shader, drawables()};
        cell_mesh.setParent(&root());
        // cell_phong->deactivate();
        // cell_phong->activate_edges();

        update_points();
        point_flat =
            new mtao::opengl::MeshDrawable<Magnum::Shaders::VertexColor3D>{
                point_mesh, vertex_color_shader, drawables()};
        point_flat->point_size = 10;
        input_mesh.setParent(&root());
        point_mesh.setParent(&root());
        point_flat->deactivate();
        point_flat->activate_points();
        // point_phong = new
        // mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{point_mesh,phong_shader,
        // drawables()}; point_mesh.setParent(&root());
    }
    void gui() override {
        if (ImGui::Checkbox("Single particle mode", &single_particle_mode)) {
            update_points();
        }

        if (ImGui::Checkbox("Nearest mode", &nearest_mode)) {
            update_points();
        }
        if (ImGui::Checkbox("Triangle-cut-face", &triangle_cut_face_mode)) {
            update_points();
        }
        if (input_phong) {
            input_phong->gui("Input Phong");
        }
        if (cell_phong) {
            cell_phong->gui("Cell Phong");
        }
        if (point_flat) {
            point_flat->gui("Point Flat");
        }
        auto&& io = ImGui::GetIO();

        {
            if (ImGui::Button("Random point")) {
                update_points();
                update_cell();
            }
        }

        {
            bool dirty = ImGui::SliderInt("Cell position", &current_ccm_cell, 0,
                                          ccm.cell_size() - 1);
            if (dirty) {
                update_cell();
            }
        }
    }
    void draw() override {
        Magnum::GL::Renderer::disable(
            Magnum::GL::Renderer::Feature::FaceCulling);
        Magnum::GL::Renderer::enable(
            Magnum::GL::Renderer::Feature::PolygonOffsetFill);
        Magnum::GL::Renderer::setPolygonOffset(2.f, 1.f);
        Magnum::GL::Renderer::setPointSize(10.);
        Window3::draw();

        Magnum::GL::Renderer::setPolygonOffset(0, 0);
    }

   private:
    Magnum::Shaders::Phong phong_shader;
    Magnum::Shaders::VertexColor3D vertex_color_shader;
    Magnum::Shaders::Flat3D _flat_shader;
    mtao::opengl::objects::Mesh<3> input_mesh;
    mtao::opengl::objects::Mesh<3> cell_mesh;
    mtao::opengl::objects::Mesh<3> point_mesh;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* input_phong = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* cell_phong = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* point_color = nullptr;

    mtao::opengl::MeshDrawable<Magnum::Shaders::VertexColor3D>* point_flat =
        nullptr;
    mtao::opengl::objects::BoundingBox<3> bbox;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>* bbox_drawable =
        nullptr;
};

MAGNUM_APPLICATION_MAIN(MeshViewer)
