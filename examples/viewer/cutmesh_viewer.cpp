#include <Corrade/Utility/Arguments.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/Phong.h>
#include <mtao/eigen/stack.h>
#include <mtao/opengl/drawables.h>
#include <mtao/opengl/objects/mesh.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <mtao/cmdline_parser.hpp>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/types.hpp>

#include "imgui.h"
#include "mandoline/mesh3.hpp"
#include "mtao/opengl/Window.h"
using namespace mtao::logging;

class MeshViewer : public mtao::opengl::Window3 {
   public:
    int current_ccm_cell = 0;
    int current_ccm_face = 0;
    std::optional<std::string> filename;

    mandoline::CutCellMesh<3> ccm;

    void update_face() {
        if (ccm.faces().empty()) {
            return;
        }
        current_ccm_face =
            std::clamp<int>(current_ccm_face, 0, ccm.faces().size() - 1);
        auto [V, F] = ccm.compact_triangulated_face(current_ccm_face);

        face_mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
    }

    void update_cell() {
        if (ccm.cells().empty()) {
            return;
        }
        current_ccm_cell =
            std::clamp<int>(current_ccm_cell, 0, ccm.cell_size() - 1);

        auto [V, F] = ccm.compact_triangulated_cell(current_ccm_cell);

        cell_mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
    }

    MeshViewer(const Arguments& args) : Window3(args) {
        Corrade::Utility::Arguments myargs;
        myargs.addArgument("filename").parse(args.argc, args.argv);
        std::string filename = myargs.value("filename");

        ccm = mandoline::CutCellMesh<3>::from_proto(filename);
        ccm.triangulate_faces();
        auto bbox = ccm.bbox();

        edges_mesh.setEdgeBuffer((ccm.vertices()).cast<float>(),
                                 ccm.edges().cast<unsigned int>());
        input_mesh.setTriangleBuffer((ccm.origV()).cast<float>(),
                                     ccm.origF().cast<unsigned int>());
        face_centroid_mesh.setVertexBuffer(
            ccm.face_centroids().cast<float>().eval());
        cell_centroid_mesh.setVertexBuffer(
            ccm.cell_centroids().cast<float>().eval());

        input_phong = new mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{
            input_mesh, phong_shader, drawables()};
        input_wireframe =
            new mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>{
                input_mesh, wireframe_shader, wireframe_drawables};
        input_mesh.setParent(&root());

        update_cell();
        update_face();
        cell_phong = new mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{
            cell_mesh, phong_shader, drawables()};
        cell_wireframe =
            new mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>{
                cell_mesh, wireframe_shader, wireframe_drawables};
        cell_mesh.setParent(&root());

        face_phong = new mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>{
            face_mesh, phong_shader, drawables()};
        face_wireframe =
            new mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>{
                face_mesh, wireframe_shader, wireframe_drawables};
        face_mesh.setParent(&root());

        cell_flat = new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>{
            cell_centroid_mesh, flat_shader, drawables()};
        cell_centroid_mesh.setParent(&root());
        face_flat = new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>{
            face_centroid_mesh, flat_shader, drawables()};
        face_centroid_mesh.setParent(&root());

        face_flat->deactivate();
        cell_flat->deactivate();
        face_flat->activate_points();
        cell_flat->activate_points();

        edges_mesh.setParent(&root());
        edges_flat = new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>{
            edges_mesh, flat_shader, drawables()};
        edges_flat->deactivate();
        edges_flat->activate_edges();

        face_flat->data().color[0] = 1;
        face_flat->data().color[1] = 0;
        face_flat->data().color[2] = 0;
        face_flat->data().color[3] = 1;
        cell_flat->data().color[0] = 0;
        cell_flat->data().color[1] = 1;
        cell_flat->data().color[2] = 0;
        cell_flat->data().color[3] = 1;
        for (auto&& r : {input_wireframe, cell_wireframe, face_wireframe}) {
            r->data().color[3] = 0;
            r->data().color[0] = 0;
            r->data().color[1] = 0;
            r->data().color[2] = 0;
        }
    }
    void gui() override {
        if (input_phong) {
            input_phong->gui("Input Phong");
        }
        if (input_wireframe) {
            input_wireframe->gui("Input Wireframe");
        }
        if (cell_phong) {
            cell_phong->gui("Cell Phong");
        }
        if (face_phong) {
            face_phong->gui("Face Phong");
        }
        if (cell_wireframe) {
            cell_wireframe->gui("Cell Wireframe");
        }
        if (face_wireframe) {
            face_wireframe->gui("Face Wireframe");
        }

        if (cell_flat) {
            cell_flat->gui("Cell Flat");
        }
        if (face_flat) {
            face_flat->gui("Face Flat");
        }
        auto&& io = ImGui::GetIO();

        {
            bool dirty = ImGui::SliderInt("Face position", &current_ccm_face, 0,
                                          ccm.faces().size() - 1);
            if (dirty) {
                update_face();
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
        camera().draw(wireframe_drawables);
    }

   private:
    Magnum::SceneGraph::DrawableGroup3D wireframe_drawables;
    Magnum::Shaders::Phong phong_shader;
    Magnum::Shaders::Flat3D flat_shader;
    Magnum::Shaders::MeshVisualizer3D wireframe_shader{
        supportsGeometryShader()
            ? Magnum::Shaders::MeshVisualizer3D::Flag::Wireframe
            : Magnum::Shaders::MeshVisualizer3D::Flag{}};
    mtao::opengl::objects::Mesh<3> input_mesh;
    mtao::opengl::objects::Mesh<3> edges_mesh;
    mtao::opengl::objects::Mesh<3> cell_mesh;
    mtao::opengl::objects::Mesh<3> face_mesh;
    mtao::opengl::objects::Mesh<3> cell_centroid_mesh;
    mtao::opengl::objects::Mesh<3> face_centroid_mesh;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* input_phong = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* cell_phong = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Phong>* face_phong = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>*
        input_wireframe = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>*
        cell_wireframe = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>*
        face_wireframe = nullptr;

    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>* edges_flat = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>* cell_flat = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>* face_flat = nullptr;
};

MAGNUM_APPLICATION_MAIN(MeshViewer)
