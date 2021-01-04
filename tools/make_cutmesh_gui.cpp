#include "mtao/opengl/Window.h"
#include <iostream>
#include <Eigen/Geometry>
#include "imgui.h"
#include <memory>
#include <algorithm>
#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/mesh/read_obj.hpp"
#include "mtao/geometry/bounding_box.hpp"
#include "mtao/geometry/mesh/write_obj.hpp"
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
#include <mandoline/construction/construct_gui_widget.hpp>
#include <mandoline/construction/remesh_self_intersections.hpp>
#include <mandoline/construction/preprocess_mesh.hpp>
#include <mandoline/mesh3.hpp>
#include <mandoline/tools/cutmesh_info.hpp>
#include "validation/cutmesh_validation.hpp"
#include <spdlog/spdlog.h>


using namespace mtao::opengl;


class MeshViewer : public mtao::opengl::Window3 {
  public:
    enum class Mode : int { Smoothing,
                            LSReinitialization };
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
    Eigen::AlignedBox<float, 3> bbox, orig_bbox;
    Eigen::SparseMatrix<float> L;
    bool use_cube = false;
    int adaptive_level = 0;
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    std::vector<std::array<int, 2>> edges;
    std::map<std::array<int, 2>, std::vector<std::array<int, 2>>> mapped_edges;
    std::optional<std::array<int, 2>> edge_choice = {};

    std::optional<mandoline::CutCellMesh<3>> ccm;


    MeshViewer(const Arguments &args) : Window3(args), output_view(output_filename, mtao::types::container_size(output_filename)), constructor(_flat_shader, drawables()) {
        Corrade::Utility::Arguments myargs;
        myargs.addArgument("filename").parse(args.argc, args.argv);
        std::string filename = myargs.value("filename");
        Eigen::MatrixXd VV;
        Eigen::MatrixXi FF;
        igl::read_triangle_mesh(filename, VV, FF);


        //std::tie(V,F) = mtao::geometry::mesh::read_objD(filename);
        V = VV.transpose();
        F = FF.transpose();
        std::cout << "V/E/F " << V.cols() << "/" << mtao::geometry::mesh::boundary_facets(F).cols() << "/" << F.cols() << std::endl;
        mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
        orig_bbox = bbox = mtao::geometry::bounding_box(V.cast<float>().eval());
        mtao::Vec3f trans_e = -((bbox.min() + bbox.max()) / 2).cast<float>();
        Magnum::Math::Vector3<float> trans(trans_e.x(), trans_e.y(), trans_e.z());

        std::tie(V, F) = mtao::geometry::prune(V, F, 0);
        std::tie(V, F) = mandoline::construction::preprocess_mesh(V, F);

        {
            auto bbox = mtao::geometry::bounding_box(V);
            bbox = mtao::geometry::expand_bbox(bbox, 1.1);
            auto s = bbox.sizes().eval();

            for (int i = 0; i < 3; ++i) {
                if (s(i) < 1e-5) {
                    bbox.min()(i) -= 1e-2;
                    bbox.max()(i) += 1e-2;
                }
            }
            auto g = mtao::geometry::grid::StaggeredGrid3d::from_bbox(bbox, constructor.N);
            constructor.set_mesh(V,F,{});
            constructor.update_mesh_and_bbox(V,F,1.1);

        }

        float s = 1. / bbox.sizes().maxCoeff();
        mesh.scale(Magnum::Math::Vector3<float>(s));
        edge_mesh.scale(Magnum::Math::Vector3<float>(s));
        vertex_mesh.scale(Magnum::Math::Vector3<float>(s));
        constructor.scale(Magnum::Math::Vector3<float>(s));
        mv_drawable = new mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D>{ mesh, _wireframe_shader, drawables() };



        point_drawable = new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>{ vertex_mesh, _flat_shader, edge_drawables };
        point_drawable->deactivate();
        point_drawable->activate_points();

        edge_boundary_drawable = new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D>{ edge_mesh, _flat_shader, edge_drawables };
        edge_boundary_drawable->set_visibility(false);
        edge_boundary_drawable->activate_triangles({});
        edge_boundary_drawable->activate_edges();

        edge_boundary_drawable->data().color = Magnum::Color4(1, 0, 0, 1);

        //center_point.translate(trans);

        vertex_mesh.setVertexBuffer(mtao::Vec3f(1, 1, 1));

        center_point.setParent(&root());
        edge_mesh.setParent(&center_point);
        vertex_mesh.setParent(&center_point);
        constructor.setParent(&center_point);

        mesh.setParent(&center_point);
        update();
    }
    void update() {
        //mtao::geometry::grid::Grid3f g(std::array<int,3>{{NI,NJ,NK}});
        constructor.update_grid();
        constructor.update_vertices(V);
    }
    void update_edges() {
        mtao::ColVecs3f VV = ccm->vertices().cast<float>();
        if (edge_choice) {
            if (auto it = mapped_edges.find(*edge_choice); it != mapped_edges.end()) {
                auto &&edges = it->second;
                if (edges.size() > 0) {
                    auto E = mtao::eigen::stl2eigen(edges);

                    edge_mesh.setEdgeBuffer(E.cast<unsigned int>().eval());
                    std::set<int> Vidx;
                    std::copy(E.data(), E.data() + E.size(), std::inserter(Vidx, Vidx.end()));
                    mtao::ColVecs3f V(3, Vidx.size());
                    for (auto &&[i, j] : mtao::iterator::enumerate(Vidx)) {
                        V.col(i) = VV.col(j).cast<float>();
                    }
                    vertex_mesh.setVertexBuffer(V);

                    edge_boundary_drawable->set_visibility(true);
                    return;
                }
            }
        }
        if (edges.size() > 0) {
            auto E = mtao::eigen::stl2eigen(edges);
            vertex_mesh.setVertexBuffer(VV.cast<float>());

            edge_mesh.setEdgeBuffer(E.cast<unsigned int>().eval());

            edge_boundary_drawable->set_visibility(true);
            return;
        }
        edge_boundary_drawable->set_visibility(false);
    }

    void make_ccm() {
        ccm = constructor.generate();
        if (ccm) {
            spdlog::info("Made CCM!!");
            mandoline::tools::print_region_info(*ccm);

#if defined(MANDOLINE_HANDLE_SELF_INTERSECTIONS)
            auto regions_vec = input_mesh_regions(*ccm);
            std::set<int> regions;
            std::copy(regions_vec.data(), regions_vec.data() + regions_vec.size(), std::inserter(regions, regions.end()));
            std::vector<std::set<std::array<int,3>>> Fs(regions.size());
            for(int i = 0; i < regions_vec.cols(); ++i) {
                auto R = regions_vec.col(i);
                auto f = ccm->origF().col(i);
                for(int j = 0; j < regions_vec.rows(); ++j) {
                    Fs[R(j)].emplace(std::array<int,3>{{f(0),f(1),f(2)}});
                }
            }

            spdlog::info("Input mesh had {} regions", regions.size());
            std::copy(regions.begin(),regions.end(),std::ostream_iterator<int>(std::cout,","));
            std::cout << std::endl;
            for(auto&& [idx,fs]: mtao::iterator::enumerate(Fs)) {
                std::stringstream name;
                name << "objfile" << idx << ".obj";
                mtao::geometry::mesh::write_objD(ccm->origV(),mtao::eigen::stl2eigen(fs), name.str());
            }
#endif
        }
        {// check exterior grid valences
    spdlog::warn("gui side Exterior grid face size: {}", ccm->exterior_grid().num_faces());
            if(!exterior_cell_valence_counts(*ccm)) {
                std::cout << "Bad exterior cell valences!" << std::endl;
            }
        }
        {
            auto g = ccm->Base::vertex_grid();
            auto s = g.shape();
        }

        std::vector<mtao::ColVecs3i> Fs;
        mtao::ColVecs3f V = ccm->vertices().cast<float>();
        for (auto &&f : ccm->faces()) {
            if (f.is_mesh_face()) {
                Fs.push_back(f.triangulate_fan());
            }
        }
        if (Fs.size() > 0) {
            auto F = mtao::eigen::hstack_iter(Fs.begin(), Fs.end()).cast<unsigned int>().eval();
            mesh.setTriangleBuffer(V, F);
            edge_mesh.setVertexBuffer(V);


            using E = std::array<int, 2>;
            edges.clear();
            mapped_edges.clear();
            for (int i = 0; i < ccm->cut_edge_size(); ++i) {
                auto e = ccm->cut_edge(i);
                auto va = ccm->masked_vertex(e.indices[0]);
                auto vb = ccm->masked_vertex(e.indices[1]);
                auto mask = va.mask() & vb.mask();
                if (mask.active()) {
                    edges.emplace_back(e.indices);
                    for (int i = 0; i < 3; ++i) {
                        if (mask[i]) {
                            mapped_edges[E{ { i, *mask[i] } }].emplace_back(edges.back());
                        }
                    }
                }
            }
            update_edges();
        }
    }
    void gui() override {

        if (mv_drawable) {
            mv_drawable->gui();
        }
        if (ImGui::InputFloat("Boundingbox Offset", &bbox_offset)) {
            mtao::Vec3f C = (orig_bbox.min() + orig_bbox.max()) / 2;
            mtao::Vec3f s = orig_bbox.sizes() / 2;
            bbox.min() = C - (1 + bbox_offset) * s;
            bbox.max() = C + (1 + bbox_offset) * s;
            update();
        }
        if (ImGui::Button("Reset BBox")) {
            bbox = orig_bbox;
            update();
        }
        if (ImGui::CollapsingHeader("Cutmesh Generation")) {
            if (constructor.gui()) {
                update();
            }
            if (ImGui::Button("Make CCM")) {
                make_ccm();
            }
        }
        {
            bool active = bool(edge_choice);
            if (ImGui::Checkbox("Edge choice", &active)) {
                std::cout << "Checked!" << std::endl;
                if (active) {
                    edge_choice = std::array<int, 2>{ { 0, 0 } };
                } else {
                    edge_choice.reset();
                    update_edges();
                }
            }
            if (active) {
                bool changed = false;
                changed |= ImGui::InputInt("Axis", edge_choice->begin());
                changed |= ImGui::InputInt("Plane", edge_choice->begin() + 1);
                if (changed) {
                    auto &&[a, b] = *edge_choice;
                    a = std::clamp<int>(a, 0, 2);
                    b = std::clamp<int>(b, 0, ccm->vertex_shape()[a] - 1);
                    update_edges();
                }
            }
        }
        ImGui::InputText("Filename", output_filename, 128);
        if (ImGui::Button("Save")) {
            if (ccm) {
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
    Magnum::Shaders::MeshVisualizer3D _wireframe_shader{ supportsGeometryShader() ? Magnum::Shaders::MeshVisualizer3D::Flag::Wireframe : Magnum::Shaders::MeshVisualizer3D::Flag{} };
    Magnum::Shaders::Flat3D _flat_shader;
    Magnum::SceneGraph::DrawableGroup3D edge_drawables, dummy;
    mtao::opengl::Object3D center_point;
    mandoline::construction::CutmeshGeneratorGui constructor;
    mtao::opengl::objects::Mesh<3> mesh;
    mtao::opengl::objects::Mesh<3> edge_mesh;
    mtao::opengl::objects::Mesh<3> vertex_mesh;
    mtao::opengl::MeshDrawable<Magnum::Shaders::MeshVisualizer3D> *mv_drawable = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D> *edge_boundary_drawable = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat3D> *point_drawable = nullptr;
};


MAGNUM_APPLICATION_MAIN(MeshViewer)
