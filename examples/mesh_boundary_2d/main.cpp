#include <Corrade/Utility/Arguments.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <Magnum/GL/Renderer.h>
#include <mtao/eigen/stack.h>
#include <mtao/opengl/drawables.h>
#include <mtao/opengl/objects/bbox.h>
#include <mtao/opengl/objects/grid.h>
#include <mtao/opengl/objects/mesh.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/geometry/mesh/stack_meshes.hpp>
#include <mtao/solvers/linear/preconditioned_conjugate_gradient.hpp>
#include <mtao/types.hpp>
#include <thread>

#include "imgui.h"
#include "mandoline/construction/tools/json_to_cutmesh2.hpp"
#include "mandoline/mesh3.hpp"
#include "mandoline/operators/diffgeo2.hpp"
#include "mandoline/operators/region_boundaries2.hpp"
#include "mtao/opengl/Window.h"
using namespace mtao::logging;
using namespace Magnum::Math::Literals;

class MeshViewer : public mtao::opengl::Window2 {
   public:
    std::optional<mandoline::CutCellMesh<2>> ccm;
    bool show_multi = false;
    int index = 0;
    std::optional<int> active_region_index;

    mtao::opengl::Vector2 cursor;

    bool show_wireframes = false;

    float region_center_scale = 0.0;
    bool do_slice = false;
    mtao::ColVectors<double, 4> colors;
    mtao::ColVectors<double, 2> cachedV;
    mtao::ColVectors<int, 2> cachedE;
    mtao::Vec2f harmonic_dir = mtao::Vec2f::UnitX();

    int face_index = -1;

    mtao::Vec2f origin = mtao::Vec2f::Zero(), direction = mtao::Vec2f::Unit(1);

    int size() const { return colors.cols(); }

    MeshViewer(const Arguments &args) : Window2(args) {
        cutcell_drawable =
            new mtao::opengl::MeshDrawable<Magnum::Shaders::VertexColor2D>{
                cutcell_mesh, vcolor_shader, background_drawgroup};
        cutcell_mesh.setParent(&root());
        cutcell_drawable->deactivate();

        edge_drawable = new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>{
            edge_mesh, _flat_shader, curve_drawgroup};
        edge_mesh.setParent(&root());
        edge_drawable->deactivate();

        cutcell_face_drawable =
            new mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>{
                cutcell_face_mesh, _flat_shader, curve_drawgroup};
        cutcell_face_mesh.setParent(&root());
        cutcell_face_drawable->deactivate();

        Corrade::Utility::Arguments myargs;
        if (args.argc > 1) {
            myargs.addArgument("filename").parse(args.argc, args.argv);
            std::string filename = myargs.value("filename");
            spdlog::warn("Reading file [{}]", filename);
            ccm = mandoline::construction::tools::json_to_cutmesh2(filename);
            update_ccm();
            update_edges();
        }
    }

    void draw() override;
    void gui() override;
    void update_edges();
    void update_faces();
    void update_ccm();
    void update_colors();  // doesnt upload the colors to anything

    void mouseMoveEvent(MouseMoveEvent &event) override;
    void mousePressEvent(MouseEvent &event) override;
    void update_face(int idx);

   private:
    Magnum::SceneGraph::DrawableGroup2D background_drawgroup, curve_drawgroup;

    Magnum::Shaders::Flat2D _flat_shader;
    Magnum::Shaders::VertexColor2D vcolor_shader;
    mtao::opengl::objects::Mesh<2> edge_mesh;
    mtao::opengl::objects::Mesh<2> cutcell_mesh;
    mtao::opengl::objects::Mesh<2> cutcell_face_mesh;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D> *edge_drawable =
        nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::VertexColor2D>
        *cutcell_drawable = nullptr;
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D> *cutcell_face_drawable =
        nullptr;
};
void MeshViewer::draw() {
    Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::DepthTest);
    Magnum::GL::Renderer::disable(Magnum::GL::Renderer::Feature::FaceCulling);
    Magnum::GL::Renderer::setPointSize(10);

    camera().draw(background_drawgroup);
    Window2::draw();
    camera().draw(curve_drawgroup);
}

void MeshViewer::mouseMoveEvent(MouseMoveEvent &event) {
    Window2::mouseMoveEvent(event);
    cursor = localPosition(event.position());
}
void MeshViewer::mousePressEvent(MouseEvent &event) {
    Window2::mousePressEvent(event);
}

void MeshViewer::gui() {
    {
        bool per_region_boundary = bool(active_region_index);
        if (ImGui::Checkbox("Per Region Boundary", &per_region_boundary)) {
            if (per_region_boundary) {
                active_region_index = 0;
            } else {
                active_region_index = {};
            }
            update_edges();
        }
        if (per_region_boundary) {
            if (ImGui::InputInt("Region index", &*active_region_index)) {
                update_edges();
            }
        }
    }
    {
        if (ImGui::Button("Update Color")) {
            update_faces();
        }
    }

    if (ccm) {
        ImGui::Text("Cursor position (%f,%f) cell: %d", cursor.x(), cursor.y(),
                    face_index);
    } else {
        ImGui::Text("Cursor position (%f,%f)", cursor.x(), cursor.y());
    }
}

void MeshViewer::update_edges() {
    if (!ccm) {
        spdlog::info("Cutmesh not found, can't update cutmesh vis");
        return;
    }
    std::set<int> boundaries;
    if (!active_region_index) {
        spdlog::info("Getting all region boundaries");
        boundaries = mandoline::operators::region_boundaries(*ccm);
    } else {
        spdlog::info("Getting region {} boundaries", *active_region_index);
        auto b =
            mandoline::operators::region_boundaries(*ccm, *active_region_index);
        std::transform(b.begin(), b.end(),
                       std::inserter(boundaries, boundaries.end()),
                       [](auto &&v) { return std::get<0>(v); });
    }

    auto E = ccm->edges();
    mtao::ColVectors<uint32_t, 2> EE(2, boundaries.size());

    for (auto [idx, eidx] : mtao::iterator::enumerate(boundaries)) {
        EE.col(idx) = E.col(eidx).cast<uint32_t>();
    }

    // set_colors(ccm.active_grid_cell_mask);

    {
        auto V = ccm->vertices();
        edge_mesh.setEdgeBuffer(V.cast<float>().eval(), EE);
        spdlog::info("activating edge drawable");
        edge_drawable->activate_edges();
        spdlog::info("activating edge darwable done");
    }
}

// void MeshViewer::update_grid() {
//    auto grid =
//    mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox.cast<double>(), N);
//    constructor.update_grid(grid);
//
//}
void MeshViewer::update_ccm() {
    if (!ccm) {
        spdlog::info("Cutmesh not found, can't update cutmesh vis");
    }

    // set_colors(ccm.active_grid_cell_mask);

    {
        auto V = ccm->vertices();
        auto E = ccm->cut_edges_eigen();
        // mtao::ColVectors<float, 3> ccg_cols_edge(3, E.cols());
        // ccg_cols_edge.setZero();
        // for (int i = 0; i < ccg_cols_edge.cols(); ++i) {
        //    ccg_cols_edge(i % 3, i) = 1;
        //}
        // mtao::ColVectors<float, 3> ccg_cols(3, V.cols());
        // ccg_cols.setZero();
        // for (int i = 0; i < ccm->StaggeredGrid::vertex_size(); ++i) {
        //    ccg_cols.col(i).setConstant(.7);
        //}
        cutcell_mesh.setEdgeBuffer(V.cast<float>().eval(),
                                   E.cast<unsigned int>().eval());
        cutcell_drawable->deactivate();

        spdlog::warn("Making gird draable");
        // ccm->faces();
        auto r = ccm->regions();
        spdlog::warn("Writing regions");
        std::set<int> regions;
        for (int i = 0; i < r.size(); ++i) {
            regions.insert(r(i));
        }
        std::cout << "Region count: " << regions.size() << std::endl;
    }
    update_faces();
}

void MeshViewer::update_faces() {
    if (!ccm) return;

    spdlog::warn("updating faces");

    update_colors();
    std::vector<std::tuple<mtao::ColVecs3i, mtao::Vec4d>> FCs;
    FCs.resize(ccm->num_cells());
    spdlog::warn("making meshes");
    auto verts = ccm->vertices();
    // tbb::parallel_for(tbb::blocked_range<size_t>(0, ccm->num_cells()),
    // [&](const tbb::blocked_range<size_t> &range) {
    //    for (size_t idx = range.begin(); idx != range.end(); ++idx) {
    for (int idx = 0; idx < ccm->num_cells(); ++idx) {
        // spdlog::info("Getting cell {}", idx);
        auto c = ccm->cell(idx);
        if (c.size() > 0) {
            auto &&d = *c.begin();
            // std::copy(d.begin(),d.end(),std::ostream_iterator<int>(std::cout,","));
            // std::cout << std::endl;

            // std::cout << idx << " => " << mtao::geometry::curve_volume(verts,
            // d) << std::endl;
            auto F = mtao::geometry::mesh::earclipping(verts, d);
            FCs[idx] = std::tuple<mtao::ColVecs3i, mtao::Vec4d>{
                std::move(F), colors.col(idx)};
        }
        // spdlog::info("{} / {}", idx, ccm->num_cells());
    }
    //});

    spdlog::warn("stacking meshes");
    auto [V, F, C] = mtao::geometry::mesh::stack_meshes(ccm->vertices(), FCs);
    spdlog::error("nV/nF/nC: {} {} {}", V.cols(), F.cols(), C.cols());

    cutcell_face_mesh.setVertexBuffer(V.cast<float>().eval());
    cutcell_mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
    cutcell_mesh.setColorBuffer(C.cast<float>().eval());
    cutcell_drawable->activate_triangles();
}
void MeshViewer::update_colors() {
    spdlog::warn("updating colors");
    if (!ccm) return;
    colors.resize(4, ccm->num_cells());
    auto region_colors = [&]() {
        spdlog::warn("computing region colors");
        mtao::VecXi R = ccm->regions();
        int num_regions = R.maxCoeff() + 1;
        mtao::ColVecs4d C(4, num_regions);
        C.setRandom();
        C.row(3).setConstant(1);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, R.size()),
                          [&](const tbb::blocked_range<size_t> &range) {
                              for (size_t idx = range.begin();
                                   idx != range.end(); ++idx) {
                                  colors.col(idx) = C.col(R(idx));
                              }
                          });
    };
    auto random_colors = [&]() {
        spdlog::warn("computing random colors");
        colors.setRandom();
    };

    // region_colors();
    random_colors();
    colors.setConstant(.1);
    spdlog::warn("Picked out colors");
    colors.noalias() = (colors.array() > 0).select(colors, 0);
    colors.row(3).setConstant(1);
}

MAGNUM_APPLICATION_MAIN(MeshViewer)

