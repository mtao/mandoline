#include "mtao/geometry/mesh/boundary_facets.h"
#include "mtao/geometry/bounding_box.hpp"
#include <mtao/types.hpp>
#include <tbb/parallel_for.h>
#include "mtao/opengl/Window.h"
#include <iostream>
#include "imgui.h"
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/stack.h>
#include <memory>
#include <algorithm>
#include <mtao/opengl/drawables.h>
#include <Magnum/EigenIntegration/Integration.h>
#include <mtao/opengl/objects/mesh.h>
#include <mtao/opengl/objects/bbox.h>
#include <mtao/opengl/objects/grid.h>
#include <Corrade/Utility/Arguments.h>
#include "plcurve2.hpp"
#include "mandoline/construction/generator2.hpp"
#include "mandoline/construction/construct_imgui2.hpp"
#include "mandoline/construction/face_collapser.hpp"
#include <Magnum/GL/Renderer.h>
#include <mtao/geometry/mesh/stack_meshes.hpp>
#include <mtao/solvers/linear/preconditioned_conjugate_gradient.hpp>
#include "mandoline/operators/diffgeo2.hpp"

#include <thread>
#include "mandoline/mesh3.hpp"
using namespace mtao::logging;
using namespace Magnum::Math::Literals;


class MeshViewer : public mtao::opengl::Window2 {
  public:
      using Base = CutmeshGenerator2_Imgui;
    std::optional<mandoline::CutCellMesh<2>> ccm;
    std::optional<mandoline::construction::CutmeshGenerator2_Imgui> constructor;
    bool show_multi = false;
    int index = 0;
    PLCurve2 curve;


    mtao::opengl::Vector2 cursor;

    bool show_wireframes = false;

    float region_center_scale = 0.0;
    bool do_slice = false;
    mtao::ColVectors<double, 4> colors;
    mtao::ColVectors<double, 2> cachedV;
    mtao::ColVectors<int, 2> cachedE;
    mtao::Vec2f harmonic_dir = mtao::Vec2f::UnitX();

    int face_index = -1;
    enum ColorType : char {
        Random,
        Regions,
        Harmonic,
        Harmonic_RHS
    };
    enum InterfaceMode : char {
        CurveEdit,
        BBoxResize,
        Browse
    };

    InterfaceMode ui_mode = CurveEdit;
    ColorType color_type = ColorType::Random;

    mtao::Vec2f origin = mtao::Vec2f::Zero(), direction = mtao::Vec2f::Unit(1);


    int size() const { return colors.cols(); }


    MeshViewer(const Arguments &args) : Window2(args) {
        //mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
        mtao::logging::make_logger("profiler").set_level(mtao::logging::Level::Off);
        //Corrade::Utility::Arguments myargs;
        //myargs.addArgument("filename").parse(args.argc,args.argv);
        //std::string filename = myargs.value("filename");

        update_bbox();


        bbox_drawable = mtao::opengl::make_drawable(bbox_mesh, _flat_shader, drawables());
        bbox_mesh.setParent(&root());
        bbox_drawable->deactivate();
        bbox_drawable->activate_edges();

        bbox_drawable->data().color = 0xffffff_rgbf;


        curve_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>{ curve_mesh, _flat_shader, curve_drawgroup };
        curve_mesh.setParent(&root());
        curve_drawable->deactivate();
        cutcell_drawable = new mtao::opengl::Drawable<Magnum::Shaders::VertexColor2D>{ cutcell_mesh, vcolor_shader, background_drawgroup };
        cutcell_mesh.setParent(&root());
        cutcell_drawable->deactivate();

        cutcell_face_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>{ cutcell_face_mesh, _flat_shader, curve_drawgroup };
        cutcell_face_mesh.setParent(&root());
        cutcell_face_drawable->deactivate();


        grid_mesh.setParent(&root());
        grid_drawable = new mtao::opengl::Drawable<Magnum::Shaders::Flat2D>(grid_mesh, _flat_shader, background_drawgroup);
        grid_drawable->deactivate();

        reset_curve();

        Corrade::Utility::Arguments myargs;
        if (args.argc > 1) {
            myargs.addArgument("filename").parse(args.argc, args.argv);
            std::string filename = myargs.value("filename");
            spdlog::warn("Reading file [{}]", filename);
            curve.load(filename);
            auto p = curve.points();
            {
                using namespace mandoline::construction;
                auto edges = curve.edges().cast<int>().eval();
                FaceCollapser fc(edges);
                fc.bake(p, false);
                auto f = fc.faces_no_holes();
                std::cout << "input facecollapser size: " << f.size() << std::endl;
            }

            bbox.min() = bbox.min().cwiseMin(p.rowwise().minCoeff().cast<float>());
            bbox.max() = bbox.max().cwiseMax(p.rowwise().maxCoeff().cast<float>());
            bbox.min().array() -= .1f;
            bbox.max().array() += .1f;
            std::cout << "Updating curve" << std::endl;
            update_curve();
            std::cout << "Updating bbox" << std::endl;
            update_bbox();
            ui_mode = InterfaceMode::Browse;
        }
    }

    void update_bbox() {
        bbox_mesh.set_bbox(bbox);
    }
    void draw() override;
    void gui() override;
    void update_edges();
    void update_curve();
    void reset_curve();
    void clear_curve();
    void update_faces();
    void update_ccm();
    void update_colors();// doesnt upload the colors to anything

    void mouseMoveEvent(MouseMoveEvent &event) override;
    void mousePressEvent(MouseEvent &event) override;
    void update_face(int idx);

  private:
    Magnum::SceneGraph::DrawableGroup2D background_drawgroup, curve_drawgroup;


    Magnum::Shaders::Flat2D _flat_shader;
    Magnum::Shaders::VertexColor2D vcolor_shader;
    mtao::opengl::objects::Mesh<2> curve_mesh;
    mtao::opengl::objects::Mesh<2> cutcell_mesh;
    mtao::opengl::objects::Mesh<2> cutcell_face_mesh;
    mtao::opengl::objects::Grid<2> grid_mesh;
    mtao::opengl::objects::BoundingBox<2> bbox_mesh;
    mtao::opengl::Drawable<Magnum::Shaders::Flat2D> *bbox_drawable = nullptr;
    mtao::opengl::Drawable<Magnum::Shaders::Flat2D> *grid_drawable = nullptr;
    mtao::opengl::Drawable<Magnum::Shaders::Flat2D> *curve_drawable = nullptr;
    mtao::opengl::Drawable<Magnum::Shaders::VertexColor2D> *cutcell_drawable = nullptr;
    mtao::opengl::Drawable<Magnum::Shaders::Flat2D> *cutcell_face_drawable = nullptr;
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
    if (ui_mode == InterfaceMode::CurveEdit) {
    } else if (ui_mode == InterfaceMode::Browse) {

        int ofi = face_index;
        face_index = ccm->cell_index(mtao::Vec2d(cursor.x(), cursor.y()));
        if (ofi != face_index) {
            update_face(face_index);
        }
    }
}
void MeshViewer::mousePressEvent(MouseEvent &event) {
    Window2::mousePressEvent(event);
    if (!ImGui::GetIO().WantCaptureMouse) {
        if (ui_mode == InterfaceMode::CurveEdit) {
            if (event.button() == MouseEvent::Button::Left) {
                mtao::Vec2d p(cursor.x(), cursor.y());
                curve.add_point(p);
                update_curve();
            }
        } else {
            int idx = ccm->cell_index(mtao::Vec2d(cursor.x(), cursor.y()));
            if (idx >= 0 && idx < ccm->num_cells()) {
                auto c = ccm->cell(idx);
                if (c.size() > 0) {
                    auto d = *c.begin();
                    std::cout << mtao::eigen::stl2eigen(d).transpose() << std::endl;
                }
            }
        }
    }
}

void MeshViewer::clear_curve() {
    curve.clear();
}
void MeshViewer::reset_curve() {

    curve.clear();
    if (!curve.is_closed()) {
        curve.toggle_closed();
    }
    curve.add_point(mtao::Vec2d(-.3, -.3));
    curve.add_point(mtao::Vec2d(-.3, .3));
    curve.add_point(mtao::Vec2d(.3, .3));
    curve.add_point(mtao::Vec2d(.3, -.3));
}

void MeshViewer::gui() {
    {
        const char *items[] = { "CurveEdit", "Browse" };
        int current_item = int(ui_mode);
        if (ImGui::Combo("Interface Mode", &current_item, items, 2)) {
            ui_mode = static_cast<InterfaceMode>(current_item);
        }
    }
    {


        const char *items[] = { "Random", "Regions", "Harmonic", "Harmonid_RHS" };
        int current_item = int(color_type);
        if (ImGui::Combo("Color Type", &current_item, items, 4)) {
            color_type = static_cast<ColorType>(current_item);
            update_faces();
        }
    }


    if (ImGui::InputInt2("N", N.data())) {
        //update_bbox();
    }

    if (ImGui::InputInt("Face index", &face_index)) {
        if (ccm) {
            int nf = ccm->num_faces();
            face_index = std::clamp<int>(face_index, 0, nf - 1);
        }
        update_face(face_index);
    }
    if (ImGui::SliderFloat2("min", bbox.min().data(), -20, 20)) {
        bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(), bbox.max());
        update_bbox();
    }
    if (ImGui::SliderFloat2("max", bbox.max().data(), -20, 20)) {
        bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(), bbox.max());
        update_bbox();
    }

    {
        bool value = curve.is_closed();
        if (ImGui::Checkbox("Closed", &value)) {
            curve.toggle_closed();
            update_curve();
        }
    }
    {
        if (ImGui::Button("Make CCM")) {
            update_ccm();

            update_curve();
        }
        if (ImGui::Button("Clear Curve")) {
            clear_curve();

            update_curve();
        }
        if (ImGui::Button("Reset Curve")) {
            reset_curve();
            update_curve();
        }
    }

    if (ImGui::SliderFloat2("Harmonid cir", harmonic_dir.data(), -1, 1)) {
        update_faces();
    }
    if (ccm) {
        ImGui::Text("Cursor position (%f,%f) cell: %d", cursor.x(), cursor.y(), face_index);
    } else {
        ImGui::Text("Cursor position (%f,%f)", cursor.x(), cursor.y());
    }
}

void MeshViewer::update_edges() {
}
void MeshViewer::update_curve() {

    auto points = curve.points();
    auto V = curve.points().cast<float>().eval();
    auto E = curve.edges();
    if (V.size() == 0) {
        curve_drawable->deactivate();
    } else {
        if (E.size() > 0) {
            curve_mesh.setEdgeBuffer(V, E);
            curve_drawable->activate_edges();
        } else {
            curve_mesh.setVertexBuffer(V);
            curve_drawable->activate_points();
        }
    }
}

void MeshViewer::update_ccm() {
    auto E = curve.edges();
    auto stlp = curve.stl_points();


    auto grid = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox.cast<double>(), N);
    mandoline::construction::CutCellGenerator<2> ccg(stlp, grid);
    ccg.add_boundary_elements(E.cast<int>());
    ccg.bake();

    ccm = ccg.generate();


    //set_colors(ccm.active_grid_cell_mask);

    {
        auto V = ccm->vertices();
        auto E = ccm->cut_edges_eigen();
        //mtao::ColVectors<float, 3> ccg_cols_edge(3, E.cols());
        //ccg_cols_edge.setZero();
        //for (int i = 0; i < ccg_cols_edge.cols(); ++i) {
        //    ccg_cols_edge(i % 3, i) = 1;
        //}
        //mtao::ColVectors<float, 3> ccg_cols(3, V.cols());
        //ccg_cols.setZero();
        //for (int i = 0; i < ccm->StaggeredGrid::vertex_size(); ++i) {
        //    ccg_cols.col(i).setConstant(.7);
        //}
        cutcell_mesh.setEdgeBuffer(V.cast<float>().eval(), E.cast<unsigned int>().eval());
        cutcell_drawable->deactivate();

        spdlog::warn("Making gird draable");
        //ccm->faces();
        grid_mesh.set(ccm->vertex_grid());
        grid_drawable->activate_edges();
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


void MeshViewer::update_face(int idx) {
    if (!ccm) return;

    if (idx < 0 || idx >= ccm->num_cells()) return;
    auto c = ccm->cell(idx);
    if (c.size() == 0) { return; }
    auto d = *c.begin();
    if (d.size() < 3) return;
    auto V = ccm->vertices();
    auto F = mtao::geometry::mesh::earclipping(V, d);
    cutcell_face_mesh.setTriangleBuffer(V.cast<float>(), F.cast<unsigned int>());
    //cutcell_face_mesh.setTriangleBuffer(F.cast<unsigned int>());
    cutcell_face_drawable->activate_triangles();
}

void MeshViewer::update_faces() {
    if (!ccm) return;

    spdlog::warn("updating faces");

    update_colors();
    std::vector<std::tuple<mtao::ColVecs3i, mtao::Vec4d>> FCs;
    FCs.resize(ccm->num_cells());
    spdlog::warn("making meshes");
    auto verts = ccm->vertices();
    //tbb::parallel_for(tbb::blocked_range<size_t>(0, ccm->num_cells()), [&](const tbb::blocked_range<size_t> &range) {
    //    for (size_t idx = range.begin(); idx != range.end(); ++idx) {
    for(int idx = 0; idx < ccm->num_cells(); ++idx) {
        //spdlog::info("Getting cell {}", idx);
            auto c = ccm->cell(idx);
            if (c.size() > 0) {
                auto &&d = *c.begin();
                //std::copy(d.begin(),d.end(),std::ostream_iterator<int>(std::cout,","));
                //std::cout << std::endl;

                //std::cout << idx << " => " << mtao::geometry::curve_volume(verts, d) << std::endl;
                auto F = mtao::geometry::mesh::earclipping(verts, d);
                FCs[idx] = std::tuple<mtao::ColVecs3i, mtao::Vec4d>{ std::move(F), colors.col(idx) };
            }
            //spdlog::info("{} / {}", idx, ccm->num_cells());
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

        tbb::parallel_for(tbb::blocked_range<size_t>(0, R.size()), [&](const tbb::blocked_range<size_t> &range) {
            for(size_t idx = range.begin(); idx != range.end(); ++idx) {
            colors.col(idx) = C.col(R(idx));
            } });
    };
    auto random_colors = [&]() {
        spdlog::warn("computing random colors");
        colors.setRandom();
    };

    auto get_rhs = [&]() {
        spdlog::warn("computing rhs");
        mtao::Vec2d dir = harmonic_dir.cast<double>();
        mtao::VecXd u(ccm->num_edges());
        u.setZero();
        for (auto &&[eidx, edge] : mtao::iterator::enumerate(ccm->cut_edges())) {
            if (edge.external_boundary) {
                auto [oc, b] = *edge.external_boundary;
                if (oc == -2) {
                    u(eidx) = (b ? -1 : 1) * dir(edge.unbound_axis());
                } else {
                    u(eidx) = 0;
                }
            } else {
                if (edge.is_axial_edge()) {
                    u(eidx) = 0;
                } else {// boundary faces are not allowed to emit anything
                    u(eidx) = 0;
                }
            }
        }
        for (auto &&[idx, bfp] : mtao::iterator::enumerate(ccm->exterior_grid.boundary_facet_pairs())) {
            if (ccm->exterior_grid.is_boundary_facet(idx)) {
                int axis = ccm->exterior_grid.boundary_facet_axes()[idx];
                u(idx + ccm->cut_edges().size()) = dir(axis);
            } else {
                u(idx + ccm->cut_edges().size()) = 0;// grid domain boundaries are not allowed to emit anything
            }
        }
        Eigen::SparseMatrix<double> D = mandoline::operators::divergence(*ccm, true);
        mtao::VecXd b = D * u;
        return b;
    };
    auto harmonic_rhs_colors = [&]() {
        spdlog::warn("harmonic rhs colors");
        Eigen::SparseMatrix<double> L = mandoline::operators::laplacian(*ccm);

        mtao::VecXd x = get_rhs();
        x /= x.cwiseAbs().maxCoeff();
        colors.row(0) = x;
        colors.row(2) = -x;
        colors.row(1).array() = 0;
    };
    auto harmonic_colors = [&]() {
        spdlog::warn("harmonic colors");
        Eigen::SparseMatrix<double> L = mandoline::operators::laplacian(*ccm);
        //std::cout << L << std::endl;

        mtao::VecXd b = get_rhs();
        mtao::VecXd x(b.rows());
        x.setZero();
        mtao::solvers::linear::CholeskyPCGSolve(L, b, x);
        x /= x.cwiseAbs().maxCoeff();
        colors.row(0) = x;
        colors.row(2) = -x;
        colors.row(1).array() = 0;
    };
    switch (color_type) {
    case ColorType::Regions:
        region_colors();
        break;
    case ColorType::Random:
        random_colors();
        break;
    case ColorType::Harmonic:
        harmonic_colors();
        break;
    case ColorType::Harmonic_RHS:
        harmonic_rhs_colors();
        break;
    }
    spdlog::warn("Picked out colors");
    colors.noalias() = (colors.array() > 0).select(colors, 0);
    colors.row(3).setConstant(1);
}

MAGNUM_APPLICATION_MAIN(MeshViewer)

