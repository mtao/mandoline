#include "mandoline/construction/construct2_gui_widget.hpp"
#include <spdlog/spdlog.h>
#include <mtao/geometry/bounding_box.hpp>
namespace mandoline::construction {

CutmeshGenerator2Gui::CutmeshGenerator2Gui(Magnum::Shaders::Flat2D &shader, Magnum::SceneGraph::DrawableGroup2D &group) : mtao::opengl::objects::Grid<2>{std::array<int,2>{{2,2}}}, mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>{ *this, shader, group } {
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>::deactivate();
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>::activate_edges();
    update_grid();
}
CutmeshGenerator2Gui::CutmeshGenerator2Gui(Magnum::Shaders::Flat2D &shader, Magnum::SceneGraph::DrawableGroup2D &group, const balsa::eigen::ColVecs2d &V, const balsa::eigen::ColVecs2i &F, const mtao::geometry::grid::StaggeredGrid2d &grid, std::optional<double> threshold)
  : DeformingGeometryConstructor2(V, F, grid, threshold), mtao::opengl::objects::Grid<2>{ grid.vertex_grid() }, mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>{ *this, shader, group }, bbox(grid.bbox()), N(grid.vertex_shape()), use_cube(false), threshold(threshold) {
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>::deactivate();
    mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D>::activate_edges();
    update_grid();
}

CutmeshGenerator2Gui::~CutmeshGenerator2Gui() {}

void CutmeshGenerator2Gui::initialize_visualization() {
    using namespace Magnum::Math::Literals;
    deactivate();
    activate_edges();
    data().color = 0xffffff_rgbf;
}

CutmeshGenerator2Gui CutmeshGenerator2Gui::create(Magnum::Shaders::Flat2D &shader, Magnum::SceneGraph::DrawableGroup2D &group, const balsa::eigen::ColVecs2d &V, const balsa::eigen::ColVecs2i &F, double bbox_scale, const std::array<int, 2> &N, std::optional<double> threshold) {
    auto bbox = mtao::geometry::bounding_box(V);
    bbox = mtao::geometry::expand_bbox(bbox, bbox_scale);
    auto s = bbox.sizes().eval();

    for (int i = 0; i < 2; ++i) {
        if (s(i) < 1e-5) {
            bbox.min()(i) -= 1e-2;
            bbox.max()(i) += 1e-2;
        }
    }
    auto g = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox, N);
    return CutmeshGenerator2Gui(shader, group, V, F, g, threshold);
}

void CutmeshGenerator2Gui::set_bbox(const Eigen::AlignedBox<float, 2> &bbox_, float scale) {
    bbox = mtao::geometry::expand_bbox(bbox_, scale);

    update_grid();
}

bool CutmeshGenerator2Gui::gui() {
    bool dirty = false;

    if (ImGui::InputInt2("N", N.data())) {
        dirty = true;
    }
    if (ImGui::InputFloat2("min", bbox.min().data())) {
        bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(), bbox.max());
        dirty = true;
    }
    if (ImGui::InputFloat2("max", bbox.max().data())) {
        bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(), bbox.max());
        dirty = true;
    }
    if (ImGui::Checkbox("Cubes", &use_cube)) {
        dirty = true;
    }
    if (dirty) {
        update_grid();
    }
    return dirty;
}

void CutmeshGenerator2Gui::update_grid() {
    int min = *std::min_element(N.begin(), N.end());
    if (min < 2) {
        spdlog::warn("Grid shape is too small ({} {}) (min must be >= 2)", N[0], N[1]);
        return;
    }
    if(bbox.sizes().minCoeff() <= 0) {
        spdlog::warn("Grid edge lenghts must be positive: {} {} => {} {} (lengths {} {})",
                bbox.min().x(),
                bbox.min().y(),
                bbox.max().x(),
                bbox.max().y(),
                bbox.sizes().x(),
                bbox.sizes().y());
        return;
    }
    auto sg = staggered_grid();
    update_grid(sg);
}

void CutmeshGenerator2Gui::update_grid(const mtao::geometry::grid::StaggeredGrid2d &grid) {
    DeformingGeometryConstructor2::update_grid(grid);
    if (valid()) {
        mtao::opengl::objects::Grid<2>::set(grid.vertex_grid());
    }
}

mtao::geometry::grid::StaggeredGrid2d
  CutmeshGenerator2Gui::staggered_grid() const {
    return mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox.cast<double>(), N, use_cube);
}

mandoline::CutCellMesh<2> CutmeshGenerator2Gui::generate() {
    if (valid()) {
        bake();
        return emit();
    } else {
        spdlog::warn("Generated from an invalid input");
        return {};
    }
}
void CutmeshGenerator2Gui::update_vertices_and_bbox(const balsa::eigen::ColVecs2d &V, double bbox_scale, const std::optional<double> &threshold) {
    update_vertices(V, threshold);

    bbox = mtao::geometry::bounding_box(V).cast<float>();
    bbox = mtao::geometry::expand_bbox(bbox, float(bbox_scale));
    auto s = bbox.sizes().eval();

    for (int i = 0; i < 2; ++i) {
        if (s(i) < 1e-5) {
            bbox.min()(i) -= 1e-2;
            bbox.max()(i) += 1e-2;
        }
    }
    auto g = mtao::geometry::grid::StaggeredGrid2d::from_bbox(bbox.cast<double>(), N);
    update_grid();
}
void CutmeshGenerator2Gui::update_mesh_and_bbox(const balsa::eigen::ColVecs2d &V, const balsa::eigen::ColVecs2i &F, double scale, const std::optional<double> &threshold) {

    update_vertices_and_bbox(V, scale, threshold);
    update_topology(F);
}
}// namespace mandoline::construction
