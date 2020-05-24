#include "mandoline/construction/construct_gui_widget.hpp"
#include <mtao/geometry/bounding_box.hpp>
#include <spdlog/spdlog.h>
namespace mandoline::construction {

CutmeshGeneratorGui::CutmeshGeneratorGui(Magnum::Shaders::Flat3D &shader, Magnum::SceneGraph::DrawableGroup3D &group) : mtao::opengl::objects::Grid<3>{}, mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{ *this, shader, group } {
    mtao::opengl::Drawable<Magnum::Shaders::Flat3D>::deactivate();
    mtao::opengl::Drawable<Magnum::Shaders::Flat3D>::activate_edges();
}

CutmeshGeneratorGui::CutmeshGeneratorGui(Magnum::Shaders::Flat3D &shader, Magnum::SceneGraph::DrawableGroup3D &group, const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const mtao::geometry::grid::StaggeredGrid3d &grid, int adaptive_level, std::optional<double> threshold)
  : DeformingGeometryConstructor(V, F, grid, adaptive_level, threshold), mtao::opengl::objects::Grid<3>{ grid.vertex_grid() }, mtao::opengl::Drawable<Magnum::Shaders::Flat3D>{ *this, shader, group }, bbox(grid.bbox()), N(grid.vertex_shape()), use_cube(false), adaptive_level(adaptive_level), threshold(threshold) {
    mtao::opengl::Drawable<Magnum::Shaders::Flat3D>::deactivate();
    mtao::opengl::Drawable<Magnum::Shaders::Flat3D>::activate_edges();
}
CutmeshGeneratorGui CutmeshGeneratorGui::create(Magnum::Shaders::Flat3D &shader, Magnum::SceneGraph::DrawableGroup3D &group, const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, double bbox_scale, const std::array<int, 3> &N, int adaptive_level, std::optional<double> threshold) {
    auto bbox = mtao::geometry::bounding_box(V);
    bbox = mtao::geometry::expand_bbox(bbox, bbox_scale);
    auto s = bbox.sizes().eval();

    for (int i = 0; i < 3; ++i) {
        if (s(i) < 1e-5) {
            bbox.min()(i) -= 1e-2;
            bbox.max()(i) += 1e-2;
        }
    }
    auto g = mtao::geometry::grid::StaggeredGrid3d::from_bbox(bbox, N);
    return CutmeshGeneratorGui(shader, group, V, F, g, adaptive_level, threshold);
}

void CutmeshGeneratorGui::set_bbox(const Eigen::AlignedBox<float, 3> &bbox_, float scale) {
    bbox = mtao::geometry::expand_bbox(bbox_, scale);

    update_grid();
}

bool CutmeshGeneratorGui::gui() {
    bool dirty = false;

    if (ImGui::InputInt3("N", N.data())) {
        dirty = true;
    }
    if (ImGui::InputFloat3("min", bbox.min().data())) {
        bbox.min() = (bbox.min().array() < bbox.max().array()).select(bbox.min(), bbox.max());
        dirty = true;
    }
    if (ImGui::InputFloat3("max", bbox.max().data())) {
        bbox.max() = (bbox.min().array() > bbox.max().array()).select(bbox.min(), bbox.max());
        dirty = true;
    }
    if (ImGui::Checkbox("Cubes", &use_cube)) {
        dirty = true;
    }
    if (ImGui::InputInt("Adaptive level", &adaptive_level)) {

        set_adaptivity(adaptive_level);
        dirty = true;
    }
    if (dirty) {
        update_grid();
    }
    return dirty;
}


void CutmeshGeneratorGui::update_grid(const mtao::geometry::grid::StaggeredGrid3d &grid) {
    DeformingGeometryConstructor::update_grid(grid);
    if (valid()) {
        mtao::opengl::objects::Grid<3>::set(grid.vertex_grid());
    }
}
void CutmeshGeneratorGui::update_grid() {
    int min = *std::min_element(N.begin(), N.end());
    if (min < 2) {
        spdlog::warn("Grid shape is too small ({} {} {}) (min must be >= 2)", N[0], N[1], N[2]);
        return;
    }
    auto sg = staggered_grid();
    update_grid(sg);
}

mtao::geometry::grid::StaggeredGrid3d
  CutmeshGeneratorGui::staggered_grid() const {
    return mtao::geometry::grid::StaggeredGrid3d::from_bbox(bbox.cast<double>(), N, use_cube);
}

mandoline::CutCellMesh<3> CutmeshGeneratorGui::generate() {
    if (valid()) {
        bake();
        return emit();
    } else {
        spdlog::warn("Generated from an invalid input");
        return {};
    }
}
void CutmeshGeneratorGui::update_vertices_and_bbox(const mtao::ColVecs3d &V, double bbox_scale, const std::optional<double> &threshold) {
    update_vertices(V, threshold);

    bbox = mtao::geometry::bounding_box(V).cast<float>();
    bbox = mtao::geometry::expand_bbox(bbox, float(bbox_scale));
    auto s = bbox.sizes().eval();

    for (int i = 0; i < 3; ++i) {
        if (s(i) < 1e-5) {
            bbox.min()(i) -= 1e-2;
            bbox.max()(i) += 1e-2;
        }
    }
    auto g = mtao::geometry::grid::StaggeredGrid3d::from_bbox(bbox.cast<double>(), N);
    update_grid();
}
void CutmeshGeneratorGui::update_mesh_and_bbox(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, double scale, const std::optional<double> &threshold) {

    update_vertices_and_bbox(V, scale, threshold);
    update_topology(F);
}
}// namespace mandoline::construction
