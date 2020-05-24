#pragma once
#include <imgui.h>

#include "mandoline/construction/construct.hpp"
#include "mtao/opengl/drawables.h"
#include "mtao/opengl/objects/grid.h"
#include <Magnum/Shaders/Flat.h>

namespace mandoline::construction {

//Imgui wrapper around DeformingGeometryConstructor
struct CutmeshGeneratorGui : public DeformingGeometryConstructor, public mtao::opengl::objects::Grid<3>, public mtao::opengl::Drawable<Magnum::Shaders::Flat3D> {
    CutmeshGeneratorGui(Magnum::Shaders::Flat3D& shader, Magnum::SceneGraph::DrawableGroup3D& group);
    CutmeshGeneratorGui(Magnum::Shaders::Flat3D& shader, Magnum::SceneGraph::DrawableGroup3D& group, const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const mtao::geometry::grid::StaggeredGrid3d &grid, int adaptive_level = 0, std::optional<double> threshold = {});
    static CutmeshGeneratorGui create(Magnum::Shaders::Flat3D& shader, Magnum::SceneGraph::DrawableGroup3D& group, const mtao::ColVecs3d &V = {}, const mtao::ColVecs3i &F = {}, double bbox_scale = 1.1, const std::array<int, 3> &N = std::array<int, 3>{ { 5, 5, 5 } }, int adaptive_level = 0, std::optional<double> threshold = {});
    CutmeshGeneratorGui(CutmeshGeneratorGui &&) = default;
    CutmeshGeneratorGui &operator=(CutmeshGeneratorGui &&) = default;
    using DeformingGeometryConstructor::update_grid;
    using DeformingGeometryConstructor::update_vertices;
    using DeformingGeometryConstructor::bake;
    using DeformingGeometryConstructor::emit;

    void update_vertices_and_bbox(const mtao::ColVecs3d &V, double scale = 1.1, const std::optional<double> &threshold = -1);
    void update_mesh_and_bbox(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, double scale = 1.1, const std::optional<double> &threshold = -1);
    //Cutcell mesh parameters. all floats but will be cast to double
    Eigen::AlignedBox<float, 3> bbox;
    std::array<int, 3> N{ { 5, 5, 5 } };
    bool use_cube = false;
    int adaptive_level = 0;
    std::optional<double> threshold = -1;


    bool gui();
    void set_bbox(const Eigen::AlignedBox<float, 3> &bbox, float scale = 1.1);
    void update_grid(const mtao::geometry::grid::StaggeredGrid3d& grid);
    void update_grid();
    mandoline::CutCellMesh<3> generate();
    mtao::geometry::grid::StaggeredGrid3d staggered_grid() const;
};

}// namespace mandoline::construction
