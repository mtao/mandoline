#pragma once
#include <imgui.h>

#include "mandoline/construction/construct2.hpp"
#include "mtao/opengl/drawables.h"
#include "mtao/opengl/objects/grid.h"
#include <Magnum/Shaders/Flat.h>

namespace mandoline::construction {



//Imgui wrapper around DeformingGeometryConstructor
struct CutmeshGenerator2Gui : public DeformingGeometryConstructor2, public mtao::opengl::objects::Grid<2>, public mtao::opengl::MeshDrawable<Magnum::Shaders::Flat2D> {
    CutmeshGenerator2Gui(Magnum::Shaders::Flat2D& shader, Magnum::SceneGraph::DrawableGroup2D& group);
    CutmeshGenerator2Gui(Magnum::Shaders::Flat2D& shader, Magnum::SceneGraph::DrawableGroup2D& group,const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, const mtao::geometry::grid::StaggeredGrid2d &grid, std::optional<double> threshold = {});
    static CutmeshGenerator2Gui create(Magnum::Shaders::Flat2D& shader, Magnum::SceneGraph::DrawableGroup2D& group,const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, double bbox_scale = 1.1, const std::array<int, 2> &N = std::array<int, 2>{ { 5, 5 } }, std::optional<double> threshold = {});
    CutmeshGenerator2Gui(CutmeshGenerator2Gui &&) = default;
    CutmeshGenerator2Gui &operator=(CutmeshGenerator2Gui &&) = default;
    ~CutmeshGenerator2Gui();
    using DeformingGeometryConstructor2::update_grid;
    using DeformingGeometryConstructor2::update_vertices;
    using DeformingGeometryConstructor2::bake;
    using DeformingGeometryConstructor2::emit;

    void update_vertices_and_bbox(const mtao::ColVecs2d &V, double scale = 1.1, const std::optional<double> &threshold = -1);
    void update_mesh_and_bbox(const mtao::ColVecs2d &V, const mtao::ColVecs2i &F, double scale = 1.1, const std::optional<double> &threshold = -1);
    //Cutcell mesh parameters. all floats but will be cast to double
    Eigen::AlignedBox<float, 2> bbox = Eigen::AlignedBox<float,2>(Eigen::Vector2f::Constant(-1.), Eigen::Vector2f::Constant(1.));;
    std::array<int, 2> N{ { 5, 5 } };
    bool use_cube = false;
    std::optional<double> threshold = -1;


    // If dirty is emitted the user MUST pass vertices back in
    bool gui();
    void draw();
    void initialize_visualization();
    void set_bbox(const Eigen::AlignedBox<float, 2> &bbox, float scale = 1.1);
    void update_grid();
    void update_grid(const mtao::geometry::grid::StaggeredGrid2d &g);
    mandoline::CutCellMesh<2> generate();
    mtao::geometry::grid::StaggeredGrid2d staggered_grid() const;
    private:

    //void update_bbox_visualization();
    void update_grid_visualization();




};

}// namespace mandoline::construction
