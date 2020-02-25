#pragma once
#include <imgui.h>

#include "mandoline/construction/construct.hpp"

namespace mandoline::construction {

//Imgui wrapper around DeformingGeometryConstructor
struct CutmeshGenerator_Imgui : public DeformingGeometryConstructor {
    CutmeshGenerator_Imgui(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, const mtao::geometry::grid::StaggeredGrid3d &grid, int adaptive_level = 0, std::optional<double> threshold = -1);
    static CutmeshGenerator_Imgui create(const mtao::ColVecs3d &V, const mtao::ColVecs3i &F, double bbox_scale = 1.1, const std::array<int, 3> &N = std::array<int, 3>{ { 5, 5, 5 } }, int adaptive_level = 0, std::optional<double> threshold = -1);
    CutmeshGenerator_Imgui(const CutmeshGenerator_Imgui &) = default;
    CutmeshGenerator_Imgui(CutmeshGenerator_Imgui &&) = default;
    CutmeshGenerator_Imgui &operator=(CutmeshGenerator_Imgui &&) = default;
    CutmeshGenerator_Imgui &operator=(const CutmeshGenerator_Imgui &) = default;
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
    void update_grid();
    mandoline::CutCellMesh<3> generate();
    mtao::geometry::grid::StaggeredGrid3d staggered_grid() const;
};

}// namespace mandoline::construction
