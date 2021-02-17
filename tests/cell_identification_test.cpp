#include <igl/WindingNumberAABB.h>

#include <catch2/catch.hpp>
#include <iostream>
#include <iterator>
#include <mandoline/construction/generator3.hpp>
#include <mandoline/mesh3.hpp>
#include <mtao/geometry/mesh/sphere.hpp>
#include <random>

using E = std::array<int, 2>;
using namespace mandoline::construction;

TEST_CASE("Sphere cell access", "[cell_access]") {
    auto sg = mtao::geometry::grid::StaggeredGrid<double, 3>::from_bbox(
        {mtao::Vec3d::Constant(-2), mtao::Vec3d::Constant(2)},
        std::array<int, 3>{{10, 10, 10}});

    auto [V, F] = mtao::geometry::mesh::shapes::sphere<double>(4);
    auto VV = V.transpose().eval();
    auto FF = F.transpose().eval();

    igl::WindingNumberAABB<mtao::Vec3d, decltype(VV), decltype(FF)> aabb(VV,
                                                                         FF);
    aabb.set_mesh(VV, FF);
    aabb.init();

    CutCellGenerator<3> ccg(V, sg);
    ccg.add_boundary_elements(F);
    ccg.bake();
    auto ccm = ccg.generate();

    auto bb = sg.bbox();
    auto sample_bb = bb;
    sample_bb.min().array() -= .2;
    sample_bb.max().array() += .2;
    srand(0);
    auto regions = ccm.regions();

    std::mt19937 gen(1);

    // normal distribution around the radius to increase the hits on non-trivial
    // cutcells
    std::normal_distribution<> d{1, .2};
    mtao::ColVecs3d pts(3, 10000);

    spdlog::info("Making points");
    for (int j = 0; j < pts.cols(); ++j) {
        auto p = pts.col(j);
        p = sample_bb.sample();
        p.normalize();
        double rad = d(gen);
        p *= rad;
    }
    spdlog::info("Computing cells");
    auto I = ccm.get_cell_indices(pts);
    spdlog::info("Running checks");
    for (int j = 0; j < pts.cols(); ++j) {
        auto p = pts.col(j);
        int cell_index = I(j);
        CHECK(ccm.get_cell_index(p) == cell_index);
        if (!bb.contains(p)) {
            CHECK(cell_index == -2);
            continue;
        }
        CHECK(cell_index >= 0);
        CHECK(cell_index < ccm.cell_size());
        CHECK(ccm.is_in_cell(p, cell_index));

        int region = regions[cell_index];
        double w = igl::winding_number(VV, FF, p);
        // spdlog::info("{}/{} {} {}", aabb.max_abs_winding_number(p), w,
        // p.norm(),
        //             region);
        // if (aabb.inside(p)) {
        if (w > .5) {
            CHECK(region == 1);
        } else {
            CHECK(region == 0);
        }
    }
}

