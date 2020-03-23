#include <iostream>
#include <mandoline/construction/face_collapser.hpp>
#include <catch2/catch.hpp>
#include <mandoline/construction/cutdata.hpp>
#include <mtao/geometry/mesh/face_normals.hpp>
#include <mandoline/construction/generator3.hpp>

using namespace mtao::geometry::trigonometry;

TEST_CASE("Polygon angles", "[trigonometry]") {
    mtao::ColVecs2d V(2, 5);
    V.col(0) = mtao::Vec2d(0.0, 0.0);
    V.col(1) = mtao::Vec2d(1.0, 0.0);
    V.col(2) = mtao::Vec2d(1.0, 1.0);
    V.col(3) = mtao::Vec2d(0.0, 1.0);
    V.col(4) = mtao::Vec2d(0.5, 0.5);


    for (int n = 3; n < 5; ++n) {
        std::stringstream ss;
        ss << "n = " << n;
        SECTION(ss.str().c_str()) {
            std::vector<int> a(n);
            std::iota(a.begin(), a.end(), 0);
            std::vector<int> b(n);
            std::copy(a.rbegin(), a.rend(), b.begin());
            REQUIRE(Approx(M_PI * 2 * n) == interior_angle_sum(n) + exterior_angle_sum(n));
            REQUIRE(angle_sum(V, a) == Approx(interior_angle_sum(n)));
            REQUIRE(angle_sum(V, b) == Approx(exterior_angle_sum(n)));
        }
    }
}
using namespace mandoline::construction;
TEST_CASE("Cutdata Orientation", "[cutdata]") {


    mtao::ColVecs3d V(3, 3);
    V.col(0) << .8, .5, 1;
    V.col(1) << .2, .5, .5;
    V.col(2) << .8, .5, 1.5;
    mtao::ColVecs3i F(3, 1);
    F.col(0) << 0, 1, 2;


    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 2, 2, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);

    ccg.set_boundary_elements(F);
    ccg.bake();
    std::cout << ccg.origN << std::endl;
    auto ccm = ccg.generate();

    V = ccm.vertices();
    for (auto &&f : ccm.faces()) {
        if (f.is_mesh_face()) {
            auto &&inds = *f.indices.begin();
            double N = mtao::geometry::mesh::face_normals(V, mtao::eigen::stl2eigen(inds)).col(0).dot(f.N);
            REQUIRE(N > 0);
            ;
        }
    }
}
