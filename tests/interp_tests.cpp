#include <iostream>
#include <mandoline/exterior_grid.hpp>
#include <catch2/catch.hpp>
#include <mandoline/construction/generator2.hpp>
#include <mandoline/operators/interpolation2.hpp>
#include <iterator>
#include <random>

using E = std::array<int, 2>;
using namespace mandoline::construction;

TEST_CASE("2D", "[interpolation]") {
    mtao::ColVecs2d V(2, 4);

    mtao::ColVecs2i E(2, 4);

    V.col(0) = mtao::Vec2d(1.5, 1.5);
    V.col(1) = mtao::Vec2d(1.5, 2.5);
    V.col(2) = mtao::Vec2d(2.5, 2.5);
    V.col(3) = mtao::Vec2d(2.5, 1.5);

    E.col(0) = mtao::Vec2i(0, 1);
    E.col(1) = mtao::Vec2i(1, 2);
    E.col(2) = mtao::Vec2i(2, 3);
    E.col(3) = mtao::Vec2i(3, 0);

    auto sg = mtao::geometry::grid::StaggeredGrid<double, 2>::from_bbox({ mtao::Vec2d::Zero(), mtao::Vec2d::Constant(4) }, std::array<int, 2>{ { 5, 5 } });
    CutCellGenerator<2> ccg(V, sg);
    ccg.add_boundary_elements(E);
    ccg.bake();
    auto ccm = ccg.generate();


    {
        Eigen::SparseMatrix<double> E = mandoline::operators::edge_grid_volume_matrix(ccm);
        REQUIRE(E.rows() == ccm.num_edges());
        REQUIRE(E.cols() == ccm.Base::form_size<1>());
        auto CS = Eigen::MatrixXd(E).colwise().sum();
        double min = CS.minCoeff();
        double max = CS.maxCoeff();
        REQUIRE(min == Approx(1));
        REQUIRE(max == Approx(1));
    }
    {
        Eigen::SparseMatrix<double> E = mandoline::operators::face_grid_volume_matrix(ccm);
        std::cout << E << std::endl;
        REQUIRE(E.rows() == ccm.num_faces());
        REQUIRE(E.cols() == ccm.Base::form_size<2>());
        auto CS = Eigen::MatrixXd(E).colwise().sum();
        double min = CS.minCoeff();
        double max = CS.maxCoeff();
        REQUIRE(min == Approx(1));
        REQUIRE(max == Approx(1));
    }
}
