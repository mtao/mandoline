#include <iostream>
#include <mandoline/exterior_grid.hpp>
#include <catch2/catch.hpp>
#include <mandoline/construction/generator2.hpp>
#include <mandoline/operators/diffgeo2.hpp>
#include <mandoline/operators/volume2.hpp>
#include <iterator>

using E = std::array<int, 2>;
using namespace mandoline::construction;

TEST_CASE("2D", "[boundary,exterior_grid]") {
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
        std::cout << mandoline::operators::edge_lengths(ccm).transpose() << std::endl;
        std::cout << mandoline::operators::dual_edge_lengths(ccm).transpose() << std::endl;
        std::cout << mandoline::operators::dual_hodge1(ccm).transpose() << std::endl;
        auto D = mandoline::operators::divergence(ccm);
        auto L = mandoline::operators::laplacian(ccm);
        std::cout << "Divergence: \n"
                  << D << std::endl;
        std::cout << "Laplacian: \n"
                  << L << std::endl;

        Eigen::MatrixXd LL = L;
        Eigen::MatrixXd LSL = LL - LL.transpose();
        Eigen::VectorXd B = L * Eigen::VectorXd::Ones(L.cols());

        std::cout << B.transpose() << std::endl;
        REQUIRE((LSL).norm() == Approx(0));
        REQUIRE(B.norm() == Approx(0));
    }
}
