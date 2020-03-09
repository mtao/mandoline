#include <iostream>
#include <mandoline/exterior_grid.hpp>
#include <catch2/catch.hpp>
#include <mandoline/construction/generator2.hpp>
#include <mandoline/operators/diffgeo2.hpp>
#include <mandoline/operators/volume2.hpp>
#include <iterator>
#include <random>

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
        //std::cout << mandoline::operators::edge_lengths(ccm).transpose() << std::endl;
        //std::cout << mandoline::operators::dual_edge_lengths(ccm).transpose() << std::endl;
        //std::cout << mandoline::operators::dual_hodge1(ccm).transpose() << std::endl;
        Eigen::SparseMatrix<double> D = mandoline::operators::divergence(ccm);
        Eigen::SparseMatrix<double> L = mandoline::operators::laplacian(ccm);
        //std::cout << "Divergence: \n"
        //          << D << std::endl;
        //std::cout << "Laplacian: \n"
        //          << L << std::endl;

        Eigen::MatrixXd LL = L;
        Eigen::MatrixXd LSL = LL - LL.transpose();
        Eigen::VectorXd B = L * Eigen::VectorXd::Ones(L.cols());

        REQUIRE((LSL).norm() == Approx(0));
        REQUIRE(B.norm() < 1e-5);

        std::mt19937 gen(0);
        std::uniform_real_distribution<> dis(-10., 10.);


        // some example pressure solves
        for (int k = 0; k < 50; ++k) {
            mtao::VecXd u(ccm.num_edges());
            for (auto &&[eidx, edge] : mtao::iterator::enumerate(ccm.cut_edges())) {
                if (edge.is_axial_edge()) {
                    u(eidx) = dis(gen);
                } else {// boundary faces are not allowed to emit anything
                    u(eidx) = 0;
                }
            }
            for (auto &&[idx, bfp] : mtao::iterator::enumerate(ccm.exterior_grid.boundary_facet_pairs())) {
                if (!ccm.exterior_grid.is_boundary_facet(idx)) {
                    u(idx + ccm.cut_edges().size()) = dis(gen);
                } else {
                    u(idx + ccm.cut_edges().size()) = 0;// grid domain boundaries are not allowed to emit anything
                }
            }

            //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double, Eigen::Upper | Eigen::Lower>> cg;
            //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::IncompleteCholesky<double>> cg;
            {
                mtao::VecXd b = D * u;
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> cg;
                cg.compute(L);
                mtao::VecXd x = cg.solve(b);
                REQUIRE(cg.info() == Eigen::Success);

                mtao::VecXd pg = ccm.boundary(false) * x;
                double err = (D * (u - pg)).norm();
                REQUIRE(err == Approx(0.).margin(1e-6));
            }
            //{
            //    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower> cg;
            //    Eigen::SparseMatrix<double> B = L + mandoline::operators::surface_laplacian(ccm);
            //    cg.compute(B);
            //    mtao::VecXd x = cg.solve(b);
            //    REQUIRE(cg.info() == Eigen::Success);

            //    mtao::VecXd pg = ccm.boundary(false) * x;
            //    double err = (D * (u - pg)).norm();
            //    REQUIRE(err == Approx(0.).margin(1e-6));
            //}
        }
    }
}
//switch (cg.info()) {
//case Eigen::Success:
//    std::cout << "Success" << std::endl;
//    break;
//case Eigen::NumericalIssue:
//    std::cout << "NumericalIssue" << std::endl;
//    break;
//case Eigen::NoConvergence:
//    std::cout << "NoConvergence" << std::endl;
//    break;
//case Eigen::InvalidInput:
//    std::cout << "InvalidInput" << std::endl;
//    break;
//}
