#include <iostream>
#include <mandoline/exterior_grid.hpp>
#include <catch2/catch.hpp>
#include <mandoline/construction/generator2.hpp>
#include <iterator>

using E = std::array<int, 2>;
using namespace mandoline::construction;

TEST_CASE("2D", "[boundary,exterior_grid]") {
    mtao::ColVecs2d V(2, 4);

    mtao::ColVecs2i E(2, 4);

    auto make_test = [&]() {
        auto sg = mtao::geometry::grid::StaggeredGrid<double, 2>::from_bbox({ mtao::Vec2d::Zero(), mtao::Vec2d::Constant(4) }, std::array<int, 2>{ { 5, 5 } });
        CutCellGenerator<2> ccg(V, sg);
        ccg.add_boundary_elements(E);
        ccg.bake();
        auto ccm = ccg.generate();
        /*
    std::cout << "Faces: " << std::endl;
    for(auto&& [i,f]: mtao::iterator::enumerate(ccm.m_faces)) {
        std::cout << i << " ";
        for(auto&& c: f.indices) {
            std::copy(c.begin(),c.end(),std::ostream_iterator<int>(std::cout,","));
            std::cout << " ";
        }
        std::cout << std::endl;
        
    }
    std::cout << "Edges: " << std::endl;
    for(auto&& [i,e]: mtao::iterator::enumerate(ccm.cut_edges())) {
        std::cout << i << " ";
        std::copy(e.indices.begin(),e.indices.end(),std::ostream_iterator<int>(std::cout,","));
        if(e.external_boundary) {
            auto [bc,s] = *e.external_boundary;
            std::cout << "[" << bc << "(" << s << ")]";
        }
        
        std::cout << std::endl;
        
    }
    */


        {
            auto B = ccm.boundary(true);
            auto Bd = ccm.boundary(false);
            auto br = Eigen::MatrixXd(B).cwiseAbs().rowwise().sum().eval();
            auto bdr = Eigen::MatrixXd(Bd).cwiseAbs().rowwise().sum().eval();
            auto bc = Eigen::MatrixXd(B).cwiseAbs().colwise().sum().eval();
            for (int i = 0; i < bc.cols(); ++i) {
                size_t facets = 0;
                auto f = ccm.cell(i);
                for (auto &&c : f) {
                    facets += c.size();
                }
                REQUIRE(facets == bc(i));
            }
            for (int i = 0; i < br.rows(); ++i) {
                auto a = br(i);
                auto b = bdr(i);
                REQUIRE((a == 0 || a == 1 || a == 2));
                REQUIRE(((b == 0) || (b == 2)));
                REQUIRE((!(b == 0) || (a == 1)));
            }
            REQUIRE(Bd * mtao::VecXd::Ones(Bd.cols()) == mtao::VecXd::Zero(Bd.rows()));
        }
    };

    V.col(0) = mtao::Vec2d(1.5, 1.5);
    V.col(1) = mtao::Vec2d(1.5, 2.5);
    V.col(2) = mtao::Vec2d(2.5, 2.5);
    V.col(3) = mtao::Vec2d(2.5, 1.5);

    E.col(0) = mtao::Vec2i(0, 1);
    E.col(1) = mtao::Vec2i(1, 2);
    E.col(2) = mtao::Vec2i(2, 3);
    E.col(3) = mtao::Vec2i(3, 0);
    make_test();
    E.col(0) = mtao::Vec2i(1, 0);
    E.col(1) = mtao::Vec2i(2, 1);
    E.col(2) = mtao::Vec2i(3, 2);
    E.col(3) = mtao::Vec2i(0, 3);
    make_test();
}
