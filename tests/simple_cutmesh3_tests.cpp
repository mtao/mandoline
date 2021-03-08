#define __gmp_const const
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gmpq.h>
#include <igl/copyleft/cgal/extract_cells.h>
#include <mtao/types.hpp>

#include <mtao/geometry/mesh/triangle_fan.hpp>
#include <iostream>
#include <mtao/eigen/stl2eigen.hpp>
#include <mtao/eigen/stack.h>
#include <memory>
#include <algorithm>

#include <mtao/eigen_utils.h>
#include <mtao/logging/logger.hpp>
#include <mtao/geometry/grid/triangulation.hpp>
#include <mtao/geometry/bounding_box.hpp>
#include <mtao/geometry/mesh/sphere.hpp>
#include <mtao/geometry/bounding_box.hpp>
#include <igl/read_triangle_mesh.h>
#include <mtao/geometry/mesh/boundary_facets.h>
#include <mtao/geometry/prune_vertices.hpp>
#include "debug_cutface.hpp"

#include <mtao/geometry/mesh/shapes/cube.hpp>

#include <catch2/catch.hpp>
#include <mandoline/construction/generator3.hpp>
#include <mandoline/construction/tools/preprocess_mesh.hpp>
using namespace mtao::logging;


using namespace mandoline::construction;

template<typename GridB>
void print_gridb(const GridB &g) {

    if constexpr (GridB::D == 3) {
        for (int i = 0; i < g.shape(0); ++i) {
            for (int j = 0; j < g.shape(1); ++j) {
                for (int k = 0; k < g.shape(2); ++k) {
                    if (g(i, j, k)) {
                        std::cout << "o";
                    } else {
                        std::cout << ".";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}
template<typename GridB>
void print_grid(const GridB &g) {

    if constexpr (GridB::D == 3) {
        for (int i = 0; i < g.shape(0); ++i) {
            for (int j = 0; j < g.shape(1); ++j) {
                for (int k = 0; k < g.shape(2); ++k) {
                    std::cout << g(i, j, k) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}


TEST_CASE("3D Cube", "[ccm3]") {

    auto [V, F] = mtao::geometry::mesh::shapes::cube<double>();

    V.array() += .5;
    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 3, 3, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();

    std::cout << ccg.all_V() << std::endl;

    for (auto &&eisects : ccg.data().edge_intersections()) {
        REQUIRE(eisects.intersections.size() == 1);
        REQUIRE(eisects.intersections.front().edge_coord == Approx(.5));
    }

    for (auto &&fisects : ccg.data().triangle_intersections()) {
        int tri_size = 0;
        int quad_size = 0;
        for (auto &&f : fisects.vptr_faces()) {
            if (f.size() == 3) {
                tri_size++;
                std::array<int, 4> counts{ { 0, 0, 0, 0 } };
                for (auto &&p : f) {
                    //std::cout << std::string(*p) << std::endl;
                    counts[p->mask().count()]++;
                }
                REQUIRE(counts[0] == 1);
                REQUIRE(counts[1] == 1);
                REQUIRE(counts[2] == 1);
                REQUIRE(counts[3] == 0);
            } else if (f.size() == 4) {
                quad_size++;
                std::array<int, 4> counts{ { 0, 0, 0, 0 } };
                for (auto &&p : f) {
                    counts[p->mask().count()]++;
                }
                REQUIRE(counts[0] == 1);
                REQUIRE(counts[1] == 2);
                REQUIRE(counts[2] == 1);
                REQUIRE(counts[3] == 0);
            }
        }
        REQUIRE(tri_size == 2);
        REQUIRE(quad_size == 1);
    }


    auto ccm = ccg.generate();


    {
        auto V = ccm.vertices();
        //for(auto&& f: ccm.faces()) {
        //    std::cout << std::string(f) << std::endl;
        //}
        int mesh_faces = 0;
        int axial_faces = 0;
        for (auto &&[idx, f] : mtao::iterator::enumerate(ccm.faces())) {

            if (f.is_mesh_face()) {
                mesh_faces++;
            } else {
                axial_faces++;
            }
            if (!(f.is_axial_face() && f.indices.size() == 1 && f.indices.begin()->size() == 6)) {
                CHECK(check_convex_face_normal(V, f));
            }
        }
        REQUIRE(mesh_faces == 36);
        REQUIRE(axial_faces == 48);
    }
    //std::cout << "Vertices: \n";
    //for(int i = 0; i < ccm.vertices().cols(); ++i) {
    //    std::cout << i << ")) " << ccm.vertices().col(i).transpose() << std::endl;
    //}
    ////std::cout  << ccg.V() << std::endl;
    //std::cout  << ccg.origV().size() << std::endl;
    //for(auto&& p: ccg.origV()) {
    //    std::cout << ">" << p.transpose() << std::endl;
    //}

    ///*
    //std::cout << "CCG active cells: " << std::endl;
    //for(auto&& c: ccg.active_cells()) {
    //    std::cout << c[0] << "," << c[1] << "," << c[2] << " ";
    //}
    //std::cout << std::endl;

    //print_gridb(ccm.active_cell_mask());
    //for(auto&& c: ccm.active_cells()) {
    //    std::cout << c[0] << "," << c[1] << "," << c[2] << " ";
    //}
    //std::cout << std::endl;
    //*/

    //std::cout << "Num cells: " << ccm.cell_size() << std::endl;
    //for(auto&& c: ccm.cells()) {
    //    for(auto&& [f,s]: c) {
    //        std::cout << "[" << f <<":" << s << "]";

    //    }
    //    std::cout << std::endl;
    //}
}

TEST_CASE("3D Tet", "[ccm3]") {

    mtao::ColVecs3d V(3, 4);
    mtao::ColVecs3i F(3, 4);
    V.col(0) << 1, 1, 1.5;
    V.col(1) << .5, .5, .5;
    V.col(2) << 1.5, 1.5, .5;
    V.col(3) << .5, 1.5, .5;
    ;

    F.col(0) << 0, 1, 2;
    F.col(1) << 0, 1, 3;
    F.col(2) << 0, 2, 3;
    F.col(3) << 1, 2, 3;


    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 3, 3, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();


    auto ccm = ccg.generate();

    auto R = ccm.regions();
    int regions = *std::max_element(R.begin(), R.end()) + 1;
    REQUIRE(regions == 2);
}

TEST_CASE("3D 2Tets", "[ccm3]") {

    mtao::ColVecs3d V(3, 5);
    mtao::ColVecs3i F(3, 7);
    V.col(0) << 1, 1, 1.5;
    V.col(1) << .5, .5, 1;
    V.col(2) << 1.5, 1.5, 1;
    V.col(3) << .5, 1.5, 1;
    ;
    V.col(4) << 1, 1, .5;
    ;

    F.col(0) << 0, 1, 2;
    F.col(1) << 0, 1, 3;
    F.col(2) << 0, 2, 3;
    F.col(3) << 1, 2, 3;
    F.col(4) << 4, 1, 2;
    F.col(5) << 4, 1, 3;
    F.col(6) << 4, 2, 3;


    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 3, 3, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();


    auto ccm = ccg.generate();

    {
        for (auto &&[idx, f] : mtao::iterator::enumerate(ccm.faces())) {

            if (f.mask()[2] && *f.mask()[2] == 1) {
                std::cout << std::string(f) << std::endl;
                if (f.is_mesh_face()) {
                    //mesh_faces++;
                } else {
                    //axial_faces++;
                }
            }
        }
    }

    auto R = ccm.regions();

    int regions = *std::max_element(R.begin(), R.end()) + 1;
    REQUIRE(regions == 3);
}

TEST_CASE("3D IntersectingTets", "[ccm3]") {

    mtao::ColVecs3d V(3, 4);
    mtao::ColVecs3i F(3, 4);
    V.col(0) << 0, 0, 1;
    V.col(1) << .5, .5, 0;
    V.col(2) << 1.5, 1.5, 0;
    V.col(3) << .5, 1.5, 0;
    ;

    F.col(0) << 0, 1, 2;
    F.col(1) << 0, 1, 3;
    F.col(2) << 0, 2, 3;
    F.col(3) << 1, 2, 3;
    //std::tie(V,F) = mtao::geometry::mesh::shapes::cube<double>();

    F = mtao::eigen::hstack(F, F.array() + V.cols());
    V = mtao::eigen::hstack(V.array() + .25, V.colwise() + mtao::Vec3d(.25, .25, .1));
    std::tie(V, F) = mandoline::construction::tools::preprocess_mesh(V, F);
    //for (int i = 0; i < V.cols(); ++i) {
    //    std::cout << "v " << V.col(i).transpose() << std::endl;
    //}

    //for (int i = 0; i < F.cols(); ++i) {
    //    std::cout << "f " << (F.col(i).array() + 1).transpose() << std::endl;
    //}

    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 3, 3, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();


    auto ccm = ccg.generate();

    auto R = ccm.regions();

    int regions = *std::max_element(R.begin(), R.end()) + 1;
    REQUIRE(regions == 4);
}
TEST_CASE("3D IntersectingCubes", "[ccm3]") {

    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    std::tie(V, F) = mtao::geometry::mesh::shapes::cube<double>();
    std::cout << V << std::endl;

    F = mtao::eigen::hstack(F, F.array() + V.cols());
    V = mtao::eigen::hstack(V.colwise() + mtao::Vec3d(.5, .5, .3), V.colwise() + mtao::Vec3d(.5, .5, .7));
    std::tie(V, F) = mandoline::construction::tools::preprocess_mesh(V, F);
    std::cout << V << std::endl;

    //std::cout << V << std::endl;
    //std::cout << F << std::endl;

    //for (int i = 0; i < V.cols(); ++i) {
    //    std::cout << "v " << V.col(i).transpose() << std::endl;
    //}

    //for (int i = 0; i < F.cols(); ++i) {
    //    std::cout << "f " << (F.col(i).array() + 1).transpose() << std::endl;
    //}

    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 3, 3, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();


    auto ccm = ccg.generate();

    auto R = ccm.regions();

    int regions = *std::max_element(R.begin(), R.end()) + 1;


    {
        Eigen::MatrixXi C;
        Eigen::MatrixXd VV = V.transpose();
        Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Eigen::Dynamic, Eigen::Dynamic> VVV = VV.cast<CGAL::Lazy_exact_nt<CGAL::Gmpq>>();
        Eigen::MatrixXi FF = F.transpose();
        igl::copyleft::cgal::extract_cells(VV, FF, C);
        std::set<int> regions;
        std::copy(C.data(), C.data() + C.size(), std::inserter(regions, regions.end()));

        spdlog::info("Input mesh had {} regions", regions.size());
    }


    REQUIRE(regions == 4);
}
TEST_CASE("3D TwoBoxes", "[ccm3]") {

    mtao::ColVecs3d V(3, 12);
    V.col(0) << .2, .2, .2;
    V.col(1) << .8, .2, .2;
    V.col(2) << .8, .8, .2;
    V.col(3) << .2, .8, .2;
    V.col(4) << .2, .2, 1;
    V.col(5) << .8, .2, 1;
    V.col(6) << .8, .8, 1;
    V.col(7) << .2, .8, 1;
    V.col(8) << .2, .2, 1.8;
    V.col(9) << .8, .2, 1.8;
    V.col(10) << .8, .8, 1.8;
    V.col(11) << .2, .8, 1.8;

    mtao::ColVecs3i F(3, 6 + 16);

    F.col(0) << 0, 1, 2;
    F.col(1) << 0, 2, 3;

    F.col(2) << F.col(0).array() + 4;
    F.col(3) << F.col(1).array() + 4;
    F.col(4) << F.col(0).array() + 8;
    F.col(5) << F.col(1).array() + 8;

    F.col(6) << 0, 1, 5;
    F.col(7) << 0, 5, 4;

    F.col(8) << 1, 2, 6;
    F.col(9) << 1, 6, 5;

    F.col(10) << 2, 3, 7;
    F.col(11) << 2, 7, 6;

    F.col(12) << 3, 0, 4;
    F.col(13) << 3, 4, 7;

    F.col(14) << F.col(6).array() + 4;
    F.col(15) << F.col(7).array() + 4;

    F.col(16) << F.col(8).array() + 4;
    F.col(17) << F.col(9).array() + 4;

    F.col(18) << F.col(10).array() + 4;
    F.col(19) << F.col(11).array() + 4;

    F.col(20) << F.col(12).array() + 4;
    F.col(21) << F.col(13).array() + 4;


    //std::cout << V << std::endl;
    //std::cout << F << std::endl;

    for (int i = 0; i < V.cols(); ++i) {
        std::cout << "v " << V.col(i).transpose() << std::endl;
    }

    for (int i = 0; i < F.cols(); ++i) {
        std::cout << "f " << (F.col(i).array() + 1).transpose() << std::endl;
    }

    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int, 3>{ { 2, 2, 3 } }, mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V, vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();


    auto ccm = ccg.generate();
    for (auto &&[fidx,f] : mtao::iterator::enumerate(ccm.faces())) {
        std::cout << fidx << "))) "<< std::string(f.mask()) << "  =  " << std::string(f)  << ": " << f.N.transpose() << std::endl;
        auto m = f.mask();
        if(m[2] && *m[2] == 1) {
            for(auto&& s: f.indices) {
                for(auto&& v: s) {
                    std::cout << std::string(ccm.masked_vertex(v)) << " ";
                }
                std::cout << " ";
            }
                std::cout << std::endl;
        }
    }
    for(auto&& [cidx,c]: mtao::iterator::enumerate(ccm.cells())) {
        std::cout << cidx << ")" << std::string(c) << std::endl;
    }

    auto R = ccm.regions();

    int regions = *std::max_element(R.begin(), R.end()) + 1;


    {
        Eigen::MatrixXi C;
        Eigen::MatrixXd VV = V.transpose();
        Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Eigen::Dynamic, Eigen::Dynamic> VVV = VV.cast<CGAL::Lazy_exact_nt<CGAL::Gmpq>>();
        Eigen::MatrixXi FF = F.transpose();
        igl::copyleft::cgal::extract_cells(VV, FF, C);
        std::set<int> regions;
        std::copy(C.data(), C.data() + C.size(), std::inserter(regions, regions.end()));

        spdlog::info("Input mesh had {} regions", regions.size());
    }


    REQUIRE(regions == 3);
}
