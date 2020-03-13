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

#include <mtao/geometry/mesh/shapes/cube.hpp>

#include <catch2/catch.hpp>
#include <mandoline/construction/generator3.hpp>
using namespace mtao::logging;



using namespace mandoline::construction;

template <typename GridB>
void print_gridb(const GridB& g) {

    if constexpr(GridB::D == 3) 
    {
        for(int i = 0; i < g.shape(0); ++i) {
            for(int j = 0; j < g.shape(1); ++j) {
                for(int k = 0; k < g.shape(2); ++k) {
                    if(g(i,j,k)) {
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
template <typename GridB>
void print_grid(const GridB& g) {

    if constexpr(GridB::D == 3) 
    {
        for(int i = 0; i < g.shape(0); ++i) {
            for(int j = 0; j < g.shape(1); ++j) {
                for(int k = 0; k < g.shape(2); ++k) {
                    std::cout << g(i,j,k) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

}




TEST_CASE("3D Cube", "[ccm3]") {

    auto [V,F] = mtao::geometry::mesh::shapes::cube<double>();

    V.array() += .5;
    //V.colwise() += mtao::Vec3d(.5,.4,.6);
    mtao::geometry::grid::Grid3d vertex_grid(std::array<int,3>{{3,3,3}},mtao::Vec3d::Ones());
    auto ccg = CutCellGenerator<3>(V,vertex_grid);
    ccg.set_boundary_elements(F);

    ccg.bake();

    std::cout << ccg.all_V() << std::endl;

    for(auto&& eisects: ccg.data().edge_intersections()) {
        //std::cout << eisects.edge_index  << ")";
        //for(auto&& e: eisects.vptr_edge) {
        //    std::cout << std::string(*e) << " ";
        //}
        //std::cout << " ===> ";
        REQUIRE(eisects.intersections.size() == 1);
        REQUIRE(eisects.intersections.front().edge_coord == Approx(.5));
        //for(auto&& p: eisects.intersections) {
        //    std::cout << std::string(p) << " ";
        //}
        //std::cout << std::endl;
    }

    for(auto&& fisects: ccg.data().triangle_intersections()) {
        int tri_size = 0;
        int quad_size = 0;
        for(auto&& f: fisects.vptr_faces()) {
            if(f.size() == 3) {
                tri_size++;
                std::array<int,4> counts{{0,0,0,0}};
                for(auto&& p: f) {
                    //std::cout << std::string(*p) << std::endl;
                    counts[p->mask().count()]++;
                }
                REQUIRE(counts[0]==1);
                REQUIRE(counts[1]==1);
                REQUIRE(counts[2]==1);
                REQUIRE(counts[3]==0);
            } else if(f.size() == 4) {
                quad_size++;
                std::array<int,4> counts{{0,0,0,0}};
                for(auto&& p: f) {
                    counts[p->mask().count()]++;
                }
                REQUIRE(counts[0]==1);
                REQUIRE(counts[1]==2);
                REQUIRE(counts[2]==1);
                REQUIRE(counts[3]==0);
            }
        }
        REQUIRE(tri_size == 2);
        REQUIRE(quad_size == 1);
        /*
        std::cout << fisects.triangle_index  << ")";
        for(auto&& e: fisects.vptr_tri) {
            std::cout << std::string(*e) << " ";
        }
        std::cout << " ===> ";
        for(auto&& p: fisects.intersections) {
            std::cout << std::string(p) << " ";
        }
        std::cout << std::endl;
        std::cout << "==============" << std::endl;
        for(auto&& [ve,eisects]: fisects.edge_intersections) {
            for(auto&& e: ve) {
                std::cout << std::string(*e) << " ";
            }
            std::cout << " ===> " << std::endl;
            std::cout << eisects.edge_index  << ")";
            for(auto&& e: eisects.vptr_edge) {
                std::cout << std::string(*e) << " ";
            }
            std::cout << " ===> ";
            for(auto&& p: eisects.intersections) {
                std::cout << std::string(p) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "==============" << std::endl;
        std::cout << "==============" << std::endl;
        */
    }



    auto ccm = ccg.generate();
    
    /*
    std::cout  << ccg.V() << std::endl;
    std::cout  << ccg.origV().size() << std::endl;
    for(auto&& p: ccg.origV()) {
        std::cout << ">" << p.transpose() << std::endl;
    }

    std::cout << "CCG active cells: " << std::endl;
    for(auto&& c: ccg.active_cells()) {
        std::cout << c[0] << "," << c[1] << "," << c[2] << " ";
    }
    std::cout << std::endl;

    print_gridb(ccm.active_cell_mask());
    for(auto&& c: ccm.active_cells()) {
        std::cout << c[0] << "," << c[1] << "," << c[2] << " ";
    }
    std::cout << std::endl;
    std::cout << "Num faces: " << ccm.face_size() << std::endl;
    for(auto&& f: ccm.faces()) {
        std::cout << std::string(f) << std::endl;
    }

    std::cout << "Num cells: " << ccm.cell_size() << std::endl;
    for(auto&& c: ccm.cells()) {
        for(auto&& [f,s]: c) {
            std::cout << "[" << f <<":" << s << "]";

        }
        std::cout << std::endl;
    }
    */

}


