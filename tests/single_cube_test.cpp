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
#include "../tools/make_cutmesh_generator_from_cmdline.hpp"
#include "../tools/make_cutmesh_from_cmdline.hpp"
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


int main(int argc, char * argv[]) {

    //active_loggers["default"].set_level(Level::Error);
    auto ccg = CutCellGenerator<3>(.5 * mtao::Vec3d::Ones(),std::array<int,3>{{2,2,2}});
    ccg.bake();
    auto ccm = ccg.generate();
    
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


    return 0;
}


