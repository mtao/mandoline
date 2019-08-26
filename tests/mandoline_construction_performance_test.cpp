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






int main(int argc, char * argv[]) {

    active_loggers["default"].set_level(Level::Error);
    auto&& log = make_file_logger("profiler","profiler.log",mtao::logging::Level::Fatal,true);
    auto&& log_remesh = make_file_logger("remesh_profiler","remesh_profiler.log",mtao::logging::Level::Fatal,true);
    auto clp = make_cutmesh_clparser();
    clp.add_option("output-cutmesh",false);
    if(clp.parse(argc,argv)) {
        std::set<int> NIs;
        bool output_file = clp.optT<bool>("output-cutmesh");
        std::optional<std::string> output_prefix;
        int size_offset = 1;
        if(output_file) {
            output_prefix = clp.arg(1);
            size_offset++;
        }

        for(int i = size_offset; i < clp.args().size(); ++i) {
            std::string filename = clp.arg(i);
            NIs.insert(std::atoi(clp.arg(i).c_str()));
        }
        std::string obj_filename = clp.arg(0);

        int NI = clp.optT<int>("NI");

        auto [V,F] = read_mesh_input(clp);



        auto run =[&](int M) {
            clp.set_option("N",M);

            mtao::logging::profiler::log_all();
            auto ccg = make_generator(V,F,clp);
            auto ccm = make_cutmesh(ccg,clp);
            auto&& dur = mtao::logging::profiler::durations();
            for(auto&& [pr,times]: dur) {
                auto&& [name,level] = pr;
                static const std::string pname = "profiler"; if(name == pname) {
                    static const std::string generator_bake = "generator_bake";
                    static const std::string ccm_generation = "ccm_generation";
                    static const std::string datastructure_initialization= "creating facets";
                    static const std::string mesh_bake = "mesh bake";
                    static const std::string grid_cells= "grid bake cells";
                    static const std::string grid_faces= "grid bake faces";
                    static const std::string grid_edges= "grid bake edges";
                    static const std::string grid_vertices = "grid bake vertices";
                    static const std::string mesh_faces= "mesh face bake";
                    static const std::string vertices = "mesh vertex bake";

                    auto get_time = [&](const std::string& str) -> int {
                        auto dms = std::chrono::duration_cast<std::chrono::milliseconds>(times.at(str).first);
                        auto time = dms.count();
                        return time;
                    };
                    log.write(mtao::logging::Level::Fatal,  obj_filename , ",", M, ", " , ccg.data().cut_vertex_size() , "," , ccg.data().edge_intersections().size() , "," , ccg.data().triangle_intersections().size(), ", ",
                            ccm.faces().size(), ",",
                            ccm.cells().size(), ",",
                            ccm.active_cell_count(), ", ",
                            get_time(generator_bake), ",",
                            get_time(ccm_generation), ",",
                            get_time(mesh_bake), ", ",

                            get_time(vertices), ",",
                            get_time(mesh_faces), ", ",

                            get_time(datastructure_initialization), ",",
                            get_time(grid_vertices), ",",
                            get_time(grid_edges), ",",
                            get_time(grid_faces), ",",
                            get_time(grid_cells), ",");


                }
            }
            mtao::logging::profiler::clear();

            if(output_prefix) {
                ccm.write(*output_prefix);
            }
        };
        for(auto&& ni: NIs) {
            run(ni);
        }

    }

    return 0;
}


