#include "make_cutmesh_from_cmdline.hpp"


int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);

    auto clp = make_cutmesh_clparser();

    clp.add_options()
    ("write-regions","write per-region objs",cxxopts::value<bool>()->default_value("false"))
    ("write-monolithic","write a single obj with all info",cxxopts::value<bool>()->default_value("false"))
    ("write-separate","write per-cell objs",cxxopts::value<bool>()->default_value("false"))
    ("write-flaps","write non-manifold components in objs",cxxopts::value<bool>()->default_value("false"))
    ("write-mesh","write input mesh as obj",cxxopts::value<bool>()->default_value("false"));

    auto result = clp.parse(argc,argv);

        bool write_regions = result["write-regions"].as<bool>(); 
        bool write_monolithic = result["write-monolithic"].as<bool>();
        bool write_separate = result["write-separate"].as<bool>();
        bool write_flaps = result["write-flaps"].as<bool>();
        bool write_mesh = result["write-mesh"].as<bool>();
        std::string output_prefix = result["mesh_file"].as<std::string>();
        auto ccm = make_cutmesh(result);

        ccm.triangulate_faces();
        if(write_regions   ) {  
            ccm.write_obj_regions(output_prefix);
        }
        if(write_monolithic) {  
            ccm.write_obj(output_prefix);
        }
        if(write_separate  ) { 
            ccm.write_obj_separate(output_prefix,false);
        }
        if(write_flaps ) { 
            ccm.write_obj_flaps(output_prefix);
        }
        if(write_mesh ) { 
            ccm.write_mesh_surface_obj(output_prefix);
            ccm.write_mesh_obj_separate(output_prefix);
        }

    return 0;
}


