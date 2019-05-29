#include "make_cutmesh_from_cmdline.hpp"


int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);

    auto clp = make_cutmesh_clparser();
    clp.add_option("write-regions",false);
    clp.add_option("write-monolithic",false);
    clp.add_option("write-separate",false );
    clp.add_option("write-flaps",false);
    clp.add_option("write-mesh",false);

    if(clp.parse(argc, argv)) {

        bool write_regions = clp.optT<bool>("write-regions");
        bool write_monolithic = clp.optT<bool>("write-monolithic");
        bool write_separate = clp.optT<bool>("write-separate");
        bool write_flaps = clp.optT<bool>("write-flaps");
        bool write_mesh = clp.optT<bool>("write-mesh");
        std::string output_prefix = clp.arg(1);
        auto ccm = make_cutmesh(clp);

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
    }

    return 0;
}


