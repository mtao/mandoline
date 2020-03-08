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
#include "mandoline/construction/generator.hpp"
#include <mtao/geometry/bounding_box.hpp>
#include <mtao/geometry/mesh/sphere.hpp>
#include <mtao/geometry/bounding_box.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#include <igl/read_triangle_mesh.h>
#pragma GCC diagnostic pop
#include <mtao/geometry/mesh/boundary_facets.h>
#include <mtao/geometry/prune_vertices.hpp>
#include "mandoline/construction/remesh_self_intersections.hpp"
using namespace mtao::logging;
#include "make_cutmesh_from_cmdline.hpp"
#include "make_cutmesh_generator_from_cmdline.hpp"
#include "mandoline/tools/cutmesh_info.hpp"


using namespace mandoline;



cxxopts::Options make_cutmesh_clparser() {
    cxxopts::Options options("make_cutmesh", "make a cutcell-mesh using mandoline");

    options.add_options()
        ("mesh_file", "mesh grid file",cxxopts::value<std::string>())
        ("output", "output cutmesh file",cxxopts::value<std::string>())
        ("N,shape", "output shape as a triplet of csv NI,NJ,NV" ,cxxopts::value<std::string>()->default_value("5,5,5"))
        ("p,prescaled", "Whether the mesh was already scaled to grid index space",cxxopts::value<bool>()->default_value("false"))
        ("r,rsi", "remove self intersections (may be slow)",cxxopts::value<bool>()->default_value("false"))
        ("a,adaptivity_level", "Number of grid resolutions",cxxopts::value<int>()->default_value("0"))
        ("c,checks", "Do some quality checks on the resulting ccm",cxxopts::value<bool>()->default_value("false"))
        ("n,normalize", "Normalize data to a unit cube",cxxopts::value<bool>()->default_value("false"))
        ("i,info", "show extra info after creating the ccm",cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");
    options.parse_positional({"mesh_file","output"});
    options.positional_help({"<mesh_file> <output_filename>"});
    return options;
}

CutCellMesh<3> make_cutmesh(int argc, char * argv[]) {
    auto clp = make_cutmesh_clparser();
    auto result = clp.parse(argc,argv);
    return make_cutmesh(result);
}
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> read_mesh_input(const cxxopts::ParseResult& result) {

    if(!bool(result.count("mesh_file"))) {
        fatal() << "No input mesh filename!";
        return {};
    }
    std::string obj_filename = result["mesh_file"].as<std::string>();
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    {
        Eigen::MatrixXd VV;
        Eigen::MatrixXi FF;
        igl::read_triangle_mesh(obj_filename,VV,FF);
        V = VV.transpose();
        F = FF.transpose();
        bool si = result["rsi"].as<bool>();
        if(si) {
            auto t = mtao::logging::profiler("remesh_time",false,"remesh_profiler");
            std::tie(V,F) = mtao::geometry::prune(V,F,0);
            std::tie(V,F) = construction::remesh_self_intersections(V,F);
        }
        auto&& dur = mtao::logging::profiler::durations();
        for(auto&& [pr,times_]: dur) {
            auto& times = times_;
            auto&& [name,level] = pr;
            static const std::string pname = "remesh_profiler"; if(name == pname) {
                static const std::string rmt = "remesh_time";
                auto get_time = [&](const std::string& str) -> int {
                    auto dms = std::chrono::duration_cast<std::chrono::milliseconds>(times.at(str).first);
                    auto time = dms.count();
                    return time;
                };
                auto logger = mtao::logging::get_logger("remesh_time", mtao::logging::Level::Fatal);

                logger << obj_filename << ", "  << get_time(rmt);
            }
        }
    }
    return {V,F};
}
construction::CutCellGenerator<3> make_generator(const cxxopts::ParseResult& result) {
    auto [V,F] = read_mesh_input(result);
    return make_generator(V,F,result);

}
construction::CutCellGenerator<3> make_generator(const mtao::ColVecs3d& VV, const mtao::ColVecs3i& FF, const cxxopts::ParseResult& result) {
    auto&& log = make_logger("profiler",mtao::logging::Level::All);

    if(!bool(result.count("mesh_file"))) {
        fatal() << "No input mesh filename!";
        return {};
    }
    std::string obj_filename = result["mesh_file"].as<std::string>();
    int adaptive_level = result["adaptivity_level"].as<int>();
    bool adaptive = adaptive_level >= 0;

    std::string Ns = result["shape"].as<std::string>();
    std::array<int,3> N;
    {
        auto it = Ns.begin();
        auto it2 = Ns.begin();
        size_t pos = 0;
        size_t comma_count = std::count_if(Ns.begin(),Ns.end(), [](std::string::value_type c) -> bool {
                return c == ',';
                });
        if(comma_count != 2) {
            fatal() << "Shape should have 2 commas";
            return {};
        }
        auto step = [&]() {
            it2 = std::find(it,Ns.end(),',');
            std::string ret(it,it2);
            it = it2+1;
            return ret;

        };
        N[0] = std::stoi(step());
        N[1] = std::stoi(step());
        N[2] = std::stoi(step());
    }
    bool do_checks = result["checks"].as<bool>();
    bool prescaled = result["prescaled"].as<bool>();
    bool normalize = result["normalize"].as<bool>();



    mtao::ColVecs3d V = VV;
    mtao::ColVecs3i F = FF;


    //auto [V,F] = mtao::geometry::mesh::read_objD(obj_filename);
    if(!prescaled) {
        auto bbox = mtao::geometry::bounding_box(V);
        V.colwise() -= bbox.center();
    }
    if(normalize) {
        auto bbox = mtao::geometry::bounding_box(V);
        double maxsize = bbox.sizes().maxCoeff();
        V.colwise() -= bbox.center();
        V.array() /= maxsize;
    }
    mtao::vector<mtao::Vec3d> stlp(V.cols());
    for(auto&& [i,v]: mtao::iterator::enumerate(stlp)) {
        v = V.col(i);
    }

    mtao::Vec3d dx = 1.0 / (mtao::eigen::stl2eigen(N).cast<double>().array()+1);
    mtao::Vec3d origin = mtao::Vec3d::Zero();

    construction::CutCellGenerator<3> ccg;
    if(prescaled) {
        auto vg = mtao::geometry::grid::Grid3d(N,mtao::Vec3d::Ones());
        auto sg = CutCellMesh<3>::StaggeredGrid(vg);
        ccg = construction::CutCellGenerator<3>(stlp,sg, {});
        auto bbox = mtao::geometry::bounding_box(V);
    } else {
        auto bbox = mtao::geometry::bounding_box(V);

        double bbox_offset = .1;

        mtao::Vec3d C = (bbox.min() + bbox.max()) / 2;
        mtao::Vec3d s = bbox.sizes() / 2;
        bbox.min() -= bbox_offset * s;
        bbox.max() += bbox_offset * s;


        auto sg = CutCellMesh<3>::StaggeredGrid::from_bbox(bbox,N,true);
        ccg = construction::CutCellGenerator<3>(stlp,sg, {});
    }
    if(adaptive) {
        ccg.adaptive = adaptive;
        if(adaptive_level >= 0) {
            ccg.adaptive_level = adaptive_level;
        }
    }
    {
        auto t = mtao::logging::profiler("generator_bake",false,"profiler");
        ccg.add_boundary_elements(F);
        ccg.bake();
    }
    if(do_checks) {
        if(!ccg.check_cell_containment()) {
            std::cout << "Some faces are too large!" << std::endl;
        }
        if(!ccg.check_face_utilization()) {
            std::cout << "NOT ALL FACES USED!" << std::endl;
        }
    }
    return ccg;
}
CutCellMesh<3> make_cutmesh(const construction::CutCellGenerator<3>& ccg, const cxxopts::ParseResult& result) {

    CutCellMesh<3> ccm;
    {
        auto t = mtao::logging::profiler("ccm_generation",false,"profiler");
        ccm = ccg.generate();
    }
    bool normalize = result["normalize"].as<bool>();

    bool info = result["info"].as<bool>();
    if(info) {
        tools::print_all_info(ccm);
    }
    //ccm.normalize_output = normalize;
    return ccm;
}

CutCellMesh<3> make_cutmesh(const cxxopts::ParseResult& result) {

    auto t = mtao::logging::timer("Total construction time");
    return make_cutmesh(make_generator(result),result);
}
