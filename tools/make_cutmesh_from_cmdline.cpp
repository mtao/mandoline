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


using namespace mandoline;



mtao::CommandLineParser make_cutmesh_clparser() {
    mtao::CommandLineParser clp;
    clp.add_option("adaptive",true);
    clp.add_option("normalize");
    clp.add_option("checks");
    clp.add_option("prescaled");
    clp.add_option("N",-1);
    clp.add_option("NI",3);
    clp.add_option("NJ",3);
    clp.add_option("NK",3);
    clp.add_option("adaptive-level",-1);
    clp.add_option("remove-self-intersections",true);
    return clp;
}

CutCellMesh<3> make_cutmesh(int argc, char * argv[]) {
    auto clp = make_cutmesh_clparser();
    if(clp.parse(argc,argv)) {
        return make_cutmesh(clp);
    } else {
        return {};
    }
}
std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> read_mesh_input(const mtao::CommandLineParser& clp) {

    std::string obj_filename = clp.arg(0);
    if(clp.args().size() < 1) {
        fatal() << "No input mesh filename!";
        return {};
    }
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    {
        Eigen::MatrixXd VV;
        Eigen::MatrixXi FF;
        igl::read_triangle_mesh(obj_filename,VV,FF);
        V = VV.transpose();
        F = FF.transpose();
        bool si = clp.optT<bool>("remove-self-intersections");
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
construction::CutCellGenerator<3> make_generator(const mtao::CommandLineParser& clp) {
    auto [V,F] = read_mesh_input(clp);
    return make_generator(V,F,clp);

}
construction::CutCellGenerator<3> make_generator(const mtao::ColVecs3d& VV, const mtao::ColVecs3i& FF, const mtao::CommandLineParser& clp) {
    auto&& log = make_logger("profiler",mtao::logging::Level::All);
    std::string obj_filename = clp.arg(0);
    int adaptive_level = clp.optT<int>("adaptive-level");
    bool adaptive = clp.optT<bool>("adaptive") || adaptive_level >= 0;

    int NI = clp.optT<int>("NI");
    int NJ = clp.optT<int>("NJ");
    int NK = clp.optT<int>("NK");
    int N = clp.optT<int>("N");
    bool do_checks = clp.optT<bool>("checks");
    bool prescaled = clp.optT<bool>("prescaled");
    bool normalize = clp.optT<bool>("normalize");
    std::cout << "N: " << N << std::endl;
    if(N > 0) {
        NI = NJ = NK = N;
    }

    assert(NI>0);


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

    std::array<int,3> NN{{NI,NJ,NK}};
    mtao::Vec3d dx = 1.0 / (mtao::eigen::stl2eigen(NN).cast<double>().array()+1);
    mtao::Vec3d origin = mtao::Vec3d::Zero();

    construction::CutCellGenerator<3> ccg;
    if(prescaled) {
        auto sg = CutCellMesh<3>::StaggeredGrid(NN,mtao::Vec3d::Ones());
        ccg = construction::CutCellGenerator<3>(stlp,sg, {});
        auto bbox = mtao::geometry::bounding_box(V);
    } else {
        auto bbox = mtao::geometry::bounding_box(V);

        double bbox_offset = .1;

        mtao::Vec3d C = (bbox.min() + bbox.max()) / 2;
        mtao::Vec3d s = bbox.sizes() / 2;
        bbox.min() = C - (1 + bbox_offset) * s;
        bbox.max() = C + (1 + bbox_offset) * s;

        auto sg = CutCellMesh<3>::StaggeredGrid::from_bbox(bbox,NN,true);
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
CutCellMesh<3> make_cutmesh(const construction::CutCellGenerator<3>& ccg, const mtao::CommandLineParser& clp) {

    CutCellMesh<3> ccm;
    {
        auto t = mtao::logging::profiler("ccm_generation",false,"profiler");
        ccm = ccg.generate();
    }
    bool normalize = clp.optT<bool>("normalize");
    //ccm.normalize_output = normalize;
    return ccm;
}

CutCellMesh<3> make_cutmesh(const mtao::CommandLineParser& clp) {

    auto t = mtao::logging::timer("Total construction time");
    return make_cutmesh(make_generator(clp),clp);
}
