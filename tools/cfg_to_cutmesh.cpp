#include <Corrade/Utility/Arguments.h>
#include <Corrade/Utility/Configuration.h>
#include <cstring>
#include <iterator>
#include <cstdlib>
#include <mtao/iterator/enumerate.hpp>
#include <mtao/logging/logger.hpp>
#include <mtao/logging/profiler.hpp>
#include <set>
#include <mtao/geometry/prune_vertices.hpp>
#include <mandoline/construction/remesh_self_intersections.hpp>
#include <mandoline/construction/construct.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
#include <igl/read_triangle_mesh.h>
#pragma GCC diagnostic pop

std::string find_file(const std::string& filename, const std::string& path) {
    //TODO: actually use path
    return path+"/"+filename;
}


mandoline::CutCellMesh<3> parse_configuration(const std::string& cfg_file, const std::string& path_file = ".", bool no_cleaning=false) {
    Corrade::Utility::Configuration cfg{cfg_file};

    std::string input_mesh = find_file(cfg.value("input_mesh"),path_file);

    auto parse_3toks = [](const std::string& str) -> std::array<std::string,3> {
        std::array<std::string,3> r;
        size_t start = -1;
        for(auto&& [i, v]: mtao::iterator::enumerate(r)) {
            start = str.find_first_not_of(' ',start+1);
            size_t end = str.find_first_of(' ',start);
            size_t size = (end == std::string::npos?str.size():end) - start;

            std::cout << str << " => " << start << "," << size << std::endl;
            v = str.substr(start,size);
            start = end;
        }
        return r;
    };
    auto parse_vec3d = [&](const std::string& str) -> mtao::Vec3d {
        auto vals = parse_3toks(str);
        mtao::Vec3d r;
        std::transform(vals.begin(),vals.end(),r.data(),[](const std::string& tok) -> double {
                return std::strtod(tok.data(),nullptr);
                });
        return r;
    };
    auto parse_vec3i = [&](const std::string& str) -> std::array<int,3> {
        auto vals = parse_3toks(str);
        std::array<int,3> r;
        std::transform(vals.begin(),vals.end(),r.begin(),[](const std::string& tok) -> int {
                return std::atoi(tok.data());
                });
        return r;
    };
    
    Eigen::AlignedBox<double,3> bbox;
    bbox.min() = parse_vec3d(cfg.value("min"));
    bbox.max() = parse_vec3d(cfg.value("max"));
    std::array<int,3> cell_shape = parse_vec3i(cfg.value("shape"));

    std::cout << "File: " << input_mesh << std::endl;
    std::cout << "BBox: " << bbox.min().transpose() << " => " << bbox.max().transpose() << std::endl;
    std::cout << "Cell shape: ";
    std::copy(cell_shape.begin(),cell_shape.end(),std::ostream_iterator<int>(std::cout," "));
    std::cout << std::endl;
    mtao::ColVecs3d V;
    mtao::ColVecs3i F;
    {
        Eigen::MatrixXd VV;
        Eigen::MatrixXi FF;
        igl::read_triangle_mesh(input_mesh,VV,FF);
        V = VV.transpose();
        F = FF.transpose();
        if(!no_cleaning) {
            auto t = mtao::logging::profiler("remesh_time",false,"remesh_profiler");
            std::tie(V,F) = mtao::geometry::prune(V,F,0);
            std::tie(V,F) = mandoline::construction::remesh_self_intersections(V,F);
        }
    }
    return mandoline::construction::from_bbox(V,F,bbox,cell_shape);
}


int main(int argc, char * argv[]) {
    Corrade::Utility::Arguments args;
    args.addArgument("filename").setHelp("filename","Path to cfg file")
        .addArgument("prefix").setHelp("prefix", "Output .cutmesh file prefix")
        .addOption('P',"path",".").setHelp("path","A semi-colon separated list of paths for objs")
        .addBooleanOption('v',"verbose").setHelp("verbose","log verbosity")
        .addBooleanOption("no-cleaning").setHelp("no-cleaning","do not clean mesh of self-intersections and redundant facets")
        .parse(argc,argv);


    if(!args.isSet("verbose")) {
        mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
        mtao::logging::make_logger("profiler").set_level(mtao::logging::Level::Off);
    } else {
        auto&& log = make_logger("profiler",mtao::logging::Level::All);
    }

    auto ccm = parse_configuration(args.value("filename"),args.value("path"), args.isSet("no-cleaning"));

    ccm.write(args.value("prefix"));


    return 0;
}


