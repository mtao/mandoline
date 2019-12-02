#include <mandoline/mesh3.hpp>
#include <mtao/cmdline_parser.hpp>
#include <iostream>

void print_file_info(const std::string& filename)
{
    using namespace mandoline;

    CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(filename);

    std::cout << "Filename: " << filename << std::endl;

    if(ccm.empty())
    {
        std::cout << "Empty cutmesh or not a cutmesh!" << std::endl;
        return;
    }

    {
        auto s = ccm.shape();
        std::cout << "Grid resolution: ";
        std::copy(s.begin(),s.end(),std::ostream_iterator<int>(std::cout," "));
        std::cout << std::endl;
    }
    if(size_t size = ccm.cell_size(); size > 0) {
        std::cout << "Number of cells: ";
        std::cout << ccm.cell_size() << " (";

        if(size_t size = ccm.cut_cell_size(); size > 0) {
            std::cout << size << " cut-cells, ";
        }
        if(size_t size = ccm.adaptive_grid().num_cells(); size > 0) {
            std::cout << size << " cubic/adaptive-cells";
        }
        std::cout << ")" << std::endl;
    } else {
        std::cout << "No cutcells!" << std::endl;
    }
    if(size_t size = ccm.face_size(); size > 0) {
        std::cout << "Number of faces: ";
        std::cout << ccm.face_size() << " (";

        if(size_t size = ccm.cut_face_size(); size > 0) {
            std::cout << size << " cut-faces, ";
        }
        if(size_t size = ccm.adaptive_grid().num_faces(); size > 0) {
            std::cout << size << " cubic/adaptive-faces";
        }
        std::cout << ")" << std::endl;
    } else {
        std::cout << "No cutfaces!" << std::endl;
    }
    std::cout << std::endl;
    {
        auto r = ccm.regions();
        std::map<int,int> region_counts;
        for(auto&& v: r)
        {
            region_counts[v]++;
        }
        std::cout << "Region information: (" << region_counts.size() << " regions found)" << std::endl;
        for(auto&& [id,count]: region_counts)
        {
            std::cout << ">>Region " << id << " has " << count << " cells" << std::endl;
        }
    }
}

int main(int argc, char * argv[]) {
    mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
    if(argc < 1)
    {
        std::cout << "cutmesh_info <filename1> <filename2> ..." << std::endl;
    }
    for(int i = 1; i < argc; ++i)
    {
        if(i > 1)
        {
            std::cout << std::endl << std::endl;
        }
        print_file_info(argv[i]);
    }


    return 0;

}
