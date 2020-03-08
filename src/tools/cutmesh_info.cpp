#include "mandoline/tools/cutmesh_info.hpp"
namespace mandoline::tools {
    void print_file_info(const std::string& filename)
    {

        CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(filename);
        if(ccm.empty())
        {
            std::cout << "Empty cutmesh or not a cutmesh!" << std::endl;
            return;
        } else {
            std::cout << "Filename: " << filename << std::endl;
            print_all_info(ccm);
        }
    }

    void print_all_info(const CutCellMesh<3>& ccm) {
            print_general_info(ccm);
            print_cell_info(ccm);
            print_face_info(ccm);
            print_region_info(ccm);
    }


    void print_general_info(const CutCellMesh<3>& ccm)
    {

        std::cout << "General Cutcell Information" << std::endl;
        std::cout << "===========================" << std::endl;


        {
            auto s = ccm.cell_shape();
            std::cout << "Grid cell shape: ";
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
    }

    void print_face_info(const CutCellMesh<3>& ccm)
    {

        std::cout << "Face Information" << std::endl;
        std::cout << "================" << std::endl;
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
    }
    void print_region_info(const CutCellMesh<3>& ccm)
    {
        std::cout << "Region Information" << std::endl;
        std::cout << "==================" << std::endl;
        auto r = ccm.regions();
        std::map<int,int> region_counts;
        for(auto&& v: r)
        {
            region_counts[v]++;
        }
        std::cout << region_counts.size() << " regions found" << std::endl;
        for(auto&& [id,count]: region_counts)
        {
            std::cout << ">>Region " << id << " has " << count << " cells" << std::endl;
        }
    }

    void print_cell_info(const CutCellMesh<3>& ccm) {
        std::cout << "Cell Information" << std::endl;
        std::cout << "================" << std::endl;
        for(auto&& c: ccm.cells()) {
            auto g = c.grid_cell;
            for(auto&& [dim,i,m]: mtao::iterator::enumerate(g,ccm.cell_shape())) {
                if(i < 0 || i >= m) {
                    std::cout << "Cell outside grid (dim=" << dim << ")" << i << "/" << m << std::endl;
                    std::cout << "Face/sign pairs: ";
                    for(auto&& [fidx,b]: c) {
                        std::cout << fidx << ":" << b << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}
