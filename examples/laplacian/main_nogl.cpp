#include "mandoline/mesh3.hpp"
#include "setup.h"

int main(int argc, char* argv[]) {
    auto ccm = read(argv[1]);
    mtao::Vec3d dir = mtao::Vec3d::Unit(1);
    std::cout << divergence(ccm,dir).transpose() << std::endl;
    auto B = boundary(ccm);
    Eigen::SparseMatrix<double> B2 = boundary(ccm).transpose();
    std::cout << B.rows() << "," << B.cols() << std::endl;
    for(int i = 0; i < B.cols(); ++i) {
        std::cout << B.col(i).nonZeros()  << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < B2.cols(); ++i) {
        std::cout << B2.col(i).nonZeros()  << " ";
    auto b = B2.col(i);
        if(b.nonZeros() == 1) {
            if(i < ccm.faces().size()) {
                auto& f = ccm.faces()[i];
                std::cout << std::string(f) << std::endl;
                int a = ccm.cell_shape()[f.as_axial_axis()];
                int v = f.as_axial_coord();
                std::cout << "ASDF" << a << ": " << v << std::endl;
            }
            for(int j = 0; j < b.size(); ++j) {
                if(b.coeff(j) != 0) {
                    if(b.coeff(j) < ccm.cells().size()) {
                        auto& c = ccm.cells()[j];
                        std::cout << j << ": " << c.region << std::endl;
                    }
                }
            }
        }
    }
    std::cout << std::endl;
    pressure(ccm,dir);
    //std::cout << std::endl;
    //std::cout << ccm.dual_edge_lengths().transpose() << std::endl << std::endl;
    //std::cout << ccm.face_volumes().transpose() << std::endl << std::endl;
    //std::cout << ccm.cell_volumes().transpose() << std::endl << std::endl;
    return 0;


}
