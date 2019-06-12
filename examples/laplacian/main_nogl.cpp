#include "mandoline/mesh3.hpp"
#include "setup.h"

int main(int argc, char* argv[]) {
    auto ccm = read(argv[1]);
    std::cout << divergence(ccm).transpose() << std::endl;
    auto B = boundary(ccm);
    Eigen::SparseMatrix<double> B2 = boundary(ccm).transpose();
    std::cout << B.rows() << "," << B.cols() << std::endl;
    for(int i = 0; i < B.cols(); ++i) {
        std::cout << B.col(i).nonZeros()  << " ";
    }
    std::cout << std::endl;
    for(int i = 0; i < B2.cols(); ++i) {
        std::cout << B2.col(i).nonZeros()  << " ";
    }
    std::cout << std::endl;
    return 0;


}
