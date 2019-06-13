#include <mandoline/adaptive_grid.hpp>


int main(int argc, char * argv[]) {
    using AG = mandoline::AdaptiveGrid;
    using AGF = mandoline::AdaptiveGridFactory;
    using GM = AGF::GridData3;

    using SG = AG::Base;

    int N = 4;
    GM mask(N,N,N);
    int NC = N*N*N;
    mask.set_constant(true);

    AGF agf(mask);
    agf.make_cells(1);
    auto ag = agf.create();

    auto g = ag.grid();
    for(auto&& [i,c]: ag.cells()) {
        std::cout << i << ": " << mtao::eigen::stl2eigen(c.corner()).transpose() << " => " << c.width() << std::endl;
    }
    for(int k = 0; k < 4; ++k) {
        for(int j = 0; j < 4; ++j) {
            for(int i = 0; i < 4; ++i) {
                std::cout << g(i,j,k) << " ";
            }
            std::cout << std::endl;
        }
            std::cout << std::endl;
    }

    Eigen::SparseMatrix<double> B(ag.num_faces(),ag.num_cells());
    auto trips = ag.boundary_triplets(0);
    B.setFromTriplets(trips.begin(),trips.end());
    std::cout << B << std::endl;

    Eigen::SparseMatrix<double> L = B.transpose() * B;
    std::cout << L << std::endl;
    std::cout << mtao::VecXd::Ones(ag.num_cells()).transpose() * L << std::endl;
    return 0;

}
