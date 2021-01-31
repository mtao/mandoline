#include <spdlog/spdlog.h>

#include "mandoline/mesh3.hpp"
#include "mandoline/operators/interpolation3.hpp"
using namespace mandoline;

bool faces_fully_utilized(const CutCellMesh<3>& ccm) {
    auto B = ccm.face_barycentric_volume_matrix();
    bool good = true;
    if (Eigen::Map<const Eigen::VectorXd>(B.valuePtr(), B.nonZeros())
            .minCoeff() < 0) {
        std::cout << "Inverted face" << std::endl;
        good = false;
    }
    auto BV = B.transpose() * mtao::VecXd::Ones(B.rows());
    for (int i = 0; i < BV.rows(); ++i) {
        if (std::abs(BV(i) - 1) > 1e-5) {
            std::cout << "Underutilized face: " << i << std::endl;
            good = false;
        }
    }
    auto G = mandoline::operators::face_grid_volume_matrix(ccm, true);
    auto GV = G.transpose() * mtao::VecXd::Ones(G.rows());
    // std::cout << G << std::endl;
    // std::cout << "Grid recovered areas" << std::endl;
    // std::cout << GV.transpose() << std::endl;
    for (int i = 0; i < GV.rows(); ++i) {
        if (std::abs(GV(i) - 1) > 1e-5) {
            auto [c, dim] = ccm.form_unindex<2>(i);

            fmt::print("Underutilized grid face {} ({}-face {})\n", i, dim,
                       fmt::join(c, ","));
            good = false;
        }
    }
    return good;
}
