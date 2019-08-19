#include "mandoline/mesh3.hpp"
using namespace mandoline;

bool faces_fully_utilized(const CutCellMesh<3>& ccm) {
    auto B = ccm.face_barycentric_volume_matrix();
    auto BV = B.transpose() * mtao::VecXd::Ones(B.rows());
    for(int i = 0; i < BV.rows(); ++i) {
        if(ccm.folded_faces().find(i) != ccm.folded_faces().end()) {
            if(std::abs(BV(i) - 1) > 1e-5) {
                return false;
            }
        }
    }
    return true;
}
