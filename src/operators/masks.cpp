#include "mandoline/operators/masks.hpp"
namespace mandoline::operators {
balsa::eigen::VecXd mesh_face_mask(const CutCellMesh<3> &ccm) {
    balsa::eigen::VecXd M = balsa::eigen::VecXd::Ones(ccm.face_size());
    for (auto &&[i, f] : mtao::iterator::enumerate(ccm.faces())) {
        if (f.is_mesh_face()) {
            M(i) = 0;
        }
    }
    return M;
}
balsa::eigen::VecXd grid_boundary_mask(const CutCellMesh<3> &ccm) {
    balsa::eigen::VecXd M = balsa::eigen::VecXd::Ones(ccm.face_size());
    for (auto &&[i, f] : mtao::iterator::enumerate(ccm.faces())) {
        if (f.is_axial_face()) {
            int cidx = f.as_axial_coord();
            if (cidx == 0 || ccm.vertex_shape()[f.as_axial_axis()] == cidx) {
                M(i) = 0;
            }
        }
    }
    M.tail(ccm.exterior_grid().num_faces()) = ccm.exterior_grid().grid_boundary_face_mask();
    return M;
}

}// namespace mandoline::operators
