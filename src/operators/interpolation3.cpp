#include "mandoline/operators/interpolation3.hpp"


namespace mandoline::operators {
using coord_type = CutCellMesh<3>::coord_type;


//mesh vertex -> cut vertex
Eigen::SparseMatrix<double> barycentric_matrix(const CutCellMesh<3> &ccm) {
    Eigen::SparseMatrix<double> A(ccm.vertex_size(), ccm.origV().cols());
    std::vector<Eigen::Triplet<double>> trips;
    std::map<std::array<int, 2>, double> mp;
    for (auto &&[fid, btf] : ccm.mesh_cut_faces()) {
        auto t = btf.sparse_matrix_entries(ccm.cut_faces()[fid], ccm.origF());
        std::copy(t.begin(), t.end(), std::inserter(mp, mp.end()));
    }
    trips.reserve(mp.size());
    for (auto &&[pr, v] : mp) {
        auto [a, b] = pr;
        trips.emplace_back(a, b, v);
    }
    A.setFromTriplets(trips.begin(), trips.end());
    mtao::VecXd sums = A * mtao::VecXd::Ones(A.cols());
    sums = (sums.array().abs() > 1e-10).select(1.0 / sums.array(), 0);
    A = sums.asDiagonal() * A;
    return A;
}

//mesh face -> cut face
Eigen::SparseMatrix<double> face_barycentric_volume_matrix(const CutCellMesh<3> &ccm) {
    int face_size = 0;
    //artifact from before i passed in m_origF
    if (ccm.origF().size() == 0) {
        for (auto &&f : ccm.cut_faces()) {
            if (f.is_mesh_face()) {
                face_size = std::max<int>(face_size, f.as_face_id());
            }
        }
        face_size++;
    } else {
        face_size = ccm.origF().cols();
    }
    Eigen::SparseMatrix<double> A(ccm.face_size(), face_size);
    std::vector<Eigen::Triplet<double>> trips;
    for (auto &&[fid, btf] : ccm.mesh_cut_faces()) {
        double vol = btf.volume() * 2;//proportion of 2 is required because barycentric coordinates live in a unit triangle
        trips.emplace_back(fid, btf.parent_fid, vol);
    }
    A.setFromTriplets(trips.begin(), trips.end());
    //for(int i = 0; i < A.cols(); ++i) {
    //    A.col(i) /= A.col(i).sum();
    //}
    return A;
}
//grid vertex -> cut vertex
Eigen::SparseMatrix<double> trilinear_matrix(const CutCellMesh<3> &ccm) {
    Eigen::SparseMatrix<double> A(ccm.vertex_size(), ccm.StaggeredGrid::vertex_size());
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(ccm.StaggeredGrid::vertex_size() + 8 * ccm.cut_vertex_size());

    //emplace the identity map for grid vertices
    for (int i = 0; i < ccm.StaggeredGrid::vertex_size(); ++i) {
        trips.emplace_back(i, i, 1);
    }

    int offset = ccm.StaggeredGrid::vertex_size();
    for (auto &&[idx, v] : mtao::iterator::enumerate(ccm.cut_vertices())) {
        int index = idx + offset;
        coord_type a = v.coord;
        for (int i = 0; i < 2; ++i) {
            a[0] = v.coord[0] + i;
            double vi = i == 0 ? (1 - v.quot(0)) : v.quot(0);
            for (int j = 0; j < 2; ++j) {
                a[1] = v.coord[1] + j;
                double vij = vi * (j == 0 ? (1 - v.quot(1)) : v.quot(1));
                for (int k = 0; k < 2; ++k) {
                    a[2] = v.coord[2] + k;
                    double vijk = vij * (k == 0 ? (1 - v.quot(2)) : v.quot(2));
                    trips.emplace_back(index, ccm.vertex_index(a), vijk);
                }
            }
        }
    }
    A.setFromTriplets(trips.begin(), trips.end());
    return A;
}
//grid face -> cut face
Eigen::SparseMatrix<double> face_grid_volume_matrix(const CutCellMesh<3> &ccm, bool include_mesh_cutfaces) {
    auto trips = ccm.adaptive_grid().grid_face_projection(ccm.cut_faces().size());
    auto FV = ccm.face_volumes();
    auto &dx = ccm.dx();
    mtao::Vec3d gfv;
    gfv(0) = dx(1) * dx(2);
    gfv(1) = dx(0) * dx(2);
    gfv(2) = dx(0) * dx(1);
    Eigen::SparseMatrix<double> A(ccm.face_size(), ccm.form_size<2>());
    for (auto &&t : trips) {
        const int row = t.row();
        const int col = t.col();
    }
    for (auto &&[i, face] : mtao::iterator::enumerate(ccm.cut_faces())) {
        // we can't assume is_axial_face() returns something reasoanble because
        // some triangles may lie on a cut-face
        bool included = include_mesh_cutfaces?(face.count() == 1) : face.is_axial_face();
        if (included) {
            int axis = face.bound_axis();
            constexpr static int maxval = std::numeric_limits<int>::max();
            coord_type c{ { maxval, maxval, maxval } };
            for (auto &&ind : face.indices) {
                for (auto &&i : ind) {
                    auto v = ccm.masked_vertex(i).coord;
                    for (auto &&[a, b] : mtao::iterator::zip(c, v)) {
                        a = std::min(a, b);
                    }
                }
            }
            const int row = i;
            const int col = ccm.staggered_index<2>(c, axis);
            double value = FV(i) / gfv(axis) * face.N(axis);//face.N(axis) should be a unit vector either facing up or down....
            trips.emplace_back(row, col, value);
        }
    }
    A.setFromTriplets(trips.begin(), trips.end());
    //mtao::VecXd sums = A * mtao::VecXd::Zero(A.cols());
    //sums = (sums.array().abs() > 1e-10).select(1.0 / sums.array(), 0);
    //A = sums.asDiagonal() * A;
    return A;
}
}// namespace mandoline::operators
