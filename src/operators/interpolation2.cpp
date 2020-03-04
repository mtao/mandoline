#include "mandoline/operators/interpolation2.hpp"
#include "mandoline/operators/volume2.hpp"
#include <iostream>


namespace mandoline::operators {
using coord_type = CutCellMesh<2>::coord_type;


//mesh vertex -> cut vertex
Eigen::SparseMatrix<double> barycentric_matrix(const CutCellMesh<2> &ccm) {
    Eigen::SparseMatrix<double> A(ccm.vertex_size(), ccm.origV().cols());
    std::vector<Eigen::Triplet<double>> trips;
    std::map<std::array<int, 2>, double> mp;
    for (auto &&[eid, btf] : ccm.mesh_edges()) {
        auto t = btf.sparse_matrix_entries(ccm.cut_edges()[eid], ccm.origE());
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
Eigen::SparseMatrix<double> edge_barycentric_volume_matrix(const CutCellMesh<2> &ccm) {
    int edge_size = 0;
    //artifact from before i passed in m_origF
    if (ccm.origE().size() == 0) {
        for (auto &&f : ccm.cut_edges()) {
            if (f.is_mesh_edge()) {
                edge_size = std::max<int>(edge_size, f.as_edge_id());
            }
        }
        edge_size++;
    } else {
        edge_size = ccm.origE().cols();
    }
    Eigen::SparseMatrix<double> A(ccm.num_edges(), edge_size);
    std::vector<Eigen::Triplet<double>> trips;
    for (auto &&[eid, btf] : ccm.mesh_edges()) {
        double vol = btf.volume();
        trips.emplace_back(eid, btf.parent_eid, vol);
    }
    A.setFromTriplets(trips.begin(), trips.end());
    //for(int i = 0; i < A.cols(); ++i) {
    //    A.col(i) /= A.col(i).sum();
    //}
    return A;
}
//grid vertex -> cut vertex
Eigen::SparseMatrix<double> trilinear_matrix(const CutCellMesh<2> &ccm) {
    Eigen::SparseMatrix<double> A(ccm.vertex_size(), ccm.StaggeredGrid::vertex_size());
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(ccm.StaggeredGrid::vertex_size() + 4 * ccm.cut_vertex_size());

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
                trips.emplace_back(index, ccm.vertex_index(a), vij);
            }
        }
    }
    A.setFromTriplets(trips.begin(), trips.end());
    return A;
}
//grid face -> cut face
Eigen::SparseMatrix<double> edge_grid_volume_matrix(const CutCellMesh<2> &ccm) {
    auto trips = ccm.exterior_grid.boundary_facet_to_staggered_grid(ccm.cut_edges().size());
    auto EV = edge_lengths(ccm);
    double vol = ccm.dx().prod();
    Eigen::SparseMatrix<double> A(ccm.num_edges(), ccm.form_size<1>());

    for (auto &&[row, edge] : mtao::iterator::enumerate(ccm.cut_edges())) {
        //std::cout << "CE: " << row << "/" << ccm.cut_edges().size() << std::endl;
        // extract the lowest coordinate
        if (edge.count() == 1) {
            int axis = edge.unbound_axis();
            auto [a, b] = edge.indices;
            auto A = ccm.masked_vertex(a);
            auto B = ccm.masked_vertex(b);
            coord_type c = edge.get_min_coord([&](int idx) { return ccm.masked_vertex(idx).coord; });
            const int col = ccm.StaggeredGrid::staggered_index<1>(c, axis);
            double value = EV(row) / vol;//face.N(axis) should be a unit vector either facing up or down....
            //std::cout << "Entry: " << row << "," << col << ": " << value<< std::endl;
            trips.emplace_back(row, col, value);
        }
    }
    //std::cout << "Full triplets" << std::endl;
    //for(auto&& t: trips) {
    //    std::cout << t.row() << "," << t.col() << "=" << t.value() << "/" << A.rows() << "," << A.cols() << std::endl;
    //}
    A.setFromTriplets(trips.begin(), trips.end());
    //mtao::VecXd sums = A * mtao::VecXd::Zero(A.cols());
    //sums = (sums.array().abs() > 1e-10).select(1.0 / sums.array(), 0);
    //A = sums.asDiagonal() * A;
    return A;
}
}// namespace mandoline::operators
