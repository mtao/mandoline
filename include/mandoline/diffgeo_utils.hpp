#pragma once
#include <Eigen/Sparse>
#include <array>
#include <numeric>

namespace mandoline {
    Eigen::SparseMatrix<double> grid_boundary(int ni, int nj, bool dirichlet_boundary=true);
    void fix_axis_aligned_indefiniteness(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& rhs);
    Eigen::VectorXd solve_SPSD_system(const Eigen::SparseMatrix<double>& E, const Eigen::VectorXd& rhs);
    Eigen::VectorXd grid_edgeratio(int ni, int nj, const Eigen::VectorXd& phi);


}
