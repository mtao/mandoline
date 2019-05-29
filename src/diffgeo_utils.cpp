#include "diffgeo_utils.hpp"
#include <mtao/logging/logger.hpp>

namespace mandoline {






    Eigen::SparseMatrix<double> grid_boundary(int ni, int nj, bool dirichlet_boundary) {
        int system_size = ni*nj;
        int usize = (ni+1) * nj;
        int vsize = (nj+1) * ni;

        auto u_ind = [ni,nj](int i, int j) {
            return i + j*(ni + 1);
        };

        auto v_ind = [ni,nj](int i, int j) {
            return i + j*ni + (ni + 1)*nj;
        };

        std::vector< Eigen::Triplet<double> > triplets;

        Eigen::SparseMatrix<double> D(usize + vsize, system_size);
        int offset = usize;
        for (int j = 0; j < nj ; ++j) {
            for (int i = 0; i < ni ; ++i) {
                int index = i + ni*j;
                if(!(dirichlet_boundary && i == 0)) {
                    triplets.push_back(Eigen::Triplet<double>( u_ind(i,j), index, 1.0));
                }
                if(!(dirichlet_boundary && i == ni-1)) {
                    triplets.push_back(Eigen::Triplet<double>( u_ind(i+1,j), index, -1.0));
                }


                if(!(dirichlet_boundary && j == 0)) {
                    triplets.push_back(Eigen::Triplet<double>(  v_ind(i,j), index, 1.0));
                }
                if(!(dirichlet_boundary && j == ni-1)) {
                    triplets.push_back(Eigen::Triplet<double>(  v_ind(i,j+1), index, -1.0));
                }

            }
        }
        D.setFromTriplets(triplets.begin(), triplets.end());
        return D;
    }
    Eigen::VectorXd grid_edgeratio(int ni, int nj, const Eigen::VectorXd& phi) {
        int system_size = ni*nj;
        int usize = (ni+1) * nj;
        int vsize = (nj+1) * ni;

        auto p_ind = [ni,nj](int i, int j) {
            return i + j*(ni + 1);
        };
        auto u_ind = [ni,nj](int i, int j) {
            return i + j*(ni + 1);
        };

        auto v_ind = [ni,nj](int i, int j) {
            return i + j*ni + (ni + 1)*nj;
        };


        auto p = [&](int i, int j) -> double{
            return phi(p_ind(i,j));
        };


        Eigen::VectorXd M(usize + vsize);
        int offset = usize;
        for (int j = 0; j < nj ; ++j) {
            for (int i = 0; i < ni+1 ; ++i) {
                M(u_ind(i,j)) = levelset::fraction_inside(p(i,j),p(i,j+1));
            }
        }
        for (int j = 0; j < nj+1 ; ++j) {
            for (int i = 0; i < ni ; ++i) {
                M(v_ind(i,j)) = levelset::fraction_inside(p(i,j),p(i+1,j));
            }
        }
        return M;
    }
    void fix_axis_aligned_indefiniteness(Eigen::SparseMatrix<double>& M, Eigen::VectorXd& rhs) {
        M = M.pruned();
        //replace empty rows/cols to make (at least) positive semi-definite
        for (int row = 0; row < M.outerSize(); ++row) {

            int count = M.innerVector(row).nonZeros();
            if (count == 0) {
                M.insert(row, row) = 1;
                rhs[row] = 0;
            }
        }
        M.makeCompressed();
    }
    Eigen::VectorXd solve_SPSD_system(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b) {
        Eigen::VectorXd x;
        auto solve = [&](auto&& solver) -> bool {
            solver.compute(A);
            if (solver.info() != Eigen::Success) {
                mtao::logging::warn()  << "Eigen factorization failed.\n";
                return false;
            }
            x = solver.solve(b);
            if (solver.info() != Eigen::Success) {
                mtao::logging::warn() << "Eigen solve failed.\n";
                return false;
            }
            return true;
        };

        if(!solve(Eigen::ConjugateGradient< Eigen::SparseMatrix<double> >()) ) {
            mtao::logging::warn() << "CG failed, falling back to LDLT!";
            if(!solve(Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> >()) ) {
                mtao::logging::error() << "Solver failed!";
            }
        }
        return x;
    }
}
