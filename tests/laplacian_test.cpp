#include <mandoline/mesh3.hpp>
#include <mtao/cmdline_parser.hpp>
#include <mtao/logging/logger.hpp>

#include <Eigen/Dense>
using namespace mandoline;







template <typename Matrix, typename Vector>
struct SparseLDLT
{
    typedef typename Matrix::Scalar Scalar;
    SparseLDLT() {}
    SparseLDLT(const Matrix & A)
    {
        // L=tril(A);
        L=A.template triangularView<Eigen::StrictlyLower>();//Don't copy the diagonal
        for(int i=0; i<L.rows(); ++i)
        {
            if(L.coeff(i,i)!=0)
                L.coeffRef(i,i)=0;
        }
        Dinv=D=A.diagonal();



        // for k=1:size(L,2)
        for(int k=0; k<A.rows(); ++k)//k is the column that we're infecting the remaining columns with
        {//L(:,k)



            //Solidify the current column values
            //==================================
            if(D(k)==0) continue;
            if(Dinv(k)<0.25*D(k))//If D has shrunk too much since it started
                Dinv(k)=1/D(k);
            else
                Dinv(k)=1/Dinv(k);
            L.innerVector(k) *= Dinv(k);

            //Add k terms to all of the following columns
            //===========================================
            for(typename Matrix::InnerIterator it(L,k); it; ++it)// -L(i,k)*D(k)*L(j,k)
            {
                int j = it.row();//j>k
                if(j<=k) continue;
                Scalar missing=0;
                Scalar multiplier=it.value();//L(j,k)*D(k)

                typename Matrix::InnerIterator k_it(L,k);
                typename Matrix::InnerIterator j_it(L,j);
                //move down teh column of L(:,k) to collect missing elements in the match with A(:,j)
                //i=k_it.row()

                while (k_it && k_it.row()<j){//L(i,k)
                    while(j_it)//L(i,j) occasionally
                    {
                        if(j_it.row() < k_it.row())
                            ++j_it;
                        else if(j_it.row() == k_it.row())//L(i,k) are L(i,j) are nonzero
                            break;
                        else
                        {
                            missing += k_it.value();//L(i,k) will fill something not in L(i,j)
                            break;
                        }
                    }
                    ++k_it;
                }


                if(k_it && j_it.row() == j)
                {
                    Dinv(j) -= it.value() * multiplier;
                }


                typename Matrix::InnerIterator j_it2(L,j);
                while(k_it && j_it2)
                {
                    if(j_it2.row() < k_it.row())
                        ++j_it2;
                    else if(j_it2.row() == k_it.row())//L(i,k) and L(i,j) are both nonzero, -=L(i,k)*L(j,k)*D(k)
                    {
                        j_it2.valueRef() -= multiplier * k_it.value() ;//k_it.value()=L(i,k)
                        ++j_it2;
                        ++k_it;
                    }
                    else
                    {
                        missing+=k_it.value();
                        ++k_it;
                    }
                }

                while(k_it)
                {
                    missing+=k_it.value();
                    ++k_it;
                }
                Dinv(j)-=0.97*missing*multiplier;
            }
        }

        /*
           std::cout << L << std::endl;
           */

    }
    void solve(const Vector & b, Vector & x)
    {
        x = L.template triangularView<Eigen::UnitLower>().solve(b);
        x.noalias() = x.cwiseProduct(Dinv);//safe beacuse it's a dot
        L.transpose().template triangularView<Eigen::UnitUpper>().solveInPlace(x);
    }
    Matrix getA()
    {
        Matrix
                A = L.template triangularView<Eigen::UnitLower>();
        A = A * D.asDiagonal();
        A = A * L.template triangularView<Eigen::UnitLower>().transpose();

        return A;
    }
private:
    Matrix L;
    Vector D,Dinv;
};

template <typename MatrixType, typename VectorType, typename Preconditioner>
struct PreconditionedConjugateGradient
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    PreconditionedConjugateGradient(const Matrix & A): A(A)
    {
        precond = Preconditioner(A);
    }
    auto solve(const Vector& b) {
        auto x = b.eval();
        x.setZero();
        Vector r = b-A*x;
        Vector z;
        precond.solve(r,z);
        Vector p = z;
        Vector Ap = A*p;
        Scalar rdz = r.dot(z);
        Scalar alpha, beta;
        auto error = [&]() { return r.template lpNorm<Eigen::Infinity>(); };

        uint iterations = 0;
        while(++iterations < 10 &&
                error() > epsilon)
        {
            alpha = (rdz)/(p.dot(Ap));
            x+=alpha * p;
            r-=alpha * Ap;
            precond.solve(r,z);
            beta=1/rdz;
            rdz = r.dot(z);
            beta*=rdz;
            p=z+beta*p;
            Ap=A*p;
        }
        return x;
    }
private:
    const Matrix & A;
    Preconditioner precond;

    Scalar epsilon = 1e-5;

};

template <typename Matrix, typename Vector>
auto ldlt_pcg_solver(const Matrix & A, const Vector& b)
{
    return PreconditionedConjugateGradient<Matrix,Vector, SparseLDLT<Matrix, mtao::Vector<typename Vector::Scalar, Vector::RowsAtCompileTime>>>(A);
    //auto solver = IterativeLinearSolver<PreconditionedConjugateGradientCapsule<Matrix,Vector, Preconditioner> >(A.rows(), 1e-5);
}
template <typename Matrix, typename Vector>
auto ldlt_pcg_solve(const Matrix & A, const Vector & b)
{
    auto solver = ldlt_pcg_solver(A,b);
    //auto solver = IterativeLinearSolver<PreconditionedConjugateGradientCapsule<Matrix,Vector, Preconditioner> >(A.rows(), 1e-5);
    return solver.solve(b);
}






using namespace mtao::logging;

double max_eigenvalue(const Eigen::SparseMatrix<double>& A) {
    mtao::VecXd xold = mtao::VecXd::Random(A.rows()).normalized();
    mtao::VecXd x = (A * xold).normalized();
    while((x - xold).norm() > 1e-5) {
        xold = x;
        x = (A * xold).normalized();
        std::cout << "Max:" << (x-xold).norm() << std::endl;
    }
    return (A*x).norm();
}

double min_eigenvalue(const Eigen::SparseMatrix<double>& A) {
    mtao::VecXd xold = mtao::VecXd::Random(A.rows()).normalized();
    std::cout << "Preparing solver" << std::endl;
    auto solver = ldlt_pcg_solver(A,xold);
    std::cout << "Done" << std::endl;

    mtao::VecXd x = (solver.solve(xold)).normalized();
    while((x - xold).norm() > 1e-5) {
        xold = x;
        x = (solver.solve(xold)).normalized();
        std::cout << "Min:" << (x-xold).norm() << std::endl;
    }
    return (A*x).norm();
}

int main(int argc, char * argv[]) {

    auto&& log = make_logger("profiler",mtao::logging::Level::All);
    mtao::CommandLineParser clp;
    clp.parse(argc, argv);

    if(clp.args().size() < 1) {
        fatal() << "No input mesh filename!";
        return {};
    }

    std::string input_cutmesh = clp.arg(0);

    CutCellMesh<3> ccm = CutCellMesh<3>::from_proto(input_cutmesh);

    auto B = ccm.boundary();
    mtao::VecXd DM = ccm.mesh_face_mask();
    B = DM.asDiagonal() * B;
    DM = mtao::VecXd::Ones(ccm.cell_size());
    for(auto&& [i,c]: mtao::iterator::enumerate(ccm.cells())) {
    }
    B = B * DM.asDiagonal();

    mtao::VecXd DH2 = ccm.dual_hodge2();
    Eigen::SparseMatrix<double> L = B.transpose() * DH2.asDiagonal() * B;


    double M = max_eigenvalue(L);
    double m = min_eigenvalue(L);
    
    std::cout << m << " < " << M  << std::endl;
    std::cout << "condition number: " << (M / m) << std::endl;
}
