#include "solveursparselu.h"

SolveurSparseLU::SolveurSparseLU() : SolveurSparse()
{
}

Eigen::MatrixXd SolveurSparseLU::calculerSolution(const Eigen::VectorXd &vecteurSecondMembre)
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    solver.analyzePattern(matricePrincipale);
    solver.factorize(matricePrincipale);
    Eigen::MatrixXd vecteurSolution = solver.solve(vecteurSecondMembre);
    double relative_error = (matricePrincipale*vecteurSolution - vecteurSecondMembre).norm() / vecteurSecondMembre.norm(); // norm() is L2 norm
    std::cout << "The relative error is:\n" << relative_error << std::endl;
    return vecteurSolution;
}
