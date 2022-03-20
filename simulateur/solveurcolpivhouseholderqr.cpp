#include "solveurcolpivhouseholderqr.h"

SolveurColPivHouseholderQR::SolveurColPivHouseholderQR() : SolveurDense()
{
}

Eigen::MatrixXd SolveurColPivHouseholderQR::calculerSolution(const Eigen::VectorXd &vecteurSecondMembre)
{
    Eigen::MatrixXd vecteurSolution = matricePrincipale.colPivHouseholderQr().solve(vecteurSecondMembre);
    double relative_error = (matricePrincipale*vecteurSolution - vecteurSecondMembre).norm() / vecteurSecondMembre.norm(); // norm() is L2 norm
    std::cout << "The relative error is:\n" << relative_error << std::endl;
    return vecteurSolution;
}
