#ifndef SOLVEURCOLPIVHOUSEHOLDERQR_H
#define SOLVEURCOLPIVHOUSEHOLDERQR_H

#include <Eigen/Eigen>

#include "solveurdense.h"

class SolveurColPivHouseholderQR : public SolveurDense
{
public:
    SolveurColPivHouseholderQR();
    Eigen::MatrixXd calculerSolution(const Eigen::VectorXd &vecteurSecondMembre);
};

#endif // SOLVEURCOLPIVHOUSEHOLDERQR_H
