#ifndef SOLVEURSPARSELU_H
#define SOLVEURSPARSELU_H

#include <Eigen/Sparse>
#include <iostream>
#include "solveursparse.h"

class SolveurSparseLU : public SolveurSparse
{
public:
    SolveurSparseLU();
    Eigen::MatrixXd calculerSolution(const Eigen::VectorXd &vecteurSecondMembre);
};

#endif // SOLVEURSPARSELU_H
