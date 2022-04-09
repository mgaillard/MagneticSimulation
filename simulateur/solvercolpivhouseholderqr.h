#pragma once

#include <Eigen/Core>

#include "solverdense.h"

class SolverColPivHouseholderQR : public SolverDense
{
public:
    SolverColPivHouseholderQR();

    Eigen::MatrixXd computeSolution(const Eigen::VectorXd& columnVector) override;
};
