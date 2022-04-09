#pragma once

#include <Eigen/Sparse>

#include "solversparse.h"

class SolverSparseLU : public SolverSparse
{
public:
    SolverSparseLU();

    Eigen::MatrixXd computeSolution(const Eigen::VectorXd & columnVector) override;
};
