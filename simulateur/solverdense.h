#pragma once

#include <Eigen/Core>

#include "solver.h"

class SolverDense : public Solver
{
public:
    SolverDense();

    void generatePrincipalMatrix(const Eigen::VectorXd& vecAlpha,
                                 const Eigen::VectorXd& vecBeta,
                                 const Eigen::VectorXd& vecGamma,
                                 const Eigen::VectorXd& vecDelta1,
                                 const Eigen::VectorXd& vecDelta2) override;

protected:

    Eigen::MatrixXd m_principalMatrix;
};
