#pragma once

#include <Eigen/SparseCore>

#include "solver.h"

class SolverSparse : public Solver
{
public:
    SolverSparse();

    void generatePrincipalMatrix(const Eigen::VectorXd& vecAlpha,
                                 const Eigen::VectorXd& vecBeta,
                                 const Eigen::VectorXd& vecGamma,
                                 const Eigen::VectorXd& vecDelta1,
                                 const Eigen::VectorXd& vecDelta2) override;

protected:

    Eigen::SparseMatrix<double> m_principalMatrix;
};
