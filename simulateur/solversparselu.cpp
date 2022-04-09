#include "solversparselu.h"

#include <spdlog/spdlog.h>

SolverSparseLU::SolverSparseLU() : SolverSparse()
{
}

Eigen::MatrixXd SolverSparseLU::computeSolution(const Eigen::VectorXd &columnVector)
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

    solver.analyzePattern(m_principalMatrix);
    solver.factorize(m_principalMatrix);

    Eigen::MatrixXd solution = solver.solve(columnVector);

    // Compute the relative error of the solution of Ax=b (norm() is L2 norm)
    const double relativeError = (m_principalMatrix * solution - columnVector).norm() / columnVector.norm();
    spdlog::debug("The relative error is: {}", relativeError);

    return solution;
}
