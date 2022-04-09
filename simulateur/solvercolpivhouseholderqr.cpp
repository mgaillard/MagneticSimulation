#include "solvercolpivhouseholderqr.h"

#include <spdlog/spdlog.h>

#include <Eigen/QR>

SolverColPivHouseholderQR::SolverColPivHouseholderQR() : SolverDense()
{
}

Eigen::MatrixXd SolverColPivHouseholderQR::computeSolution(const Eigen::VectorXd& columnVector)
{
	Eigen::MatrixXd solution = m_principalMatrix.colPivHouseholderQr().solve(columnVector);

	// Compute the relative error of the solution of Ax=b (norm() is L2 norm)
	const double relativeError = (m_principalMatrix * solution - columnVector).norm() / columnVector.norm();
	spdlog::debug("The relative error is: {}", relativeError);

	return solution;
}
