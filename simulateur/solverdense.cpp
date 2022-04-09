#include "solverdense.h"

SolverDense::SolverDense() : Solver()
{
}

void SolverDense::generatePrincipalMatrix(
	const Eigen::VectorXd& vecAlpha,
	const Eigen::VectorXd& vecBeta,
	const Eigen::VectorXd& vecGamma,
	const Eigen::VectorXd& vecDelta1,
	const Eigen::VectorXd& vecDelta2)
{
	const auto totalSize = m_rows * m_cols;
	m_principalMatrix = Eigen::MatrixXd::Constant(totalSize, totalSize, 0.0);

	for (int i = 0; i < totalSize; i++)
	{
		m_principalMatrix(i, i) = vecGamma(i);

		if (i >= 1)
		{
			m_principalMatrix(i, i - 1) = vecBeta(i - 1);
		}
		if (i >= m_cols)
		{
			m_principalMatrix(i, i - m_cols) = vecDelta1(i - m_cols);
		}
		if (i <= totalSize - 2)
		{
			m_principalMatrix(i, i + 1) = vecAlpha(i + 1);
		}
		if (i <= totalSize - 1 - m_cols)
		{
			m_principalMatrix(i, i + m_cols) = vecDelta2(i + m_cols);
		}
	}
}
