#include "solversparse.h"


SolverSparse::SolverSparse() : Solver()
{
}

void SolverSparse::generatePrincipalMatrix(
	const Eigen::VectorXd& vecAlpha,
	const Eigen::VectorXd& vecBeta,
	const Eigen::VectorXd& vecGamma,
	const Eigen::VectorXd& vecDelta1,
	const Eigen::VectorXd& vecDelta2)
{
    const auto totalSize = m_rows * m_cols;
    m_principalMatrix.resize(totalSize, totalSize);

    // The sparse matrix is composed of triplets of coefficient: coordinates (i, j) and value
    std::vector<Eigen::Triplet<double>> tripletList;

    for (int i = 0; i < totalSize; i++)
    {
	    tripletList.emplace_back(i, i, vecGamma(i));
	    if (i >= 1)
	    {
		    tripletList.emplace_back(i, i - 1, vecBeta(i - 1));
	    }
	    if (i >= m_cols)
	    {
		    tripletList.emplace_back(i, i - m_cols, vecDelta1(i - m_cols));
	    }
	    if (i <= totalSize - 2)
	    {
		    tripletList.emplace_back(i, i + 1, vecAlpha(i + 1));
	    }
	    if (i <= totalSize - 1 - m_cols)
	    {
		    tripletList.emplace_back(i, i + m_cols, vecDelta2(i + m_cols));
	    }
    }

    m_principalMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    m_principalMatrix.makeCompressed();
}
