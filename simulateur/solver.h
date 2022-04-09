#pragma once

#include <Eigen/Core>

class Solver
{
public:
    Solver();

    virtual ~Solver() = default;

    void setSize(int rows, int cols);

    /**
     * \brief Generate the matrix A that contains the coefficients of the system of linear equations to solve
     * \param vecAlpha Vector alpha (upper diagonal)
     * \param vecBeta  Vector beta (lower diagonal)
     * \param vecGamma Vector gamma (diagonal)
     * \param vecDelta1 Vector delta1 (lower band)
     * \param vecDelta2 Vector delta2 (upper band)
     */
    virtual void generatePrincipalMatrix(const Eigen::VectorXd& vecAlpha,
                                         const Eigen::VectorXd& vecBeta,
                                         const Eigen::VectorXd& vecGamma,
                                         const Eigen::VectorXd& vecDelta1,
                                         const Eigen::VectorXd& vecDelta2) = 0;

    /**
     * \brief Compute the solution of the system of linear equations Ax=b
     * \param columnVector The column vector b of the system of linear equations
     * \return The solution x of the system of linear equations
     */
    virtual Eigen::MatrixXd computeSolution(const Eigen::VectorXd& columnVector) = 0;

protected:

    int m_rows;
    int m_cols;
};
