#pragma once

#include <QString>

#include <Eigen/Eigen>

#include "scene.h"
#include "interfacexml.h"
#include "solversparselu.h"
#include "solvercolpivhouseholderqr.h"

class MagneticSimulation
{
public:
    MagneticSimulation(const QString& filename);

	bool checkValid() const;

    void runSimulation();

    void saveResults(const QString& fichierMatriceA,
                     const QString& fichierMatriceBr,
                     const QString& fichierMatriceBz,
                     const QString& fichierMatriceB,
                     const QString& fileMatrixAWithContour);
private:

    void initialize();

    void fillMatPermeability();

    void smoothMatPermeability();

    void fillMatCurrentDensity();

    void computeCurrentDerivatives();

    void convertMatricesToVectors();

    void convertSolutionToMatrix();

    void prepareSolver();

    void computeSecondColumnVector();

    void computeSolution();

    void computeMagneticField();

    Eigen::MatrixXd matrixSymmetric(const Eigen::MatrixXd& matrix, double coefficient) const;

    Scene m_scene;
    Eigen::MatrixXd m_matPermeability;
    Eigen::MatrixXd m_matCurrentDensity;
    Eigen::MatrixXd m_matDerZ;
    Eigen::MatrixXd m_matDerR;
    Eigen::VectorXd m_vecSecondColumn;
    Eigen::MatrixXd m_vecSolution;
    Eigen::MatrixXd m_matBr;
    Eigen::MatrixXd m_matBz;
    Eigen::MatrixXd m_matB;
    SolverSparseLU m_solver;
};
