#include "magneticsimulation.h"

#include <fstream>
#include <chrono>

#include <spdlog/spdlog.h>

#include "constants.h"
#include "export_image.h"
#include "export_vtk.h"

MagneticSimulation::MagneticSimulation(const QString& filename)
{
    // Load simulation definition file
    m_scene = InterfaceXml::loadFile(filename);
}

void MagneticSimulation::initialize()
{
    const auto width = m_scene.resolutionWidth();
    const auto height = m_scene.resolutionHeight();
    const auto total = m_scene.resolutionTotal();

    m_matPermeability = Eigen::MatrixXd::Constant(height, width, constants::mu_0);
    m_matCurrentDensity = Eigen::MatrixXd::Constant(height, width, 0.0);
    m_matDerZ = Eigen::MatrixXd::Constant(height, width, 0.0);
    m_matDerR = Eigen::MatrixXd::Constant(height, width, 0.0);
    m_vecSecondColumn = Eigen::MatrixXd::Constant(total, 1, 0.0);
    m_vecSolution = Eigen::MatrixXd::Constant(total, 1, 0.0);
    m_matBr = Eigen::MatrixXd::Constant(height, width, 0.0);
    m_matBz = Eigen::MatrixXd::Constant(height, width, 0.0);
    m_matB = Eigen::MatrixXd::Constant(height, width, 0.0);
    m_solver.setSize(height, width);
}

bool MagneticSimulation::checkValid() const
{
    // Check that there is at least one shape in the simulation
    if (m_scene.getShapes().empty())
    {
        return false;
    }

    return true;
}

void MagneticSimulation::runSimulation()
{
    initialize();
	fillMatPermeability();
    smoothMatPermeability();
    fillMatCurrentDensity();
    computeCurrentDerivatives();
    convertMatricesToVectors();
    computeSecondColumnVector();
    prepareSolver();
    computeSolution();
    convertSolutionToMatrix();
    computeMagneticField();
}

void MagneticSimulation::fillMatPermeability()
{
    const int rows = m_scene.resolutionHeight();
    const int cols = m_scene.resolutionWidth();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            QPoint point(j, i);

            ShapeList::iterator f;
            for (const auto& shape : m_scene.getShapes())
            {
                if (shape->isInside(point))
                {
                    // The shape contains the relative permeability
                    // We multiply it to the already existing permeability in the cell
                    m_matPermeability(i, j) *= shape->getPermeability();
                }
            }
        }
    }

    // We take the inverse of the permeability to get the eta matrix
    m_matPermeability = m_matPermeability.cwiseInverse();
}

void MagneticSimulation::smoothMatPermeability()
{
    // Fill the gaussian filter matrix
    constexpr int filterSize = 10;
    const double alpha = 0.1;
    Eigen::MatrixXd filter = Eigen::MatrixXd::Constant(filterSize, filterSize, 0.0);
    for (int i = 0; i < filterSize; i++)
    {
	    for (int j = 0; j < filterSize; j++)
	    {
		    filter(i, j) = std::exp(-alpha * (std::pow(i - filterSize / 2, 2)
		                                              + std::pow(j - filterSize / 2, 2)));
	    }
    }

    // Normalize the filter so that the sum of coefficients is 1.0
    filter /= filter.sum();

    const int rows = m_scene.resolutionHeight();
    const int cols = m_scene.resolutionWidth();

    Eigen::MatrixXd newMatPermeability = Eigen::MatrixXd::Constant(rows, cols, 0.0);

	constexpr int height = filterSize / 2;
    constexpr int width = filterSize / 2;

    // Compute a convolution with the filter
    for (int i = height; i < rows - height; i++)
    {
	    for (int j = width; j < cols - width; j++)
	    {
		    for (int c = 0; c < filterSize; c++)
		    {
			    for (int d = 0; d < filterSize; d++)
			    {
				    newMatPermeability(i, j) += m_matPermeability(i - height + c, j - width + d) * filter(c, d);
			    }
		    }
	    }
    }

    // Take care of boundaries
    for (int i = 0; i < height; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    newMatPermeability(i, j) = newMatPermeability(height, j);
	    }
    }
    for (int i = rows - height; i < rows; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    newMatPermeability(i, j) = newMatPermeability(rows - height - 1, j);
	    }
    }
    for (int i = 0; i < rows; i++)
    {
	    for (int j = 0; j < width; j++)
	    {
		    newMatPermeability(i, j) = newMatPermeability(i, width);
	    }
    }
    for (int i = 0; i < rows; i++)
    {
	    for (int j = cols - width; j < cols; j++)
	    {
		    newMatPermeability(i, j) = newMatPermeability(i, cols - width - 1);
	    }
    }

    // Replace the permeability matrix with the smoothed version
    m_matPermeability = newMatPermeability;
}

void MagneticSimulation::fillMatCurrentDensity()
{
    const int rows = m_scene.resolutionHeight();
    const int cols = m_scene.resolutionWidth();

	// This table contains the surface of all shapes in the scene
    std::vector<double> surfaces;
    for (const auto& shape : m_scene.getShapes())
    {
        // Surface in number of cells
        const auto surfaceInCells = shape->computeSurface(m_scene.getRectangleScene());
        // Surface in meters
        const auto surfaceInMeters = static_cast<double>(surfaceInCells) * m_scene.getSqStep();

        surfaces.push_back(surfaceInMeters);
    }

    // Fill the matrix of density of currents
    // For each cell of the simulation check if it is encompassed by one of the shape in the scene
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            const QPoint point(j, i);

            for (unsigned int k = 0; k < m_scene.getShapes().size(); k++)
            {
                if (m_scene.getShapes()[k]->isInside(point))
                {
                    // The current density in the cell is the current in ampere divided by the surface in meters
                    m_matCurrentDensity(i, j) += m_scene.getShapes()[k]->getCurrent() / surfaces[k];
                }
            }
        }
    }
}

void MagneticSimulation::computeCurrentDerivatives()
{
    const int rows = m_scene.resolutionHeight();
    const int cols = m_scene.resolutionWidth();
    const auto twoSteps = 2.0 * m_scene.getStep();

    for (int i = 1; i < rows - 1; i++)
    {
	    for (int j = 1; j < cols - 1; j++)
	    {
		    m_matDerZ(i, j) = (m_matPermeability(i + 1, j) - m_matPermeability(i - 1, j)) / twoSteps;
		    m_matDerR(i, j) = (m_matPermeability(i, j + 1) - m_matPermeability(i, j - 1)) / twoSteps;
	    }
    }
}


void MagneticSimulation::convertMatricesToVectors()
{
    m_matPermeability.transposeInPlace();
    m_matPermeability.resize(m_scene.resolutionTotal(), 1);
    m_matCurrentDensity.transposeInPlace();
    m_matCurrentDensity.resize(m_scene.resolutionTotal(), 1);
    m_matDerZ.transposeInPlace();
    m_matDerZ.resize(m_scene.resolutionTotal(), 1);
    m_matDerR.transposeInPlace();
    m_matDerR.resize(m_scene.resolutionTotal(), 1);
}

void MagneticSimulation::convertSolutionToMatrix()
{
    m_vecSolution.resize(m_scene.resolutionWidth(), m_scene.resolutionHeight());
    m_vecSolution.transposeInPlace();

    // Transpose back to matrices
    m_matPermeability.resize(m_scene.resolutionWidth(), m_scene.resolutionHeight());
    m_matPermeability.transposeInPlace();
    m_matCurrentDensity.resize(m_scene.resolutionWidth(), m_scene.resolutionHeight());
    m_matCurrentDensity.transposeInPlace();
}

void MagneticSimulation::prepareSolver()
{
    // Resolution constants
    const int rows = m_scene.resolutionHeight();
    const int cols = m_scene.resolutionWidth();
    const int nbCells = rows * cols;

    // Dimension constants: h: step size on axis r, k: step size on axis z
    const double h = m_scene.getStep();
    const double k = m_scene.getStep();
    const double hSq = m_scene.getSqStep();
    const double kSq = m_scene.getSqStep();

    // Initialize vectors needed for the solver
    Eigen::VectorXd vecAlpha = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecBeta = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecGamma = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecDelta1 = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecDelta2 = Eigen::VectorXd::Constant(nbCells, 0);

    for (int i = 0; i < rows; i++)
    {
	    for (int j = 1; j < cols - 1; j++)
	    {
		    const int n = cols * i + j;

		    vecAlpha(n + 1) = ((1.0 + 1.0 / (2 * j)) / hSq) + m_matDerR(n, 0) / (2.0 * h) / m_matPermeability(n, 0);
		    vecBeta(n - 1) = ((1.0 - 1.0 / (2 * j)) / hSq) - m_matDerR(n, 0) / (2.0 * h) / m_matPermeability(n, 0);
		    vecGamma(n) = (-2.0 / hSq - 2.0 / kSq - 1.0 / (j * j * hSq)) + m_matDerR(n, 0) / (j * h) / m_matPermeability(n, 0);

	    	if (j == 1)
		    {
			    vecAlpha(n) = 0.0;
			    vecBeta(n - 1) = 0.0;
			    vecGamma(n - 1) = 1.0;
		    }
	    }
    }

    for (int i = 1; i < rows; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    const int n = cols * i + j;

		    vecDelta2(n) = 1.0 / kSq + m_matDerZ(n - cols, 0) / (2 * k) / m_matPermeability(n - cols, 0);

	    	if (j == 0)
		    {
			    vecDelta2(n) = 0.0;
		    }
	    }
    }

    for (int i = 0; i < rows - 1; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    const int n = cols * i + j;

		    vecDelta1(n) = 1.0 / kSq - m_matDerZ(n + cols, 0) / (2 * k) / m_matPermeability(n + cols, 0);

	    	if (j == 0)
		    {
			    vecDelta1(n) = 0.0;
		    }
	    }
    }

    m_solver.generatePrincipalMatrix(vecAlpha, vecBeta, vecGamma, vecDelta1, vecDelta2);
}


void MagneticSimulation::computeSecondColumnVector()
{
	for (int i = 0; i < m_scene.resolutionTotal(); i++)
	{
		m_vecSecondColumn(i) = -m_matCurrentDensity(i, 0) / m_matPermeability(i, 0);
	}
}

void MagneticSimulation::computeSolution()
{
    const auto startTime = std::chrono::steady_clock::now();
    m_vecSolution = m_solver.computeSolution(m_vecSecondColumn);
    const auto endTime = std::chrono::steady_clock::now();

    spdlog::info("solving time: {} ms", std::chrono::duration<double, std::milli>(endTime - startTime).count());
}

void MagneticSimulation::computeMagneticField()
{
    // Resolution constants
    const int rows = m_scene.resolutionHeight();
    const int cols = m_scene.resolutionWidth();

    // Dimension constants: h: step size on axis r, k: step size on axis z
    const double h = m_scene.getStep();
    const double k = m_scene.getStep();

    // The center of the matrices Br and Bz are equal to B = rot(A).
    for (int i = 1; i < rows - 1; i++)
    {
	    for (int j = 1; j < cols - 1; j++)
	    {
		    m_matBz(i, j) = m_vecSolution(i, j) / (h * j) + (m_vecSolution(i, j + 1) - m_vecSolution(i, j - 1)) / (2.0 * h);
		    m_matBr(i, j) = (-m_vecSolution(i + 1, j) + m_vecSolution(i - 1, j)) / (2.0 * k);
	    }
    }

    // Take care of boundaries
    for (int i = 0; i < rows; i++)
    {
	    m_matBz(i, 0) = m_matBz(i, 1);
	    m_matBr(i, 0) = m_matBr(i, 1);
	    m_matBz(i, cols - 1) = m_matBz(i, cols - 2);
	    m_matBr(i, cols - 1) = m_matBr(i, cols - 2);
    }

    for (int j = 0; j < cols; j++)
    {
	    m_matBz(0, j) = m_matBz(1, j);
	    m_matBr(0, j) = m_matBr(1, j);
	    m_matBz(rows - 1, j) = m_matBz(rows - 2, j);
	    m_matBr(rows - 1, j) = m_matBr(rows - 2, j);
    }

    // Compute the magnitude of the magnetic field
    for (int i = 0; i < rows; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    m_matB(i, j) = std::sqrt(std::pow(m_matBr(i, j), 2) + std::pow(m_matBz(i, j), 2));
	    }
    }
}


Eigen::MatrixXd MagneticSimulation::matrixSymmetric(const Eigen::MatrixXd& matrix, double coefficient) const
{
	Eigen::MatrixXd matrixSym = Eigen::MatrixXd::Constant(m_scene.resolutionHeight(),
	                                                      2 * m_scene.resolutionWidth() - 1,
	                                                      0);

	// On the right part, we simply copy the matrix
	for (int i = 0; i < m_scene.resolutionHeight(); i++)
	{
		for (int j = m_scene.resolutionWidth() - 1; j < 2 * m_scene.resolutionWidth() - 1; j++)
		{
			matrixSym(i, j) = matrix(i, j - m_scene.resolutionWidth() + 1);
		}
	}

	// On the left part, we copy the opposite of the matrix
	for (int i = 0; i < m_scene.resolutionHeight(); i++)
	{
		for (int j = 0; j < m_scene.resolutionWidth() - 1; j++)
		{
			matrixSym(i, j) = coefficient * matrix(i, m_scene.resolutionWidth() - 1 - j);
		}
	}

	return matrixSym;
}

void MagneticSimulation::saveResults(
    const QString& fichierMatriceA,
    const QString& fichierMatriceBr,
    const QString& fichierMatriceBz,
    const QString& fichierMatriceB,
    const QString& fileMatrixAWithContour,
    const QString& fileMatrixBWithStreamlines)
{
    const auto matrixAFull = matrixSymmetric(m_vecSolution, -1.0);
    const auto matrixBrFull = matrixSymmetric(m_matBr, -1.0);
    const auto matrixBzFull = matrixSymmetric(m_matBz, 1.0);
    const auto matrixBFull = matrixSymmetric(m_matB, 1.0);

    exportScalarMatrixImage(fichierMatriceA, matrixAFull, ImageScalingStrategy::ZeroCentered);
    exportScalarMatrixImage(fichierMatriceBr, matrixBrFull, ImageScalingStrategy::ZeroCentered);
    exportScalarMatrixImage(fichierMatriceBz, matrixBzFull, ImageScalingStrategy::ZeroCentered);
    exportScalarMatrixImage(fichierMatriceB, matrixBFull, ImageScalingStrategy::ZeroMinimum);

    exportScalarMatrixWithSceneImage(fileMatrixAWithContour,
                                     m_vecSolution,
                                     ImageScalingStrategy::ZeroCentered,
                                     m_scene,
                                     0);

    exportScalarMatrixWithStreamlines(fileMatrixBWithStreamlines,
                                      m_matB,
                                      m_matBr,
                                      m_matBz,
                                      ImageScalingStrategy::ZeroMinimum,
                                      m_scene);

    // Save the values of the magnetic field on the line r=0.0
    std::ofstream file("matrixbcenter.txt", std::ios::out | std::ios::trunc);

    if(file)
    {
        for (int i = 0;i < m_scene.resolutionHeight();i++)
        {
            file << m_matB(i, 0) << std::endl;
        }
        file.close();
    }

    // Output matrices in VTK format
    exportScalarMatrixVtk("mu.vtk", m_matPermeability, m_scene.getStep(), m_scene.getStep());
    exportScalarMatrixVtk("i.vtk", m_matCurrentDensity, m_scene.getStep(), m_scene.getStep());
	exportScalarMatrixVtk("a.vtk", m_vecSolution, m_scene.getStep(), m_scene.getStep());
    exportVectorMatrixVtk("b.vtk", m_matBr, m_matBz, m_scene.getStep(), m_scene.getStep());
}
