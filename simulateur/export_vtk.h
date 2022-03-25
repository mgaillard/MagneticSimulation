#pragma once

#include <string>

#include <Eigen/Dense>

/**
 * \brief Save a Eigen matrix in VTK format 
 * \param filename Path to the file
 * \param matrix The matrix to save
 * \param h Step on the row axis in meters
 * \param k Step on the column axis in meters
 * \return True if successfully saved
 */
bool exportScalarMatrixVtk(const std::string& filename, const Eigen::MatrixXd& matrix, double h, double k);

/**
 * \brief Save a Eigen matrix in VTK format
 * \param filename Path to the file
 * \param matrixR The matrix containing the r coordinate of the vectors
 * \param matrixZ The matrix containing the z coordinate of the vectors
 * \param h Step on the row axis in meters
 * \param k Step on the column axis in meters
 * \return True if successfully saved
 */
bool exportVectorMatrixVtk(const std::string& filename,
                           const Eigen::MatrixXd& matrixR,
                           const Eigen::MatrixXd& matrixZ,
                           double h,
                           double k);
