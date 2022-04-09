#pragma once

#include <QColor>

#include <Eigen/Core>

/**
 * \brief Strategy to rescale matrices when exporting to an image file.
 */
enum class ImageScalingStrategy
{
	/**
	 * \brief Rescale so that 0.0 corresponds to 0.5 intensity
	 *        Also, 0.0 or 1.0 correspond to the maximum absolute intensity
	 */
	ZeroCentered,

	/**
	 * \brief Rescale so that 0.0 corresponds to 0.0 intensity
	 */
	ZeroMinimum
};

/**
 * \brief Save a Eigen matrix in image format 
 * \param filename Path to the file
 * \param matrix The matrix to save
 * \param strategy Strategy used to scale the matrix to a color map
 * \return True if the image was successfully saved
 */
bool exportScalarMatrixImage(const QString& filename,
                             Eigen::MatrixXd& matrix,
                             const ImageScalingStrategy& strategy);
