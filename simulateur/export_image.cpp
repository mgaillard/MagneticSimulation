#include "export_image.h"

#include <QImage>

/**
 * \brief Compute the color corresponding to an intensity (Matlab JET color map)
 * \param intensity The intensity in the range [0.0, 1.0]
 * \return A color as a QColor object
 */
QColor getColorByIntensity(const double intensity)
{
	double red = 0.0;
	double green = 0.0;
	double blue = 0.0;

	if (intensity >= 0.0 && intensity < 0.125)
	{
		red = 0.0;
		green = 0.0;
		blue = 4 * intensity + 0.5;
	}
	else if (intensity >= 0.125 && intensity < 0.375)
	{
		red = 0;
		green = 4 * intensity - 0.5;
		blue = 1;
	}
	else if (intensity >= 0.375 && intensity < 0.625)
	{
		red = 4 * intensity - 1.5;
		green = 1;
		blue = -4 * intensity + 2.5;
	}
	else if (intensity >= 0.625 && intensity < 0.875)
	{
		red = 1;
		green = -4 * intensity + 3.5;
		blue = 0;
	}
	else if (intensity >= 0.875 && intensity <= 1.0)
	{
		red = -4 * intensity + 4.5;
		green = 0;
		blue = 0;
	}

	return {
		static_cast<int>(red * 255.0),
		static_cast<int>(green * 255.0),
		static_cast<int>(blue * 255.0)
	};
}

bool exportScalarMatrixImage(
	const QString &filename,
	Eigen::MatrixXd &matrix,
	const ImageScalingStrategy& strategy)
{
	const auto cols = static_cast<int>(matrix.cols());
	const auto rows = static_cast<int>(matrix.rows());
	
	const auto minimum = matrix.minCoeff();
	const auto maximum = matrix.maxCoeff();

	double rescaleMinimum, rescaleMaximum;
	switch (strategy)
	{
	case ImageScalingStrategy::ZeroCentered:
		rescaleMaximum = std::max(std::abs(minimum), std::abs(maximum));
		rescaleMinimum = -rescaleMaximum;
		break;
	case ImageScalingStrategy::ZeroMinimum:
		rescaleMinimum = 0.0;
		rescaleMaximum = maximum;
		break;
	}

	QImage image(cols, rows, QImage::Format_RGB32);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			// Remap the value in the matrix between 0.0 and 1.0
			auto value = (matrix(i, j) - rescaleMinimum) / (rescaleMaximum - rescaleMinimum);
			value = std::clamp(value, 0.0, 1.0);

			image.setPixel(QPoint(j, i), getColorByIntensity(value).rgb());
		}
	}

	return image.save(filename);
}
