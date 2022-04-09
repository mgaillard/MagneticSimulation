#include "export_image.h"

#include <QImage>
#include <QPainter>

#include "shapecircle.h"
#include "shapepolygon.h"

// ---------------------------------------- Private functions ----------------------------------------

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

QImage generateScalarMatrixImage(const Eigen::MatrixXd& matrix, const ImageScalingStrategy& strategy)
{
	const auto cols = static_cast<int>(matrix.cols());
	const auto rows = static_cast<int>(matrix.rows());

	const auto minimum = matrix.minCoeff();
	const auto maximum = matrix.maxCoeff();

	double rescaleMinimum = 0.0;
	double rescaleMaximum = 1.0;
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

	return image;
}

// ---------------------------------------- Public functions ----------------------------------------

bool exportScalarMatrixImage(
	const QString& filename,
	const Eigen::MatrixXd& matrix,
	const ImageScalingStrategy& strategy)
{
	const auto image = generateScalarMatrixImage(matrix, strategy);

	return image.save(filename);
}

bool exportScalarMatrixWithSceneImage(
	const QString& filename,
	const Eigen::MatrixXd& matrix,
	const ImageScalingStrategy& strategy,
	const Scene& scene,
	int offsetX)
{
	// Pen width
	const int width = 3;

	// The background image is the export of the matrix
	auto image = generateScalarMatrixImage(matrix, strategy);

	QPainter painter;
	painter.begin(&image);
	
	for (const auto& shape : scene.getShapes())
	{
		if (shape->getCurrent() != 0.0)
		{
			// If the shape has current, display in yellow
			painter.setPen(QPen(QColor(255, 255, 0), width, Qt::SolidLine));
		}
		else
		{
			// If the shape has no current, display in red
			painter.setPen(QPen(QColor(255, 0, 0), width, Qt::SolidLine));
		}

		shape->draw(painter, offsetX, 0);
	}

	painter.end();

	return image.save(filename);
}
