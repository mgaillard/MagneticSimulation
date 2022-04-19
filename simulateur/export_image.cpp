#include "export_image.h"

#include <QImage>
#include <QPainter>

#include <spdlog/spdlog.h>

#include "contour.h"
#include "streamlines.h"

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
	painter.setRenderHint(QPainter::Antialiasing, true);
	
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

	// TODO: add contour definition in XML
	// Generate 10 iso-values for contours based on minimum/maximum
	Contour contours;
	const auto minimum = matrix.minCoeff();
	const auto maximum = matrix.maxCoeff();
	constexpr int nbContours = 10;
	// Exclude the minimum contour and the maximum contour, because it does not make sense to show them
	for (int i = 1; i < nbContours - 1; i++)
	{
		const auto t = static_cast<double>(i) / (nbContours - 1);
		const double k = minimum + t * (maximum - minimum);
		contours.addContour(matrix, k);
	}

	// For contours we use white pen with width 1
	painter.setPen(QPen(QColor(255, 255, 255), 1, Qt::SolidLine));

	const auto contourPoints = contours.points();
	for (const auto& curve : contours.curves())
	{
		for (unsigned int i = 1; i < curve.size(); i++)
		{
			const QPointF start(
				contourPoints[curve[i - 1]].y(),
				contourPoints[curve[i - 1]].x()
			);

			const QPointF end(
				contourPoints[curve[i]].y(),
				contourPoints[curve[i]].x()
			);

			painter.drawLine(start, end);
		}
	}

	painter.end();

	return image.save(filename);
}

bool exportScalarMatrixWithStreamlines(
	const QString& filename,
	const Eigen::MatrixXd& matB,
	const Eigen::MatrixXd& matBr,
	const Eigen::MatrixXd& matBz,
	const ImageScalingStrategy& strategy,
	const Scene& scene)
{
	// The background image is the export of the matrix
	auto image = generateScalarMatrixImage(matB, strategy);

	QPainter painter;
	painter.begin(&image);
	painter.setRenderHint(QPainter::Antialiasing, true);


	StreamlinesPlacement streamlinesPlacement(scene, matBr, matBz);
	// Middle of the scene
	const auto minimum = scene.convertFromGridCoords({
		0,
		scene.resolutionHeight() / 2
	});
	// Right of the scene
	const auto maximum = scene.convertFromGridCoords({
		scene.resolutionWidth(),
		scene.resolutionHeight() / 2
	});
	constexpr int nbStreamlines = 10;
	// Exclude the minimum streamline and the maximum streamline, because it does not make sense to show them
	for (int i = 1; i < nbStreamlines - 1; i++)
	{
		const auto t = static_cast<double>(i) / (nbStreamlines - 1);
		const auto seed = minimum + t * (maximum - minimum);
		streamlinesPlacement.placeStreamline(seed);
	}

	// For contours we use white pen with width 1
	painter.setPen(QPen(QColor(255, 255, 255), 1, Qt::SolidLine));
	
	for (int i = 0; i < streamlinesPlacement.nbStreamlines(); i++)
	{
		const auto& streamline = streamlinesPlacement.streamline(i);

		for (unsigned int j = 1; j < streamline.size(); j++)
		{
			const auto startGrid = scene.convertToGridCoords(streamline[j - 1]);
			const auto endGrid = scene.convertToGridCoords(streamline[j]);

			const QPointF start(startGrid.x(), startGrid.y());
			const QPointF end(endGrid.x(), endGrid.y());

			painter.drawLine(start, end);
		}
	}

	painter.end();

	return image.save(filename);
}
