#include "export_vtk.h"

#include <fstream>

#include <spdlog/spdlog.h>

bool exportScalarMatrixVtk(const std::string& filename, const Eigen::MatrixXd& matrix, double h, double k)
{
	std::ofstream out(filename);

	if (!out.is_open())
	{
		spdlog::error("Could not open file: {}", filename);
		return false;
	}

	out << "# vtk DataFile Version 3.0\n\n";

	out << "ASCII\n";
	out << "DATASET STRUCTURED_GRID\n";
	out << "DIMENSIONS " << matrix.cols() << " " << matrix.rows() << " 1\n";
	out << "POINTS " << matrix.cols() * matrix.rows() << " float\n";

	for (unsigned int i = 0; i < matrix.rows(); i++)
	{
		for (unsigned int j = 0; j < matrix.cols(); j++)
		{
			out << static_cast<double>(j) * k << " " << static_cast<double>(i) * h << " 0.0\n";
		}
	}

	out << "POINT_DATA " << matrix.cols() * matrix.rows() << "\n";
	out << "SCALARS value float\n";
	out << "LOOKUP_TABLE default\n";

	for (unsigned int i = 0; i < matrix.rows(); i++)
	{
		for (unsigned int j = 0; j < matrix.cols(); j++)
		{
			out << matrix(i, j) << "\n";
		}
	}

	out.close();

	return true;
}

bool exportVectorMatrixVtk(
	const std::string& filename,
	const Eigen::MatrixXd& matrixR,
	const Eigen::MatrixXd& matrixZ,
	double h,
	double k)
{
	if (matrixR.rows() != matrixZ.rows())
	{
		spdlog::error("Matrices do not have the same number of rows");
		return false;
	}

	if (matrixR.cols() != matrixZ.cols())
	{
		spdlog::error("Matrices do not have the same number of cols");
		return false;
	}

	std::ofstream out(filename);

	if (!out.is_open())
	{
		spdlog::error("Could not open file: {}", filename);
		return false;
	}

	out << "# vtk DataFile Version 3.0\n\n";

	out << "ASCII\n";
	out << "DATASET STRUCTURED_GRID\n";
	out << "DIMENSIONS " << matrixR.cols() << " " << matrixR.rows() << " 1\n";
	out << "POINTS " << matrixR.cols() * matrixR.rows() << " float\n";

	for (unsigned int i = 0; i < matrixR.rows(); i++)
	{
		for (unsigned int j = 0; j < matrixR.cols(); j++)
		{
			out << static_cast<double>(j) * k << " " << static_cast<double>(i) * h << " 0.0\n";
		}
	}

	out << "POINT_DATA " << matrixR.cols() * matrixR.rows() << "\n";
	out << "VECTORS vector float\n";

	for (unsigned int i = 0; i < matrixR.rows(); i++)
	{
		for (unsigned int j = 0; j < matrixR.cols(); j++)
		{
			out << matrixR(i, j) << " " << matrixZ(i, j) << " 0.0\n";
		}
	}

	out.close();

	return true;
}
