#pragma once

#include <QColor>

#include <Eigen/Core>

bool exportScalarMatrixImage(const QString& filename, Eigen::MatrixXd& matrix);
