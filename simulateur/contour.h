#pragma once

#include <vector>

#include <Eigen/Core>

class Contour
{
public:
    Contour() = default;

    /**
     * \brief Add an iso-contour computed from a matrix
     * \param matrix The matrix from which the contour should be drawn
     * \param k The iso-value of the contour
     */
    void addContour(const Eigen::MatrixXd& matrix, double k);

    const std::vector<Eigen::Vector2d>& points() const;

    const std::vector<std::vector<int>>& curves() const;

private:
    void addLine(const Eigen::Vector2d& a, const Eigen::Vector2d& b);

    std::tuple<bool, int> addPoint(const Eigen::Vector2d& p);

    int findCurveWithPoint(int index) const;

    bool atBeginningOfCurve(int pointIndex, int curveIndex) const;

    bool atEndOfCurve(int pointIndex, int curveIndex) const;

    std::vector<Eigen::Vector2d> m_points;
    std::vector<std::vector<int>> m_curves;
};
