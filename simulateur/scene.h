#pragma once

#include <QRect>
#include <QList>

#include <Eigen/Core>

#include "shape.h"

class Scene
{
public:
    Scene();

    void setRectangleScene(const QRect &r);
    const QRect& getRectangleScene() const;

    int resolutionHeight() const;
    int resolutionWidth() const;
    int resolutionTotal() const;

    void setStep(double s);
    const double& getStep() const;
    double getSqStep() const;

    Eigen::Vector2d lowerBoundMeters() const;
    Eigen::Vector2d higherBoundMeters() const;
    bool isInside(const Eigen::Vector2d& p) const;

    /**
     * \brief Convert coordinates in meters to coordinates in the grid (pixels)
     * \param p Coordinates in meters
     * \return Coordinates in the grid
     */
    Eigen::Vector2d convertToGridCoords(const Eigen::Vector2d& p) const;

    /**
     * \brief Convert coordinates in the grid (pixels) to coordinates in meters
     * \param p Coordinates in the grid
     * \return Coordinates in meters
     */
    Eigen::Vector2d convertFromGridCoords(const Eigen::Vector2d& p) const;

    void setShapes(const ShapeList& shapes);
    const ShapeList& getShapes() const;

private:
    QRect m_rectangleScene;
    double m_step;
    ShapeList m_shapes;
};
