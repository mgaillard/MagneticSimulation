#pragma once

#include <QRect>
#include <QList>

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

    void setShapes(const ShapeList& shapes);
    const ShapeList& getShapes() const;

private:
    QRect m_rectangleScene;
    double m_step;
    ShapeList m_shapes;
};
