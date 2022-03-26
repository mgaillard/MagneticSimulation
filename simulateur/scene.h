#pragma once

#include <QRect>
#include <QList>

#include "forme.h"

class Scene
{
public:
    Scene();

    void setRectangleScene(const QRect &r);
    const QRect& getRectangleScene() const;
    int resolutionHeight() const;
    int resolutionWidth() const;
    int resolutionTotal() const;

    void setPas(const double& p);
    const double& getPas() const;
    double getSqPas() const;

    void setShapes(const FormeList& shapes);
    const FormeList& getShapes() const;

private:
    QRect m_rectangleScene;
    double m_pas;
    FormeList m_shapes;
};
