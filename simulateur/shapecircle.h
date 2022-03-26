#pragma once

#include <QPoint>
#include <QRect>

#include "shape.h"

class ShapeCircle final : public Shape
{
public:
    ShapeCircle(const QPoint &c, double r);
    ~ShapeCircle() override = default;

    bool isInside(const QPoint &p) override;
    QRect getBounds() override;

	void setRadius(double r);
    const double& getRadius() const;

private:
    double m_radius;
};
