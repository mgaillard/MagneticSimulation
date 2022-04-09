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

    void draw(QPainter& painter, int offsetX, int offsetY) override;

private:
    double m_radius;
};
