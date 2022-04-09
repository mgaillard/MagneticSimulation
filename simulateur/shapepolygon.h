#pragma once

#include <QPoint>
#include <QRect>
#include <QPolygon>

#include "shape.h"

class ShapePolygon : public Shape
{
public:
    ShapePolygon(const QPoint &c, QPolygon p);
    ~ShapePolygon() override = default;

	bool isInside(const QPoint &p) override;
    QRect getBounds() override;

	void updateCenter();

	void setCenter(const QPoint &c) override;

    void setPolygon(const QPolygon &p);
    const QPolygon& getPolygon() const;

    void draw(QPainter& painter, int offsetX, int offsetY) override;

private:
    QPolygon m_polygon;
};
