#include "shapepolygon.h"

ShapePolygon::ShapePolygon(const QPoint &c, QPolygon p) : Shape(c), m_polygon(std::move(p))
{
    m_polygon.translate(c);
}

bool ShapePolygon::isInside(const QPoint &p)
{
    return m_polygon.containsPoint(p, Qt::OddEvenFill);
}

QRect ShapePolygon::getBounds()
{
    return m_polygon.boundingRect();
}

void ShapePolygon::updateCenter()
{
    // We define the center as the first point of the polygon
    if (!m_polygon.empty())
    {
        m_center = m_polygon.first();
    }
}

void ShapePolygon::setCenter(const QPoint &c)
{
    m_polygon.translate(c - m_center);
    m_center = c;
}

void ShapePolygon::setPolygon(const QPolygon &p)
{
    m_polygon = p;
}

const QPolygon& ShapePolygon::getPolygon() const
{
    return m_polygon;
}

void ShapePolygon::draw(QPainter& painter, int offsetX, int offsetY)
{
    painter.drawPolygon(m_polygon.translated(offsetX, offsetY));
}
