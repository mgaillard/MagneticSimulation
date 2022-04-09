#include "shapecircle.h"

#include <cmath>

ShapeCircle::ShapeCircle(const QPoint &c, double r) : Shape(c), m_radius(r)
{
}

bool ShapeCircle::isInside(const QPoint &p)
{
    const auto cx = static_cast<double>(m_center.x());
    const auto cy = static_cast<double>(m_center.y());
    const auto px = static_cast<double>(p.x());
    const auto py = static_cast<double>(p.y());

    const double distanceSq = (cx - px) * (cx - px) + (cy - py) * (cy - py);

	return distanceSq < m_radius * m_radius;
}

QRect ShapeCircle::getBounds()
{
    return {
	    m_center.x() - static_cast<int>(std::floor(m_radius)),
        m_center.y() - static_cast<int>(std::floor(m_radius)),
        2 * static_cast<int>(std::ceil(m_radius)),
        2 * static_cast<int>(std::ceil(m_radius))
    };
}

void ShapeCircle::setRadius(double r)
{
    m_radius = r;
}

const double& ShapeCircle::getRadius() const
{
    return m_radius;
}

void ShapeCircle::draw(QPainter& painter, int offsetX, int offsetY)
{
    // The rectangle starts at the center + offset - radius, and the width and height are twice the radius
    const QRectF rectangle(m_center.x() + offsetX - m_radius,
                           m_center.y() + offsetY - m_radius,
                           2.0 * m_radius,
                           2.0 * m_radius);

    // The start angle corresponds to the angle at the 3 o'clock position
    constexpr int startAngle = 0;
    // The span angle is defined in 1 / 16th of a degrees, so the full circle equals 360 * 16
    constexpr int spanAngle = 360 * 16;

    painter.drawArc(rectangle, startAngle, spanAngle);
}
