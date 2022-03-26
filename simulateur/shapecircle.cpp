#include "ShapeCircle.h"

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
