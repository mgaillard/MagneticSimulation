#include "shape.h"

Shape::Shape(const QPoint &center) : m_center(center)
{
}

int Shape::computeSurface(const QRect &rectangleScene)
{
    int nbCase = 0;

    const QRect bounds = getBounds();

    const int startX = std::max(rectangleScene.x(), bounds.x() - 1);
    const int endX = std::min(rectangleScene.x() + rectangleScene.width(), bounds.x() + bounds.width());
    const int startY = std::max(rectangleScene.y(), bounds.y() - 1);
    const int endY = std::min(rectangleScene.y() + rectangleScene.height(), bounds.y() + bounds.height());

    for (int x = startX; x < endX; x++)
    {
	    for (int y = startY; y < endY; y++)
	    {
		    if (isInside(QPoint(x, y)))
		    {
			    nbCase++;
		    }
	    }
    }

    return nbCase;
}

void Shape::setCenter(const QPoint &center)
{
    m_center = center;
}

const QPoint& Shape::getCenter() const
{
    return m_center;
}

void Shape::setPermeability(const double p)
{
    m_permeability = p;
}

const double& Shape::getPermeability() const
{
    return m_permeability;
}

void Shape::setCurrent(const double c)
{
    m_current = c;
}

const double& Shape::getCurrent() const
{
    return m_current;
}
