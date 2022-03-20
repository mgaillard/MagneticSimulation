#include "formepolygone.h"

FormePolygone::FormePolygone(const QPoint &c, const QPolygon &p) : Forme(c), polygone(p)
{
    polygone.translate(c);
}

bool FormePolygone::collision(const QPoint &p)
{
    return polygone.containsPoint(p, Qt::OddEvenFill);
}

QRect FormePolygone::getLimites()
{
    return polygone.boundingRect();
}

void FormePolygone::miseAJourCentre()
{
    //on defini le centre comme le premier point
    if (polygone.size() > 0)
    {
        centre = polygone.first();
    }
}

void FormePolygone::setCentre(const QPoint &c)
{
    polygone.translate(c - centre);
    centre = c;
}

void FormePolygone::setPolygone(const QPolygon &p)
{
    polygone = p;
}

QPolygon FormePolygone::getPolygone() const
{
    return polygone;
}
