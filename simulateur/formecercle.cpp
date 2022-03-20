#include "formecercle.h"

FormeCercle::FormeCercle(const QPoint &c, const float r) : Forme(c), rayon(r)
{
}

FormeCercle::~FormeCercle()
{
}

bool FormeCercle::collision(const QPoint &p)
{
    float distanceSq = (centre.x() - p.x())*(centre.x() - p.x()) + (centre.y() - p.y())*(centre.y() - p.y());
    return (distanceSq < rayon*rayon);
}

QRect FormeCercle::getLimites()
{
    return QRect(centre.x() - (int)floor(rayon), centre.y() - (int)floor(rayon), 2*(int)ceil(rayon), 2*(int)ceil(rayon));
}

void FormeCercle::setRayon(const float r)
{
    rayon = r;
}

float FormeCercle::getRayon() const
{
    return rayon;
}
