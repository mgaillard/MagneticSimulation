#ifndef FORMEPOLYGONE_H
#define FORMEPOLYGONE_H

#include <QPoint>
#include <QRect>
#include <QPolygon>

#include "forme.h"

class FormePolygone : public Forme
{
public:
    FormePolygone(const QPoint &c, const QPolygon &p);
    bool collision(const QPoint &p);
    QRect getLimites();
    void miseAJourCentre();
    void setCentre(const QPoint &c);
    void setPolygone(const QPolygon &p);
    QPolygon getPolygone() const;
private:
    QPolygon polygone;
};

#endif // FORMEPOLYGONE_H
