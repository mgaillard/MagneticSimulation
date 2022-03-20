#ifndef FORMECERCLE_H
#define FORMECERCLE_H

#include <QPoint>
#include <QRect>
#include <math.h>

#include "forme.h"

class FormeCercle : public Forme
{
public:
    FormeCercle(const QPoint &c, const float r);
    ~FormeCercle();
    bool collision(const QPoint &p);
    QRect getLimites();
    void setRayon(const float r);
    float getRayon() const;
private:
    float rayon;
};

#endif // FORMECERCLE_H
