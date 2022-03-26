#pragma once

#include <memory>
#include <vector>

#include <QPoint>
#include <QRect>

class Forme
{
public:
    Forme(const QPoint &c);
    virtual ~Forme();
    virtual bool collision(const QPoint &p) = 0;
    virtual QRect getLimites() = 0;
    int getSurface(const QRect &rectangleScene);
    void setCentre(const QPoint &c);
    QPoint getCentre() const;
    void setPermeabilite(const double p);
    double getPermeabilite() const;
    void setCourant(const double c);
    double getCourant() const;

protected:
    QPoint centre;
private:
    double permeabilite;
    double courant;
};

using FormePtr = std::shared_ptr<Forme>;
using FormeList = std::vector<FormePtr>;
