#include "forme.h"

Forme::Forme(const QPoint &c) : centre(c)
{
}

Forme::~Forme()
{

}

int Forme::getSurface(const QRect &rectangleScene)
{
    int nbCase = 0;
    QRect limites = getLimites();
    int debutX = std::max(rectangleScene.x(), limites.x() - 1);
    int finX = std::min(rectangleScene.x() + rectangleScene.width(), limites.x() + limites.width());
    int debutY = std::max(rectangleScene.y(), limites.y() - 1);
    int finY = std::min(rectangleScene.y() + rectangleScene.height(), limites.y() + limites.height());
    for (int x = debutX;x < finX;x++)
    {
        for (int y = debutY;y < finY;y++)
        {
            if (collision(QPoint(x, y)))
            {
                nbCase++;
            }
        }
    }
    return nbCase;
}

void Forme::setCentre(const QPoint &c)
{
    centre = c;
}

QPoint Forme::getCentre() const
{
    return centre;
}

void Forme::setPermeabilite(const double p)
{
    permeabilite = p;
}

double Forme::getPermeabilite() const
{
    return permeabilite;
}

void Forme::setCourant(const double c)
{
    courant = c;
}

double Forme::getCourant() const
{
    return courant;
}
