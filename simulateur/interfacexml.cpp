#include "interfacexml.h"

InterfaceXml::InterfaceXml()
{
}

void InterfaceXml::chargerFichier(const QString &fichier, Scene &scene, QList<Forme*> &formes)
{
    //on ouvre le fichier xml
    QDomDocument doc;
    QFile xml_doc(fichier);
    if(!xml_doc.open(QIODevice::ReadOnly))
    {
        std::cout << "Le document XML n'a pas pu être ouvert. Vérifiez que le nom est le bon et que le document est bien placé" << std::endl;
        return;
    }
    if (!doc.setContent(&xml_doc))
    {
        xml_doc.close();
        std::cout << "Le document XML n'a pas pu être attribué à l'objet QDomDocument." << std::endl;
        return;
    }

    QDomElement elScene = doc.documentElement();
    //on charge la scene
    scene = chargerScene(elScene);
    QDomElement child = elScene.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "formes")
        {
            chargerFormes(child, formes);
        }
        child = child.nextSiblingElement();
    }
    xml_doc.close();
}


void InterfaceXml::chargerFormes(QDomElement &elFormes, QList<Forme*> &formes)
{
    QDomElement child = elFormes.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "forme")
        {
            if (child.attribute("type") == "cercle")
            {
                FormeCercle* cercle = new FormeCercle(QPoint(-1, -1), -1);
                bool chargement = chargerFormeCercle(child, cercle);
                if (chargement)//si le chargement s'est passe sans erreur
                {
                    formes.append(cercle);//on ajoute le cercle aux formes
                }
            }
            else if (child.attribute("type") == "polygone")
            {
                FormePolygone* polygone = new FormePolygone(QPoint(-1, -1), QPolygon());
                bool chargement = chargerFormePolygone(child, polygone);
                if (chargement)//si le chargement s'est passe sans erreur
                {
                    formes.append(polygone);//on ajoute le polyone aux formes
                }
            }
        }
        child = child.nextSiblingElement();
    }
}

Scene InterfaceXml::chargerScene(QDomElement &elScene)
{
    Scene scene;
    float pas = -1;
    int hauteur = -1;
    int largeur = -1;
    int x = -1;
    int y = -1;
    //on lit les informations dans le xml
    QDomElement child = elScene.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "pas")
        {
            pas = child.text().toDouble();
        }
        else if (child.tagName() == "hauteurscene")
        {
            hauteur = child.text().toInt();
        }
        else if (child.tagName() == "largeurscene")
        {
            largeur = child.text().toInt();
        }
        else if (child.tagName() == "xscene")
        {
            x = child.text().toInt();
        }
        else if (child.tagName() == "yscene")
        {
            y = child.text().toInt();
        }
        child = child.nextSiblingElement();
    }
    //on verifie que toutes les informations sont correctes et on retourne une Scene
    if (pas >= 0 && hauteur >= 0 && largeur >= 0 && x >= 0 && y >= 0)
    {
        QRect rectangleScene(x, y, largeur, hauteur);
        scene.setRectangleScene(rectangleScene);
        scene.setPas(pas);
    }
    else
    {
        std::cout << "scene non conforme : xml errone." << std::endl;
    }
    return scene;
}

bool InterfaceXml::chargerForme(QDomElement &elForme, Forme* forme)
{
    //valeurs par defaut pour la forme
    QPoint centre(-1, -1);
    forme->setPermeabilite(-1);
    forme->setCourant(0);
    //on explore l'arbre pour completer les donnees de la forme
    QDomElement child = elForme.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "permeabilite")
        {
            forme->setPermeabilite(child.text().toDouble());
        }
        else if (child.tagName() == "courant")
        {
            forme->setCourant(child.text().toDouble());
        }
        else if (child.tagName() == "centrex")
        {
            centre.setX(child.text().toInt());
        }
        else if (child.tagName() == "centrey")
        {
            centre.setY(child.text().toInt());
        }
        child = child.nextSiblingElement();
    }

    if (centre.x() >= 0 && centre.y() >= 0 && forme->getPermeabilite() >= 0)
    {
        forme->setCentre(centre);
        //on signale que la forme a ete correctement chargee
        return true;
    }
    //il manque des donnees, on signale que la forme a mal ete chargee
    return false;
}

bool InterfaceXml::chargerFormeCercle(QDomElement &elCercle, FormeCercle* cercle)
{
    bool chargementForme = chargerForme(elCercle, cercle);//on commence par charger les attributs communs a toutes les formes.
    if (chargementForme)
    {
        //on charge les attributs propres aux cercles.
        QDomElement child = elCercle.firstChild().toElement();
        while (!child.isNull())
        {
            if (child.tagName() == "rayon")
            {
                cercle->setRayon(child.text().toDouble());
            }
            child = child.nextSiblingElement();
        }
        //si le rayon est correctement entre on signale que la forme est valide.
        if (cercle->getRayon() >= 0)
        {
            return true;
        }
    }
    return false;
}

bool InterfaceXml::chargerFormePolygone(QDomElement &elPolygone, FormePolygone* polygone)
{
    bool chargementForme = chargerForme(elPolygone, polygone);//on commence par charger les attributs communs a toutes les formes.
    if (chargementForme)
    {
        //on charge les attributs propres aux polygones.
        QList<int> xpoints;
        QList<int> ypoints;
        QDomElement child = elPolygone.firstChild().toElement();
        while (!child.isNull())
        {
            if (child.tagName() == "xpoints")
            {
                chargerPoints(child, QString("x"), xpoints);
            }
            else if(child.tagName() == "ypoints")
            {
                chargerPoints(child, QString("y"), ypoints);
            }
            child = child.nextSiblingElement();
        }
        //on a la liste des coordonnees sur x et y, on transforme ca en points.
        if (xpoints.size() > 2 && ypoints.size() > 2)//il faut au moins trois points pour faire un polygone
        {
            QVector<QPoint> points;
            int nombrePoints = std::min(xpoints.size(), ypoints.size());
            for (int i = 0;i < nombrePoints;i++)
            {
                points.append(QPoint(xpoints.at(i), ypoints.at(i)));
            }
            polygone->setPolygone(QPolygon(points));
            polygone->miseAJourCentre();
            return true;
        }
    }
    return false;
}

void InterfaceXml::chargerPoints(QDomElement &elPoints, const QString &nomBalise, QList<int> &points)
{
    QDomElement child = elPoints.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == nomBalise)
        {
            points.append(child.text().toInt());
        }
        child = child.nextSiblingElement();
    }
}
