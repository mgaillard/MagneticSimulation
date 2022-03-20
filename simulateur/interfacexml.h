#ifndef INTERFACEXML_H
#define INTERFACEXML_H

#include <QtXml>
#include <QList>
#include <algorithm>
#include <iostream>

#include "forme.h"
#include "formecercle.h"
#include "formepolygone.h"
#include "scene.h"

class InterfaceXml
{
public:
    InterfaceXml();
    void chargerFichier(const QString &fichier, Scene &scene, QList<Forme*> &formes);
    void chargerFormes(QDomElement &elFormes, QList<Forme*> &formes);
    Scene chargerScene(QDomElement &elScene);
    bool chargerForme(QDomElement &elForme, Forme* forme);
    bool chargerFormeCercle(QDomElement &elCercle, FormeCercle* cercle);
    bool chargerFormePolygone(QDomElement &elPolygone, FormePolygone* polygone);
    void chargerPoints(QDomElement &elPoints, const QString &nomBalise, QList<int> &points);
};

#endif // INTERFACEXML_H
