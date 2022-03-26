#pragma once

#include <QtXml>

#include "forme.h"
#include "formecercle.h"
#include "formepolygone.h"
#include "scene.h"

class InterfaceXml
{
public:
    static Scene loadFile(const QString &filename);

private:
    static Scene loadScene(const QDomElement& elScene);
	static void loadShapes(const QDomElement &elShapes, FormeList &shapes);
    static bool loadShape(const QDomElement &elShape, FormePtr shape);
    static FormePtr loadShapeCircle(const QDomElement &elCircle);
    static FormePtr loadShapePolygon(const QDomElement &elPolygon);
    static void loadPoints(const QDomElement& elPoints, const QString& tagName, QList<int>& points);
};
