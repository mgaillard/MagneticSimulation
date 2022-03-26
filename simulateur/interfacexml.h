#pragma once

#include <QtXml>

#include "shape.h"
#include "scene.h"

class InterfaceXml
{
public:
    static Scene loadFile(const QString &filename);

private:
    static Scene loadScene(const QDomElement& elScene);
	static void loadShapes(const QDomElement &elShapes, ShapeList &shapes);
    static bool loadShapeCommonAttributes(const QDomElement &elShape, ShapePtr shape);
    static ShapePtr loadShapeCircle(const QDomElement &elCircle);
    static ShapePtr loadShapePolygon(const QDomElement &elPolygon);
    static void loadPoints(const QDomElement& elPoints, const QString& tagName, QList<int>& points);
};
