#include "interfacexml.h"

#include <iostream>

Scene InterfaceXml::loadFile(const QString& filename)
{
	QDomDocument doc;
	QFile xmlDoc(filename);

	if (!xmlDoc.open(QIODevice::ReadOnly))
	{
		std::cerr << "The XML document could not be opened." << std::endl;
		return {};
	}
	if (!doc.setContent(&xmlDoc))
	{
		xmlDoc.close();
		std::cerr << "The XML document could not be attributed to the QDomDocument object." << std::endl;
		return {};
	}

	// Loading the scene
	const QDomElement elScene = doc.documentElement();
	Scene scene = loadScene(elScene);
    FormeList shapes;

	// Loading shapes in the scene
	QDomElement child = elScene.firstChild().toElement();
	while (!child.isNull())
	{
		if (child.tagName() == "formes")
		{
			loadShapes(child, shapes);
		}
		child = child.nextSiblingElement();
	}
	xmlDoc.close();

    scene.setShapes(shapes);

    return scene;
}

Scene InterfaceXml::loadScene(const QDomElement& elScene)
{
    Scene scene;

    double pas = 0.0;
    int height = -1;
    int width = -1;
    int x = -1;
    int y = -1;

    // Read information from XML file
    QDomElement child = elScene.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "pas")
        {
            pas = child.text().toDouble();
        }
        else if (child.tagName() == "hauteurscene")
        {
            height = child.text().toInt();
        }
        else if (child.tagName() == "largeurscene")
        {
            width = child.text().toInt();
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

    // Check that information is valid
    if (pas >= 0 && height >= 0 && width >= 0 && x >= 0 && y >= 0)
    {
        const QRect rectangleScene(x, y, width, height);
        scene.setRectangleScene(rectangleScene);
        scene.setPas(pas);
    }
    else
    {
        std::cerr << "XML error: The scene is not valid" << std::endl;
    }

    return scene;
}

void InterfaceXml::loadShapes(const QDomElement &elShapes, FormeList &shapes)
{
    QDomElement child = elShapes.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "forme")
        {
            if (child.attribute("type") == "cercle")
            {
                const auto circle = loadShapeCircle(child);

            	if (circle)
                {
                    // If successfully loaded, add the circle to shapes
                    shapes.push_back(std::move(circle));
                }
            }
            else if (child.attribute("type") == "polygone")
            {
                const auto polygon = loadShapePolygon(child);

            	if (polygon)
                {
                    // If successfully loaded, add the polygon to shapes
                    shapes.push_back(std::move(polygon));
                }
            }
        }
        child = child.nextSiblingElement();
    }
}

bool InterfaceXml::loadShape(const QDomElement &elShape, FormePtr shape)
{
    // Default values for a shape
    QPoint centre(-1, -1);
    shape->setPermeabilite(-1);
    shape->setCourant(0);

    // Read actual values from the XML
    QDomElement child = elShape.firstChild().toElement();
    while (!child.isNull())
    {
        if (child.tagName() == "permeabilite")
        {
            shape->setPermeabilite(child.text().toDouble());
        }
        else if (child.tagName() == "courant")
        {
            shape->setCourant(child.text().toDouble());
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

    if (centre.x() >= 0 && centre.y() >= 0 && shape->getPermeabilite() >= 0)
    {
        shape->setCentre(centre);

    	// Shape was successfully loaded
        return true;
    }

    // Shape could not be loaded properly
    return false;
}

FormePtr InterfaceXml::loadShapeCircle(const QDomElement& elCircle)
{
    auto circle = std::make_shared<FormeCercle>(QPoint(-1, -1), -1);

	// Load attributes common to all types of shapes
	const bool res = loadShape(elCircle, circle);

	if (res)
	{
		// Loading circle attributes
		QDomElement child = elCircle.firstChild().toElement();
		while (!child.isNull())
		{
			if (child.tagName() == "rayon")
			{
				circle->setRayon(child.text().toDouble());
			}
			child = child.nextSiblingElement();
		}

		// Check the radius
		if (circle->getRayon() >= 0)
		{
			return circle;
		}
	}

	return nullptr;
}

FormePtr InterfaceXml::loadShapePolygon(const QDomElement& elPolygon)
{
    auto polygon = std::make_shared<FormePolygone>(QPoint(-1, -1), QPolygon());

	// Load attributes common to all types of shapes
	const bool res = loadShape(elPolygon, polygon);

	if (res)
	{
		// Load polygon attributes
		QList<int> pointsX;
		QList<int> pointsY;

		QDomElement child = elPolygon.firstChild().toElement();
		while (!child.isNull())
		{
			if (child.tagName() == "xpoints")
			{
				loadPoints(child, "x", pointsX);
			}
			else if (child.tagName() == "ypoints")
			{
				loadPoints(child, "y", pointsY);
			}
			child = child.nextSiblingElement();
		}

		// Convert the two lists into one list of 2D points
		// At least two points are needed for a polygon
		if (pointsX.size() > 2 && pointsY.size() > 2)
		{
			QVector<QPoint> points;

			const int nbPoints = std::min(pointsX.size(), pointsY.size());

			for (int i = 0; i < nbPoints; i++)
			{
				points.append(QPoint(pointsX.at(i), pointsY.at(i)));
			}

			polygon->setPolygone(QPolygon(points));
			polygon->miseAJourCentre();

			return polygon;
		}
	}

    return nullptr;
}

void InterfaceXml::loadPoints(const QDomElement& elPoints, const QString& tagName, QList<int>& points)
{
	QDomElement child = elPoints.firstChild().toElement();

	while (!child.isNull())
	{
		if (child.tagName() == tagName)
		{
			points.append(child.text().toInt());
		}
		child = child.nextSiblingElement();
	}
}
