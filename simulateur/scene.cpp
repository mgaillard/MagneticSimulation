#include "scene.h"

Scene::Scene() :
	m_rectangleScene(0, 0, -1, -1),
	m_pas(0.0)
{
}

void Scene::setRectangleScene(const QRect &r)
{
    m_rectangleScene = r;
}

const QRect& Scene::getRectangleScene() const
{
    return m_rectangleScene;
}

int Scene::resolutionHeight() const
{
    return m_rectangleScene.height();
}

int Scene::resolutionWidth() const
{
    return m_rectangleScene.width();
}

int Scene::resolutionTotal() const
{
    return resolutionHeight() * resolutionWidth();
}

void Scene::setPas(const double& p)
{
    m_pas = p;
}

const double& Scene::getPas() const
{
    return m_pas;
}

double Scene::getSqPas() const
{
    return m_pas * m_pas;
}

void Scene::setShapes(const FormeList& shapes)
{
    m_shapes = shapes;
}

const FormeList& Scene::getShapes() const
{
    return m_shapes;
}
