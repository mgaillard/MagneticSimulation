#include "scene.h"

Scene::Scene() :
	m_rectangleScene(0, 0, -1, -1),
	m_step(0.0)
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

void Scene::setStep(double s)
{
    m_step = s;
}

const double& Scene::getStep() const
{
    return m_step;
}

double Scene::getSqStep() const
{
    return m_step * m_step;
}

void Scene::setShapes(const ShapeList& shapes)
{
    m_shapes = shapes;
}

const ShapeList& Scene::getShapes() const
{
    return m_shapes;
}
