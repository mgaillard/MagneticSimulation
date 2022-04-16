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

Eigen::Vector2d Scene::lowerBoundMeters() const
{
    return {
        0.0,
        0.0
    };
}

Eigen::Vector2d Scene::higherBoundMeters() const
{
    return {
        m_rectangleScene.width() * m_step,
        m_rectangleScene.height() * m_step
    };
}

bool Scene::isInside(const Eigen::Vector2d& p) const
{
    const auto low = lowerBoundMeters();
    const auto high = higherBoundMeters();

    return p.x() >= low.x() && p.x() <= high.x() && p.y() >= low.y() && p.y() <= high.y();
}

Eigen::Vector2d Scene::convertToGridCoords(const Eigen::Vector2d& p) const
{
    // Min X coordinate in meters: 0.0
    // Max X coordinate in meters: m_rectangleScene.width() * m_step

    // Min Y coordinate in meters: 0.0
    // Max Y coordinate in meters: m_rectangleScene.height() * m_step

    return p / m_step;
}

Eigen::Vector2d Scene::convertFromGridCoords(const Eigen::Vector2d& p) const
{
    // Min X coordinate in pixels: 0.0
    // Max X coordinate in pixels: m_rectangleScene.width()

    // Min Y coordinate in pixels: 0.0
    // Max Y coordinate in pixels: m_rectangleScene.height()

    return p * m_step;
}

void Scene::setShapes(const ShapeList& shapes)
{
    m_shapes = shapes;
}

const ShapeList& Scene::getShapes() const
{
    return m_shapes;
}
