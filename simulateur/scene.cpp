#include "scene.h"

Scene::Scene()
{
}

void Scene::setRectangleScene(const QRect &r)
{
    rectangleScene = r;
}

QRect Scene::getRectangleScene() const
{
    return rectangleScene;
}

void Scene::setPas(const double& p)
{
    pas = p;
}

const double& Scene::getPas() const
{
    return pas;
}

double Scene::getSqPas() const
{
    return pas * pas;
}
