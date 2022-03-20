#ifndef SCENE_H
#define SCENE_H

#include <QRect>

class Scene
{
public:
    Scene();
    void setRectangleScene(const QRect &r);
    QRect getRectangleScene() const;
    void setPas(const double& p);
    const double& getPas() const;
    double getSqPas() const;
private:
    QRect rectangleScene;
    double pas;
};

#endif // SCENE_H
