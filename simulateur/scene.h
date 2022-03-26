#ifndef SCENE_H
#define SCENE_H

#include <QRect>

class Scene
{
public:
    Scene();
    void setRectangleScene(const QRect &r);
    const QRect& getRectangleScene() const;
    int resolutionHeight() const;
    int resolutionWidth() const;
    int resolutionTotal() const;
    void setPas(const double& p);
    const double& getPas() const;
    double getSqPas() const;

private:
    QRect m_rectangleScene;
    double m_pas;
};

#endif // SCENE_H
