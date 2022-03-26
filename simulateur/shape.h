#pragma once

#include <memory>
#include <vector>

#include <QPoint>
#include <QRect>

class Shape
{
public:
    explicit Shape(const QPoint &center);
    virtual ~Shape() = default;

	virtual bool isInside(const QPoint &p) = 0;
    virtual QRect getBounds() = 0;

    int computeSurface(const QRect &rectangleScene);

	virtual void setCenter(const QPoint &center);
    const QPoint& getCenter() const;

	void setPermeability(double p);
    const double& getPermeability() const;

	void setCurrent(double c);
    const double& getCurrent() const;

protected:
    QPoint m_center;

private:
    double m_permeability;
    double m_current;

};

using ShapePtr = std::shared_ptr<Shape>;
using ShapeList = std::vector<ShapePtr>;
