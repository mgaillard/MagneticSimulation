#ifndef RENDUIMAGE_H
#define RENDUIMAGE_H

#include <QFile>
#include <QImage>
#include <QColor>
#include <Eigen/Eigen>

#include <iostream>

class RenduImage
{
public:
    RenduImage();
    static void enregistrerMatrice(const QString &fichier, Eigen::MatrixXd &matrice);
    static QColor getColorByIntensity(const float intensity);
};

#endif // RENDUIMAGE_H
