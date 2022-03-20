#include "renduimage.h"

RenduImage::RenduImage()
{
}

void RenduImage::enregistrerMatrice(const QString &fichier, Eigen::MatrixXd &matrice)
{
    //on trouve le min et le max
    double min = matrice(0, 0);
    double max = matrice(0, 0);
    for (int i = 0;i < matrice.rows();i++)
    {
        for (int j = 0;j < matrice.cols();j++)
        {
            double coeff = matrice(i, j);
            if (coeff > max)
            {
                max = coeff;
            }
            else if (coeff < min)
            {
                min = coeff;
            }
        }
    }

    //on genere une image
    QImage image(matrice.cols(), matrice.rows(), QImage::Format_RGB32);
    for (int i = 0;i < matrice.rows();i++)
    {
        for (int j = 0;j < matrice.cols();j++)
        {
            image.setPixel(QPoint(j, i), getColorByIntensity((float)((matrice(i, j) - min)/(max - min))).rgb());
        }
    }
    image.save(fichier);
}

QColor RenduImage::getColorByIntensity(const float intensity)
{
    float red = 0.f;
    float green = 0.f;
    float blue = 0.f;
    if (intensity >= 0.f && intensity < 0.125f) {
            red = 0.f;
            green = 0.f;
            blue = 4*intensity + 0.5f;
    } else if (intensity >= 0.125f && intensity < 0.375f) {
            red = 0;
            green = 4*intensity - 0.5f;
            blue = 1;
    } else if (intensity >= 0.375f && intensity < 0.625f) {
            red = 4*intensity - 1.5f;
            green = 1;
            blue = -4*intensity + 2.5f;
    } else if (intensity >= 0.625f && intensity < 0.875f) {
            red = 1;
            green = -4*intensity + 3.5f;
            blue = 0;
    } else if (intensity >= 0.875f && intensity <= 1.f) {
            red = -4*intensity + 4.5f;
            green = 0;
            blue = 0;
    }
    return QColor(red*255, green*255, blue*255);
}
