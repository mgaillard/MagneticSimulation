#ifndef SIMULATEUR_H
#define SIMULATEUR_H

#include <QString>
#include <Eigen/Eigen>
#include <QList>
#include <algorithm>
#include <math.h>
#include <fstream>

#include "scene.h"
#include "shape.h"
#include "ShapeCircle.h"
#include "interfacexml.h"
#include "solversparselu.h"
#include "solvercolpivhouseholderqr.h"

class Simulateur
{
public:
    Simulateur(const QString& filename);
    void initialiser();
    bool validationSimulation() const;
    void simuler();
    void remplirMatricePermittivite();
    void lissageMatricePermittivite();
    void remplirMatriceDensiteCourant();
    void calculMatricesDerivees();
    void transformerMatricesEnVecteurs();
    void transformerVecteurSolutionEnMatrice();
    void preparerSolveur();
    void calculVecteurSecondMembre();
    void calculSolution();
    void calculChampB();
    Eigen::MatrixXd symetriqueMatrice(const Eigen::MatrixXd &solution, const double coefficient);
    void enregistrerResultats(const QString& fichierMatriceA,
                              const QString& fichierMatriceBr,
                              const QString& fichierMatriceBz,
                              const QString& fichierMatriceB,
                              const QString& fileMatrixBWithScene);
private:
    Scene scene;
    //variables de la simulation
    Eigen::MatrixXd matricePermittivite;
    Eigen::MatrixXd matriceDensiteCourant;
    Eigen::MatrixXd matriceDeriveeZ;
    Eigen::MatrixXd matriceDeriveeR;
    Eigen::VectorXd vecteurSecondMembre;
    Eigen::MatrixXd vecteurSolution;
    Eigen::MatrixXd matriceBr;
    Eigen::MatrixXd matriceBz;
    Eigen::MatrixXd matriceB;
    SolverSparseLU solveur;
};

#endif // SIMULATEUR_H
