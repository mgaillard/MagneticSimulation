#ifndef SIMULATEUR_H
#define SIMULATEUR_H

#include <QString>
#include <Eigen/Eigen>
#include <QList>
#include <algorithm>
#include <math.h>
#include <fstream>

#include "scene.h"
#include "forme.h"
#include "formecercle.h"
#include "interfacexml.h"
#include "renduimage.h"
#include "solveursparselu.h"
#include "solveurcolpivhouseholderqr.h"

class Simulateur
{
public:
    Simulateur(const QString fichier);
    ~Simulateur();
    QList<Forme*>& getFormes();
    void initialiser();
    bool validationSimulation();
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
    void enregistrerResultats(const QString &fichierMatriceA, const QString &fichierMatriceBr, const QString &fichierMatriceBz, const QString &fichierMatriceB);
private:
    Scene scene;
    QList<Forme*> formes;
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
    SolveurSparseLU solveur;
};

#endif // SIMULATEUR_H
