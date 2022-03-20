#ifndef SOLVEURDENSE_H
#define SOLVEURDENSE_H

#include <Eigen/Eigen>
#include <iostream>
#include "solveur.h"

class SolveurDense : public Solveur
{
public:
    SolveurDense();
    void remplirMatricePrincipale(const Eigen::VectorXd &vecteurAlpha, const Eigen::VectorXd &vecteurBeta, const Eigen::VectorXd &vecteurGamma, const Eigen::VectorXd &vecteurDelta1, const Eigen::VectorXd &vecteurDelta2);
protected:
    Eigen::MatrixXd matricePrincipale;
};

#endif // SOLVEURDENSE_H
