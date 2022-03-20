#ifndef SOLVEURSPARSE_H
#define SOLVEURSPARSE_H

#include <Eigen/Sparse>
#include <vector>

#include "solveur.h"

class SolveurSparse : public Solveur
{
public:
    SolveurSparse();
    void remplirMatricePrincipale(const Eigen::VectorXd &vecteurAlpha, const Eigen::VectorXd &vecteurBeta, const Eigen::VectorXd &vecteurGamma, const Eigen::VectorXd &vecteurDelta1, const Eigen::VectorXd &vecteurDelta2);
protected:
    Eigen::SparseMatrix<double> matricePrincipale;
};

#endif // SOLVEURSPARSE_H
