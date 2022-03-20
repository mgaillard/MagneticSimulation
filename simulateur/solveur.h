#ifndef SOLVEUR_H
#define SOLVEUR_H

#include <Eigen/Eigen>

class Solveur
{
public:
    Solveur();
    void setSize(const int lignes, const int colonnes);
    virtual void remplirMatricePrincipale(const Eigen::VectorXd &vecteurAlpha, const Eigen::VectorXd &vecteurBeta, const Eigen::VectorXd &vecteurGamma, const Eigen::VectorXd &vecteurDelta1, const Eigen::VectorXd &vecteurDelta2) = 0;
    virtual Eigen::MatrixXd calculerSolution(const Eigen::VectorXd &vecteurSecondMembre) = 0;
protected:
    int nbLignes;
    int nbColonnes;
};

#endif // SOLVEUR_H
