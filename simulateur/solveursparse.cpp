#include "solveursparse.h"

// TODO: Remove ASAP
#include <iostream>

SolveurSparse::SolveurSparse() : Solveur()
{
}

void SolveurSparse::remplirMatricePrincipale(const Eigen::VectorXd &vecteurAlpha, const Eigen::VectorXd &vecteurBeta, const Eigen::VectorXd &vecteurGamma, const Eigen::VectorXd &vecteurDelta1, const Eigen::VectorXd &vecteurDelta2)
{
    matricePrincipale.resize(nbLignes*nbColonnes, nbLignes*nbColonnes);
    std::vector<Eigen::Triplet<double> > tripletList;
    //on rempli la matrice principale
    for (int i = 0; i < nbLignes*nbColonnes; i++) {
        tripletList.push_back(Eigen::Triplet<double>(i, i, vecteurGamma(i)));
        if (i >= 1) {
            tripletList.push_back(Eigen::Triplet<double>(i, i - 1, vecteurBeta(i - 1)));
        }
        if (i >= nbColonnes) {
            tripletList.push_back(Eigen::Triplet<double>(i, i - nbColonnes, vecteurDelta1(i - nbColonnes)));
        }
        if (i <= nbLignes*nbColonnes - 2) {
            tripletList.push_back(Eigen::Triplet<double>(i, i +1, vecteurAlpha(i + 1)));
        }
        if (i <= nbLignes*nbColonnes - 1 - nbColonnes) {
            tripletList.push_back(Eigen::Triplet<double>(i, i + nbColonnes, vecteurDelta2(i + nbColonnes)));
        }
    }
    matricePrincipale.setFromTriplets(tripletList.begin(), tripletList.end());
    matricePrincipale.makeCompressed();
}
