#include "solveurdense.h"

SolveurDense::SolveurDense() : Solveur()
{

}

void SolveurDense::remplirMatricePrincipale(const Eigen::VectorXd &vecteurAlpha, const Eigen::VectorXd &vecteurBeta, const Eigen::VectorXd &vecteurGamma, const Eigen::VectorXd &vecteurDelta1, const Eigen::VectorXd &vecteurDelta2)
{
    matricePrincipale = Eigen::MatrixXd::Constant(nbLignes*nbColonnes, nbLignes*nbColonnes, 0);
    //on rempli la matrice principale
    for (int i = 0; i < nbLignes*nbColonnes; i++) {
        matricePrincipale(i, i) = vecteurGamma(i);
        if (i >= 1) {
            matricePrincipale(i, i - 1) = vecteurBeta(i - 1);
        }
        if (i >= nbColonnes) {
            matricePrincipale(i, i - nbColonnes) = vecteurDelta1(i - nbColonnes);
        }
        if (i <= nbLignes*nbColonnes - 2) {
            matricePrincipale(i, i + 1) = vecteurAlpha(i + 1);

        }
        if (i <= nbLignes*nbColonnes - 1 - nbColonnes) {
            matricePrincipale(i, i + nbColonnes) = vecteurDelta2(i + nbColonnes);
        }
    }
}
