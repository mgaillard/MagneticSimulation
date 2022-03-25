#include "simulateur.h"

#include <chrono>

#include "constants.h"
#include "export_vtk.h"

Simulateur::Simulateur(const QString fichier)
{
    //on charge le fichier
    InterfaceXml interface;
    interface.chargerFichier(fichier, scene, formes);
    initialiser();
}

Simulateur::~Simulateur()
{
    qDeleteAll(formes);
    formes.clear();
}


QList<Forme*>& Simulateur::getFormes()
{
    return formes;
}

void Simulateur::initialiser()
{
    //on initialise les matrices
    matricePermittivite = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), constants::mu_0);
    matriceDensiteCourant = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), 0.0);
    matriceDeriveeZ = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), 0.0);
    matriceDeriveeR = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), 0.0);
    vecteurSecondMembre = Eigen::MatrixXd::Constant(scene.getRectangleScene().height()*scene.getRectangleScene().width(), 1, 0.0);
    vecteurSolution = Eigen::MatrixXd::Constant(scene.getRectangleScene().height()*scene.getRectangleScene().width(), 1, 0.0);
    matriceBr = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), 0.0);
    matriceBz = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), 0.0);
    matriceB = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), scene.getRectangleScene().width(), 0.0);
    solveur.setSize(scene.getRectangleScene().height(), scene.getRectangleScene().width());
}

bool Simulateur::validationSimulation()
{
    //pour pouvoir simuler il faut qu'il y ait des formes
    if (formes.size() <= 0)
    {
        return false;
    }
    return true;
}

void Simulateur::simuler()
{
    initialiser();
	remplirMatricePermittivite();
    lissageMatricePermittivite();
    remplirMatriceDensiteCourant();
    calculMatricesDerivees();
    transformerMatricesEnVecteurs();
    calculVecteurSecondMembre();
    preparerSolveur();
    calculSolution();
    transformerVecteurSolutionEnMatrice();
    calculChampB();
}

void Simulateur::remplirMatricePermittivite()
{
    const int nbLigne = scene.getRectangleScene().height();
    const int nbColonne = scene.getRectangleScene().width();

    for (int i = 0; i < nbLigne; i++) {
        for (int j = 0; j < nbColonne; j++) {
            QPoint point(j, i);
            QList<Forme*>::iterator f;
            for (f = formes.begin(); f != formes.end(); ++f)
            {
                if ((*f)->collision(point)) {
                    // sur chaque coefficient, on multiplie la permeabilite de ce milieu par la permeabilite relative de la forme
                    matricePermittivite(i, j) *= (*f)->getPermeabilite();
                }
            }
        }
    }
    //on inverse les coefficients de la matrice
    for (int i = 0; i < nbLigne; i++) {
        for (int j = 0; j < nbColonne; j++) {
            matricePermittivite(i, j) = 1.0/matricePermittivite(i, j);
        }
    }
}

void Simulateur::lissageMatricePermittivite()
{
    //on rempli la matrice de filtrage gaussienne
    const int tailleFiltre = 10;
    const double alpha = 0.1;
    Eigen::MatrixXd filtre = Eigen::MatrixXd::Constant(tailleFiltre, tailleFiltre, 0);
    double sommeCoeffsFiltre = 0;
    for (int i = 0;i < tailleFiltre;i++)
    {
        for (int j = 0;j < tailleFiltre;j++)
        {
            filtre(i, j) = exp(-alpha*(pow(i-tailleFiltre/2, 2)+pow(j-tailleFiltre/2, 2)));
            sommeCoeffsFiltre += filtre(i, j);
        }
    }
    //on divise les coefficients du filtre par la somme des coefficients
    filtre /= sommeCoeffsFiltre;
    //on effectue le filtrage
    int nbLigne = scene.getRectangleScene().height();
    int nbColonne = scene.getRectangleScene().width();
    Eigen::MatrixXd nouvelleMatricePermittivite = Eigen::MatrixXd::Constant(nbLigne, nbColonne, 0);
    int hauteur = tailleFiltre/2;
    int largeur = tailleFiltre/2;
    //on fait le produit de convolution
    for (int i = hauteur; i < nbLigne - hauteur; i++) {
        for (int j = largeur; j < nbColonne - largeur; j++) {
            for (int c = 0;c < tailleFiltre;c++) {
                for (int d = 0;d < tailleFiltre;d++) {
                    nouvelleMatricePermittivite(i, j) += matricePermittivite(i - hauteur + c, j - largeur + d) * filtre(c, d);
                }
            }
        }
    }
    //on rempli ce qu'on ne pourra pas calculer
    for (int i = 0; i < hauteur; i++) {
        for (int j = 0; j < nbColonne; j++) {
            nouvelleMatricePermittivite(i, j) = nouvelleMatricePermittivite(hauteur, j);
        }
    }
    for (int i = nbLigne - hauteur; i < nbLigne ; i++) {
        for (int j = 0; j < nbColonne; j++) {
            nouvelleMatricePermittivite(i, j) = nouvelleMatricePermittivite(nbLigne - hauteur - 1, j);
        }
    }
    for (int i = 0; i < nbLigne; i++) {
        for (int j = 0; j < largeur; j++) {
            nouvelleMatricePermittivite(i, j) = nouvelleMatricePermittivite(i, largeur);
        }
    }
    for (int i = 0; i < nbLigne; i++) {
        for (int j = nbColonne - largeur; j < nbColonne; j++) {
            nouvelleMatricePermittivite(i, j) = nouvelleMatricePermittivite(i, nbColonne - largeur - 1);
        }
    }
    matricePermittivite = nouvelleMatricePermittivite;
}

void Simulateur::remplirMatriceDensiteCourant()
{
    int nbLigne = scene.getRectangleScene().height();
    int nbColonne = scene.getRectangleScene().width();
    //on genere un tableau des surfaces des formes
    QVector<double> surfaces;
    QList<Forme*>::iterator f;
    for (f = formes.begin(); f != formes.end(); ++f)
    {
        // Surface in number of cells
        const auto surfaceInCells = (*f)->getSurface(scene.getRectangleScene());
        // Surface in meters
        const auto surfaceInMeters = static_cast<double>(surfaceInCells) * scene.getPas() * scene.getPas();
        surfaces.append(surfaceInMeters);
    }
    //on remplis la matrice des densite de courant
    for (int i = 0; i < nbLigne; i++) {
        for (int j = 0; j < nbColonne; j++) {
            QPoint point(j, i);
            for (f = formes.begin(); f != formes.end(); ++f)
            {
                if ((*f)->collision(point)) {
                    // The current density in the cell is the current in ampere divided by the surface in meters
                    matriceDensiteCourant(i, j) += (*f)->getCourant()/surfaces.at((int)(f - formes.begin()));
                }
            }
        }
    }
}


void Simulateur::calculMatricesDerivees()
{
    int nbLigne = scene.getRectangleScene().height();
    int nbColonne = scene.getRectangleScene().width();
    for (int i = 1; i < nbLigne - 1; i++) {
        for (int j = 1; j < nbColonne - 1; j++) {
            matriceDeriveeZ(i, j) = (matricePermittivite(i+1, j) - matricePermittivite(i-1, j))/(2*scene.getPas());
            matriceDeriveeR(i, j) = (matricePermittivite(i, j+1) - matricePermittivite(i, j-1))/(2*scene.getPas());
        }
    }
}


void Simulateur::transformerMatricesEnVecteurs()
{
    matricePermittivite.transposeInPlace();
    matricePermittivite.resize(scene.getRectangleScene().height()*scene.getRectangleScene().width(), 1);
    matriceDensiteCourant.transposeInPlace();
    matriceDensiteCourant.resize(scene.getRectangleScene().height()*scene.getRectangleScene().width(), 1);
    matriceDeriveeZ.transposeInPlace();
    matriceDeriveeZ.resize(scene.getRectangleScene().height()*scene.getRectangleScene().width(), 1);
    matriceDeriveeR.transposeInPlace();
    matriceDeriveeR.resize(scene.getRectangleScene().height()*scene.getRectangleScene().width(), 1);
}

void Simulateur::transformerVecteurSolutionEnMatrice()
{
    vecteurSolution.resize(scene.getRectangleScene().width(), scene.getRectangleScene().height());
    vecteurSolution.transposeInPlace();
}

void Simulateur::preparerSolveur()
{
    // Resolution constants
    const int rows = scene.getRectangleScene().height();
    const int cols = scene.getRectangleScene().width();
    const int nbCells = rows * cols;

    // Dimension constants: h: step size on axis r, k: step size on axis z
    const double h = scene.getPas();
    const double k = scene.getPas();
    const double hSq = scene.getSqPas();
    const double kSq = scene.getSqPas();

    //on initialise les vecteurs dont nous auront besoin pour remplir la matrice principale
    Eigen::VectorXd vecAlpha = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecBeta = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecGamma = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecDelta1 = Eigen::VectorXd::Constant(nbCells, 0);
    Eigen::VectorXd vecDelta2 = Eigen::VectorXd::Constant(nbCells, 0);

    for (int i = 0; i < rows; i++)
    {
	    for (int j = 1; j < cols - 1; j++)
	    {
		    const int n = cols * i + j;

		    vecAlpha(n + 1) = ((1.0 + 1.0 / (2 * j)) / hSq) + matriceDeriveeR(n, 0) / (2.0 * h) / matricePermittivite(n, 0);
		    vecBeta(n - 1) = ((1.0 - 1.0 / (2 * j)) / hSq) - matriceDeriveeR(n, 0) / (2.0 * h) / matricePermittivite(n, 0);
		    vecGamma(n) = (-2.0 / hSq - 2.0 / kSq - 1.0 / (j * j * hSq)) + matriceDeriveeR(n, 0) / (j * h) / matricePermittivite(n, 0);

	    	if (j == 1)
		    {
			    vecAlpha(n) = 0.0;
			    vecBeta(n - 1) = 0.0;
			    vecGamma(n - 1) = 1.0;
		    }
	    }
    }

    for (int i = 1; i < rows; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    const int n = cols * i + j;

		    vecDelta2(n) = 1.0 / kSq + matriceDeriveeZ(n - cols, 0) / (2 * k) / matricePermittivite(n - cols, 0);

	    	if (j == 0)
		    {
			    vecDelta2(n) = 0.0;
		    }
	    }
    }

    for (int i = 0; i < rows - 1; i++)
    {
	    for (int j = 0; j < cols; j++)
	    {
		    const int n = cols * i + j;

		    vecDelta1(n) = 1.0 / kSq - matriceDeriveeZ(n + cols, 0) / (2 * k) / matricePermittivite(n + cols, 0);

	    	if (j == 0)
		    {
			    vecDelta1(n) = 0.0;
		    }
	    }
    }

    solveur.remplirMatricePrincipale(vecAlpha, vecBeta, vecGamma, vecDelta1, vecDelta2);
}


void Simulateur::calculVecteurSecondMembre()
{
    for (int i = 0;i < scene.getRectangleScene().height()*scene.getRectangleScene().width();i++) {
        vecteurSecondMembre(i) = -matriceDensiteCourant(i, 0) / matricePermittivite(i, 0);
    }
}

void Simulateur::calculSolution()
{
    const auto start_time = std::chrono::steady_clock::now();
    vecteurSolution = solveur.calculerSolution(vecteurSecondMembre);
    const auto end_time = std::chrono::steady_clock::now();

    std::cout << "solving time: " << std::chrono::duration<double, std::milli>(end_time - start_time).count() << " ms" << std::endl;
}

void Simulateur::calculChampB()
{
    // Resolution constants
    const int rows = scene.getRectangleScene().height();
    const int cols = scene.getRectangleScene().width();

    // Dimension constants: h: step size on axis r, k: step size on axis z
    const double h = scene.getPas();
    const double k = scene.getPas();

    //on rempli le milieu par B = rot(A).
    for (int i = 1;i < rows - 1;i++)
    {
        for (int j = 1;j < cols - 1;j++)
        {
            matriceBz(i, j) = vecteurSolution(i, j) / (h * j) + (vecteurSolution(i, j + 1) - vecteurSolution(i, j - 1)) / (2.0 * h);
            matriceBr(i, j) = (-vecteurSolution(i + 1, j) + vecteurSolution(i - 1, j)) / (2.0 * k);
        }
    }

    //on remplis les bords.
    for (int i = 0;i < rows;i++)
    {
        matriceBz(i, 0) = matriceBz(i, 1);
        matriceBr(i, 0) = matriceBr(i, 1);
        matriceBz(i, cols - 1) = matriceBz(i, cols - 2);
        matriceBr(i, cols - 1) = matriceBr(i, cols - 2);
    }

    for (int j = 0;j < cols;j++)
    {
        matriceBz(0, j) = matriceBz(1, j);
        matriceBr(0, j) = matriceBr(1, j);
        matriceBz(rows - 1, j) = matriceBz(rows - 2, j);
        matriceBr(rows - 1, j) = matriceBr(rows - 2, j);
    }

    //on calcule le champ B avec la norme des vecteurs Br et Bz.
    for (int i = 0;i < rows;i++)
    {
        for (int j = 0;j < cols;j++)
        {
            matriceB(i, j) = sqrt(pow(matriceBr(i, j), 2) + pow(matriceBz(i, j), 2));
        }
    }
}


Eigen::MatrixXd Simulateur::symetriqueMatrice(const Eigen::MatrixXd &solution, const double coefficient)
{
    Eigen::MatrixXd solutionSymetrique = Eigen::MatrixXd::Constant(scene.getRectangleScene().height(), 2*scene.getRectangleScene().width() - 1, 0);
    //dans la partie de droite on recopie la matrice solution
    for (int i = 0;i < scene.getRectangleScene().height();i++) {
        for (int j = scene.getRectangleScene().width() - 1;j < 2*scene.getRectangleScene().width() - 1;j++) {
            solutionSymetrique(i, j) = solution(i, j-scene.getRectangleScene().width()+1);
        }
    }
    //dans la partie gauche on recopie l'oppose de la matrice
    for (int i = 0;i < scene.getRectangleScene().height();i++) {
        for (int j = 0;j < scene.getRectangleScene().width() - 1;j++) {
            solutionSymetrique(i, j) = coefficient * solution(i, scene.getRectangleScene().width() - 1 - j);
        }
    }
    return solutionSymetrique;
}

void Simulateur::enregistrerResultats(const QString &fichierMatriceA, const QString &fichierMatriceBr, const QString &fichierMatriceBz, const QString &fichierMatriceB)
{
    Eigen::MatrixXd matriceAEspace = symetriqueMatrice(vecteurSolution, -1);
    Eigen::MatrixXd matriceBrEspace = symetriqueMatrice(matriceBr, 1);
    Eigen::MatrixXd matriceBzEspace = symetriqueMatrice(matriceBz, 1);
    Eigen::MatrixXd matriceBEspace = symetriqueMatrice(matriceB, 1);
    RenduImage::enregistrerMatrice(fichierMatriceA, matriceAEspace);
    RenduImage::enregistrerMatrice(fichierMatriceBr, matriceBrEspace);
    RenduImage::enregistrerMatrice(fichierMatriceBz, matriceBzEspace);
    RenduImage::enregistrerMatrice(fichierMatriceB, matriceBEspace);

    //on enregistre le vecteur
    std::ofstream fichier("matricebcentre.txt", std::ios::out | std::ios::trunc);
    if(fichier)
    {
        for (int i = 0;i < scene.getRectangleScene().height();i++)
        {
            fichier << matriceB(i, 0) << std::endl;
        }
        fichier.close();
    }

    // TODO: Transform vectors back to matrices (matricePermittivite, matriceDensiteCourant)

    // Output matrices in VTK format
    // exportScalarMatrixVtk("mu.vtk", matricePermittivite, scene.getPas(), scene.getPas());
    // exportScalarMatrixVtk("i.vtk", matriceDensiteCourant, scene.getPas(), scene.getPas());
	exportScalarMatrixVtk("a.vtk", vecteurSolution, scene.getPas(), scene.getPas());
    exportVectorMatrixVtk("b.vtk", matriceBr, matriceBz, scene.getPas(), scene.getPas());
}
