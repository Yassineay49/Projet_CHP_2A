#include <iostream>
#include <fstream>
#include <string>
#include "read.h"


int main() {

    int Nx,kmax;
    int Ny;
    double dx;
    double dy;
    double vx;
    double vy;
    double dt;
    int space_scheme;
    int time_scheme;
    double x, y, xmin, xmax, CFL, ymin, ymax;
    int cas, Nmax;
    double t, Tf;

    int tag;
    read(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny, CFL, space_scheme, time_scheme);

    dx = (xmax - xmin) / (Nx);
    dy = (ymax - ymin) / (Ny);


    dt = 1/((1/dx)+(1/dy));
    kmax=(int(Tf/dt));


    int max_i ;
    std::cout << "Entrez le nombre de proc utilisé : ";
    std::cin >> max_i;


    for (int k =0 ;k<kmax ; ++k){

    

    // Nom du fichier de sortie
    std::string outputFileName =  "sol." + std::to_string(k) + ".dat";

    // Ouvre le fichier de sortie en mode d'écriture
    std::ofstream outputFile(outputFileName, std::ios::out | std::ios::binary);

    if (!outputFile.is_open()) {
        std::cerr << "Impossible d'ouvrir le fichier de sortie." << std::endl;
        return 1;
    }


    for (int i = 0; i <= max_i-1; ++i) {
        // Nom du fichier courant
        std::string fileName = "sol." + std::to_string(k) + ".00" + std::to_string(i) +".dat";

        // Ouvre le fichier courant en mode binaire
        std::ifstream inputFile(fileName, std::ios::binary);

        if (!inputFile.is_open()) {
            // std::cerr << "Impossible d'ouvrir le fichier " << fileName << std::endl;
            break;  
        }

        // Lit le contenu du fichier courant et écrit dans le fichier de sortie
        outputFile << inputFile.rdbuf();


        inputFile.close();
    }

    outputFile.close();

    }

    std::cout << "Concaténation des fichiers terminée." << std::endl;

    return 0;
}
