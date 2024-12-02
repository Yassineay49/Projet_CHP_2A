#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm> // Pour utiliser std::find
#include <fstream>
#include "BiCGstab.h"
#include "Matrice.h"
#include "fonctions.h"
#include "read.h"
#include <mpi.h>

using namespace std;




void writePartialSolution(int rank , int Nx , double dx , double dy , std::vector<double>& Uloc , int iEnd , int iBeg, int j) {
    // Générer le nom du fichier avec le numéro du processus
    std::string nomFichier = "sol." + std::to_string(j) +".00"+std::to_string(rank) +".dat";

    std::ofstream fichier(nomFichier);

    if (fichier.is_open()) {
        for (int s = 0; s < iEnd-iBeg+1; s++) {
                int i = s + iBeg;
                int j = i / (Nx + 1);
                int k = i % (Nx + 1);
                int m = s + Nx +1;
                double x = -1 + k * dx;
                double y = -1 + j * dy;

                fichier << x <<" " <<y<< " " << Uloc[s] << "\n";
        }
                
        fichier.close();
    }
}

double Prod_sca(vector<double> v1, vector<double> v2)
{
    double t = 0;
    double Prod_sc;
    for (int i = 0; i < v1.size(); i++) // Correction de la condition de la boucle
    {
        t += v1[i] * v2[i];
    }
    MPI_Allreduce(&t, &Prod_sc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return Prod_sc;
}

double calculateNorm(std::vector<double>& vec) {
    double squareSum = 0.0;
    double somme_total;
    for (double x : vec) {
        squareSum += x*x ;
    }
    MPI_Allreduce(&squareSum, &somme_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return std::sqrt(somme_total);
}




void charge_a(int me, int n, int Np, int &iBeg, int &iEnd)
{
    int q,r;
    q=n/Np;
    r=n-q*Np;
    if (me<r)
    {
        iBeg=me*(q+1);
        iEnd=(me+1)*(q+1)-1;
    }
    else
    {
        iBeg=me*q+r;
        iEnd=iBeg+q-1;
    }

}

//Voici Bicg pour le code parallèle il divérge dans plusieurs cas donc je l'ai commenté


// vector<double> BiCGstab(int Nx, int Ny, double dx, double dy, double vx, double vy, double dt, int space_scheme, int time_scheme, const vector<double>& b, double eProd_scailon, int Nmax , int iBeg , int iEnd , int Np, int Me , int cas)
// {
//     int size_loc = iEnd-iBeg+1 ;
//     int n=iEnd-iBeg+1;
//     vector<double> x0(n+2*(Nx+1), 0.0),x(n, 0.0) ,p2(n+2*(Nx+1),0.0),s2(n+2*(Nx+1), 0.0),r2(n+2*(Nx+1), 0.0),x0_2(n+2*(Nx+1), 0.0); 
//     vector<double> r(n), r_tilde(n), p(n), v(n), h(n), s(n),Ax_0(n),t(n);
//     double rho0, alpha, omega, rho1, beta;
//     std::vector<double> vecr1(Nx+1,0.0);
//     std::vector<double> vecr2(Nx+1,0.0); 
//     int tag=5;
//     double prod_sc1 , prod_sc1tot;
//     double prod_sc2 , prod_sc2tot;  
  

//     // Calculer r0 = b - Ax0 avec x0 = 0

    

//     Ax_0 = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, x0, iBeg, iEnd, Np,cas);



//     if (Me != 0 && Me != Np-1) {
//             MPI_Send(&Ax_0[0], Nx+1, MPI_DOUBLE, Me-1, tag,MPI_COMM_WORLD);
//             MPI_Send(&Ax_0[size_loc-Nx-1], Nx+1, MPI_DOUBLE, Me+1, tag,MPI_COMM_WORLD);

//             MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Me-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, Me+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
//         }
//         else if (Me == 0) {
//             MPI_Send(&Ax_0[0], Nx+1, MPI_DOUBLE, Np-1, tag,MPI_COMM_WORLD);
//             MPI_Send(&Ax_0[size_loc-Nx-1], Nx+1, MPI_DOUBLE, Me+1, tag,MPI_COMM_WORLD);

//             MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Np-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                           
//             MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, Me+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            
//         }
        
//         else if (Me == Np-1) {
//             MPI_Send(&Ax_0[0], Nx+1, MPI_DOUBLE, Me-1, tag,MPI_COMM_WORLD);
//             MPI_Send(&Ax_0[size_loc-Nx-1], Nx+1, MPI_DOUBLE, 0, tag,MPI_COMM_WORLD);

//             MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Me-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                            
//             MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         }

//     for (int i = 0; i < n; i++) {
//                 x0_2[i+Nx+1]=Ax_0[i];
//         }
    
    
    
//     for (int i = 0; i < Nx+1; i++) {
//                 x0[i] = vecr1[i];
//                 x0[iEnd-iBeg+Nx+i+2]=vecr2[i];
//         }
//     for ( int i=0; i<size_loc;i++){
//                 x0[i+Nx+1] = p[i];
//         }
        



//     for (int i = 0; i < n +2*(Nx+1); ++i) {
//         r2[i] = b[i] - x0[i];
//         std :: cout << b[i];
//     }

//     r_tilde = r; // Choisir r_tilde0 = r0
//     rho0 = Prod_sca(r_tilde, r);
//     p = r;

//     for (int k = 0; k < Nmax; k++) {
//         // Calculer v = A * p
//         // std::cout << 2  ;

//         v = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, p2, iBeg, iEnd,Np,cas);
//         alpha = rho0 / Prod_sca(r_tilde, v);

//         // Calculer h = x + alpha * p
//         for (int i = 0; i < n; ++i) {
//             h[i] = x[i] + alpha * p[i];
//             if (k==0){
//                 // std::cout << v[i]  ;

//             }

//         }

//         // Calculer s = r - alpha * v
//         for (int i = 0; i < n; ++i) {
//             s[i] = r[i] - alpha * v[i];
//         }



//         if (Me != 0 && Me != Np-1) {
//             MPI_Send(&s[0], Nx+1, MPI_DOUBLE, Me-1, tag,MPI_COMM_WORLD);
//             MPI_Send(&s[size_loc-Nx-1], Nx+1, MPI_DOUBLE, Me+1, tag,MPI_COMM_WORLD);

//             MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Me-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, Me+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
//         }
//         else if (Me == 0) {
//             MPI_Send(&s[0], Nx+1, MPI_DOUBLE, Np-1, tag,MPI_COMM_WORLD);
//             MPI_Send(&s[size_loc-Nx-1], Nx+1, MPI_DOUBLE, Me+1, tag,MPI_COMM_WORLD);

//             MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Np-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                           
//             MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, Me+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            
//         }
        
//         else if (Me == Np-1) {
//             MPI_Send(&s[0], Nx+1, MPI_DOUBLE, Me-1, tag,MPI_COMM_WORLD);
//             MPI_Send(&s[size_loc-Nx-1], Nx+1, MPI_DOUBLE, 0, tag,MPI_COMM_WORLD);

//             MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Me-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                            
//             MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         }

//         for (int i = 0; i < Nx+1; i++) {
//                 s2[i] = vecr1[i];
//                 s2[iEnd-iBeg+Nx+i+1]=vecr2[i];
//         }
//         for ( int i=0; i<size_loc;i++){
//                 s2[i+Nx+1] = p[i];
//         }
        



//         prod_sc1tot=Prod_sca(s, s);



//         // Vérifier si h est une solution assez précise
//         if (sqrt(prod_sc1tot) < eProd_scailon) {
//             x = h;
//             break; // Sortir de la boucle si la solution est assez précise
//         }

//         // Calculer t = A * s

//         t = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, s2, iBeg, iEnd,Np,cas);

        

//         // On calcule omega par all_reduce

//         omega = Prod_sca(t, s) / Prod_sca(t, t);
//         // std::cout << omega<< "\n"  ;

        


//         // Calculer x = h + omega * s
//         for (int i = 0; i < n; ++i) {
//             x[i] = h[i] + omega * s[i];
//         }

//         // Calculer r = s - omega * t
//         for (int i = 0; i < n; ++i) {
//             r[i] = s[i] - omega * t[i];
//             // std::cout << x[i]<< "\n"  ;

//         }

//         // Vérifier si x est une solution assez précise
//         if (sqrt(Prod_sca(r, r)) < eProd_scailon) {
//             break; // Sortir de la boucle si la solution est assez précise
//         }

//         rho1 = Prod_sca(r_tilde, r);
//         beta = (rho1 / rho0) * (alpha / omega);
//         // std::cout << beta<< "\n"  ;


//         // Calculer p = r + beta * (p - omega * v)
//         for (int i = 0; i < n; ++i) {
//             p[i] = r[i] + beta * (p[i] - omega * v[i]);
//         }


//         rho0 = rho1;
//     }

    
//     return x;
// }


//Voici Bicg pour le code séquentiel pour le tester veuillez prendre Np=1
vector<double> BiCGstab(int Nx, int Ny, double dx, double dy, double vx, double vy, double dt, int space_scheme, int time_scheme, const vector<double>& b, double eProd_scailon, int Nmax , int iBeg , int iEnd , int Np, int Me , int cas)
{
    int n = b.size();
    vector<double> x(n, 0.0); // Initialiser x avec des zéros
    vector<double> r(n), r_tilde(n), p(n), v(n), h(n), s(n), t(n);
    double rho0, alpha, omega, rho1, beta;

    // Calculer r0 = b - Ax0 avec x0 = 0
    vector<double> Ax_0 = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, x,iBeg,iEnd,Np,cas);
    for (int i = 0; i < n; ++i) {
        r[i] = b[i] - Ax_0[i];
    }

    r_tilde = r; // Choisir r_tilde0 = r0
    rho0 = Prod_sca(r_tilde, r);
    p = r;

    for (int k = 0; k < Nmax; k++) {
        // Calculer v = A * p
        v = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, p,iBeg , iEnd,Np,cas);
        alpha = rho0 / Prod_sca(r_tilde, v);

        // Calculer h = x + alpha * p
        for (int i = 0; i < n; ++i) {
            h[i] = x[i] + alpha * p[i];
        }

        // Calculer s = r - alpha * v
        for (int i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // Vérifier si h est une solution assez précise
        if (sqrt(Prod_sca(s, s)) < eProd_scailon) {
            x = h;
            break; // Sortir de la boucle si la solution est assez précise
        }

        // Calculer t = A * s
        t = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, s, iBeg ,iEnd,Np,cas);
        omega = Prod_sca(t, s) / Prod_sca(t, t);

        // Calculer x = h + omega * s
        for (int i = 0; i < n; ++i) {
            x[i] = h[i] + omega * s[i];
        }

        // Calculer r = s - omega * t
        for (int i = 0; i < n; ++i) {
            r[i] = s[i] - omega * t[i];
        }

        // Vérifier si x est une solution assez précise
        if (sqrt(Prod_sca(r, r)) < eProd_scailon) {
            break; // Sortir de la boucle si la solution est assez précise
        }

        rho1 = Prod_sca(r_tilde, r);
        beta = (rho1 / rho0) * (alpha / omega);

        // Calculer p = r + beta * (p - omega * v)
        for (int i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        rho0 = rho1;
    }

    return x;
}

int main(int argc, char** argv) {
    int Nx;
    int Ny;
    double dx;
    double dy;
    double vx;
    double vy;
    double dt;
    int space_scheme;
    int time_scheme;
    double x, y, xmin, xmax, CFL, ymin, ymax, eProd_scailon;
    int pas, cas, Nmax;
    double t, Tf;
    int indice_bas, indice_haut;

    int tag;
    read(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny, CFL, space_scheme, time_scheme);

    dx = (xmax - xmin) / (Nx);
    dy = (ymax - ymin) / (Ny);


    tag = 5;
    pas=0;

    std::vector<double> u0_total((Nx+1)*(Ny+1), 0.0); // Initialisation de u avec des zéros

    for (int j = 0; j < Ny + 1; ++j) {
        y = ymin + j * dy;
        for (int i = 0; i < Nx + 1; ++i) {
            x = xmin + i * dx;
            u0_total[pas] = u1(x, y, t, xmin, ymin,cas);
            pas = pas + 1;
        }
    }

    
    
    
    clock_t start_time = clock();

    
    
    // Initialisation MPI
    MPI_Init(&argc, &argv);
    int Me, Np;
    MPI_Comm_rank(MPI_COMM_WORLD, &Me);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    // Calcul des indices de début et de fin pour chaque processus MPI
    int iBeg, iEnd;
    int val = (Nx + 1) * (Ny + 1);
    charge_a(Me, val, Np, iBeg, iEnd);


    MPI_Status status;

    int size_loc = iEnd - iBeg + 1;


    std::vector<double> Uloc(size_loc, 0.0);
    std::vector<double> u0_loc(size_loc+2*Nx+2, 0.0);
    std::vector<double> vecr1(Nx+1,0.0);
    std::vector<double> vecr2(Nx+1,0.0);   
    std::vector<double> Uf((Nx+1)*(Ny+1),0.0);
    pas = 0;

       

    t = 0;
    int lig, col;
    int j = 0;
    double erf=0.0;
    // dt = 1/((0.5/dx)+(0.5/dy));
    dt=dy/2;
    
    int pas2;
    pas2=0;
    Nmax=1000;


    // création des vecteurs locaux, de chaque processeur
    for (int i=iBeg; i<iEnd+1; i++) {
        pas2+=1;
        Uloc[i - iBeg] = u0_total[i];

        
    }



    j=0;

    if (Np!=1){
    while (t < Tf) {


        if (Me != 0 && Me != Np-1) {
            MPI_Send(&Uloc[0], Nx+1, MPI_DOUBLE, Me-1, tag,MPI_COMM_WORLD);
            MPI_Send(&Uloc[size_loc-Nx-1], Nx+1, MPI_DOUBLE, Me+1, tag,MPI_COMM_WORLD);

            MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Me-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, Me+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        }
        else if (Me == 0) {
            MPI_Send(&Uloc[0], Nx+1, MPI_DOUBLE, Np-1, tag,MPI_COMM_WORLD);
            MPI_Send(&Uloc[size_loc-Nx-1], Nx+1, MPI_DOUBLE, Me+1, tag,MPI_COMM_WORLD);

            MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Np-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                           
            MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, Me+1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            
        }
        
        else if (Me == Np-1) {
            MPI_Send(&Uloc[0], Nx+1, MPI_DOUBLE, Me-1, tag,MPI_COMM_WORLD);
            MPI_Send(&Uloc[size_loc-Nx-1], Nx+1, MPI_DOUBLE, 0, tag,MPI_COMM_WORLD);

            MPI_Recv(&vecr1[0], Nx+1, MPI_DOUBLE, Me-1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);                            
            MPI_Recv(&vecr2[0], Nx+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

        
        for (int i = 0; i < Nx+1; i++) {
                u0_loc[i] = vecr1[i];
                u0_loc[iEnd-iBeg+Nx+i+2]=vecr2[i];
        }
        for ( int i=0; i<size_loc;i++){
                u0_loc[i+Nx+1] = Uloc[i];
        }


        
        
        
        
        
        if (time_scheme==1){
            // Uloc = BiCGstab( Nx,  Ny,  dx,  dy,  vx,  vy,  dt,  space_scheme,  time_scheme, u0_loc,  eProd_scailon,  Nmax ,  iBeg ,  iEnd,Np,Me,cas);
            std::cout << "Bicg parall ne fonctionne pas\n"  ;

        }
        else if (time_scheme==2){
            Uloc = produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, u0_loc,iBeg,iEnd,Np,cas);
        }


        writePartialSolution(Me ,  Nx ,  dx ,  dy , Uloc ,  iEnd ,  iBeg,  j) ;





            j=j+1;
            t = t + dt;
                
        }

    };

    if (Np==1){
        std::vector<double> err((Nx + 1) * (Ny + 1), 0.0); 
        std::string nomFichierer = "aerreur0.dat";
        std::ofstream fichiererreur(nomFichierer);
        while (t < Tf) {
        if (time_scheme==1){
            Uf = BiCGstab( Nx,  Ny,  dx,  dy,  vx,  vy,  dt,  space_scheme,  time_scheme, u0_total,  eProd_scailon,  Nmax ,  iBeg ,  iEnd,Np,Me,cas);

        }
        else if (time_scheme==2){
            Uf=produitmatvect(Nx, Ny, dx, dy, vx, vy, dt, space_scheme, time_scheme, u0_total,iBeg,iEnd,Np,cas);
        }    


        pas = 0;
        x = xmin;
        y = ymin;

        std::string nomFichier = "sol." + std::to_string(j) + ".dat";
        std::ofstream fichier(nomFichier);
        j = j + 1;
        for (int i = 0; i < (Nx + 1) * (Ny + 1); ++i) {
            fichier << x << " " << y << " " << u0_total[i] << endl;
            err [i]= abs(u0_total[i] - sin(2*M_PI*x-t)* sin(2*M_PI*y));
            lig += 1;
            pas += 1;
            if (pas == Nx + 1) {
                x = xmin - dx;
                y = y + dy;
                pas = 0;
            }
            
            x = x + dx;
        }
        erf = sqrt(dx*dy*calculateNorm(err));
        fichiererreur << erf << " " << t <<endl;
        t = t + dt;
        

        for (int i = 0; i < (Nx + 1) * (Ny + 1); ++i) {
            u0_total[i] = Uf[i];
        }
    }

  

    }
    



    // Finalisation MPI
    MPI_Finalize();
    clock_t endtime = clock();


    double execution_time = (double)(endtime - start_time) / CLOCKS_PER_SEC;
    std::cout << "Temps d'execution du code: " << execution_time << " secondes" << std::endl;



    return 0;
}













