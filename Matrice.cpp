#include <iostream>
#include <vector>
#include <algorithm> 
#include<fstream>
#include<mpi.h>

#include "Matrice.h"

#include "fonctions.h"

#include "read.h"


using namespace std;

// Fonction pour calculer le produit matrice-vecteur
std::vector<double> produitmatvect(int& Nx, int& Ny, double& dx, double& dy, double& xmin, double& ymin, double& dt, int& space_scheme, int& time_scheme, std::vector<double>& u , int &iBeg , int &iEnd, int &Np,int &cas) {
    int num_elements = (Nx+1) * (Ny+1);
    

   
   if (Np!=1){
    std::vector<double> produit(iEnd-iBeg+1);
    double x , y ,vx, vy ;
    if (space_scheme == 1 && time_scheme == 2) {
            for (int s = 0; s < iEnd-iBeg+1; s++) {
                int i = s + iBeg;
                int j = i / (Nx + 1);
                int k = i % (Nx + 1);
                int m = s + Nx +1;
                double x = -1 + k * dx; // nous avons pris -1 car c'est la valeur de xmin et ymin dans le cas teste 3
                double y = -1 + j * dy;
                vx = vxx(x, y, 1.0, cas);
                vy = vyy(x, y, 1.0 , cas);

                if (vy <= 0 && vx <= 0) {
                    double gamma, beta, alpha;
                    beta = -dt * (vx / dx);
                    gamma = -dt * (vy / dy);
                    alpha = 1 -  gamma -  beta;
                    if ((i+1) % (Nx+1) == 0) {
                        produit[s] = alpha * u[m] + beta * u[m-Ny-1] + gamma * u[m+Nx+1];
                        
                    }
                    else{
                        produit[s] = alpha * u[m] + beta * u[m+1] + gamma * u[m+Nx+1];

                    }

                
                }
            

                if (vy >= 0 && vx >= 0) {
                    double gamma, beta, alpha;
                    alpha = -dt * (vx / dx);
                    beta = -dt * (vy / dy);
                    gamma = 1 + alpha + beta;
            
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        

                            produit[s] = gamma * u[m] - alpha * u[m+Nx] - beta * u[m-Nx-1];
                        
                        }
                    else {
                            produit[s] = gamma * u[m] - alpha * u[m-1] - beta * u[m-Nx-1];

                    }


                }

                           
                if (vy >= 0 && vx <= 0) {
                    double gamma, beta, alpha;
                    beta = -dt * (vy / dy);
                    alpha = -dt * (vx / dx);
                    gamma = 1 - alpha +  beta;
                    if ((i+1) % (Nx+1) == 0) {
                        

                        produit[s] = - beta * u[m-Nx-1] + alpha * u[m-Ny-1] + gamma * u[m] ;
                        
                    }


                    else{
                        produit[s] = - beta * u[m-Nx-1] + alpha * u[m+1] + gamma * u[m] ;

                
                    }

                } 
 


                if (vy <= 0 && vx >= 0) {
                    double gamma, beta, alpha;
                    beta = -dt * (vy / dy);
                    alpha = -dt * (vx / dx);
                    gamma = 1 + alpha -  beta;   
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        

                        produit[s] = beta * u[m+Nx+1] -alpha * u[m+Nx]  + gamma * u[m] ;

                        
                    }                       
                    else{
                        produit[s] = beta * u[m+Nx+1] -alpha * u[m-1]  + gamma * u[m] ;

                    }


                }        
            }




    
    }

    if (space_scheme == 2 && time_scheme == 1) {
        

            for (int s = 0; s < iEnd-iBeg+1; s++) {
                int i = s + iBeg;
                int j = i / (Nx + 1);
                int k = i % (Nx + 1);
                int m = s + Nx +1;
                double x = -1 + k * dx;
                double y = -1 + j * dy;
                vx = vxx(x, y, 1.0, cas);
                vy = vyy(x, y, 1.0, cas);


                double alpha = dt * (vx / (2*dx));
                double beta = dt * (vy / (2*dy));

                if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[s] = u[m] - alpha * u[m+Nx] - beta * u[m-Nx-1] + alpha * u[m+1] + beta * u[m+Nx+1];
                        
                    }
                else if ((i+1) % (Nx+1) == 0) {
                        

                        // Calcul du produit matrice-vecteur
                        produit[s] = u[m] + alpha * u[m-Nx] + beta * u[m+Nx+1]- beta * u[m-Nx-1] - alpha*u[m-1];
                    }
                else{
                        produit[s] = u[m] -alpha * u[m-1] - beta * u[m-Nx-1] + alpha * u[m+1] + beta * u[m+Nx+1];
                    }
                
                
                
            
        }
    }
        return produit;

   }
   

   else if (Np==1){

    int num_elements = (Nx+1) * (Ny+1);
    std::vector<double> produit(num_elements);

    if (space_scheme == 1 && time_scheme == 2) {
        for (int i = 0; i < num_elements; i++) {
        double x1 = -1 + (i % (Nx +1) * dx) ;
        double y1 = -1 + (i / (Nx +1)* dy);
        double vx = vxx(x1,y1,1.0,cas);
        double vy = vyy(x1,y1,1.0,cas);


     

        
        if (vy <= 0 && vx <= 0) {
            double gamma, beta, alpha;
            beta = -dt * (vx / dx);
            gamma = -dt * (vy / dy);
            alpha = 1 -  gamma -  beta;


                if (i < (Nx+1)*Ny)
                {
                    if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = alpha * u[i] + beta * u[i-Nx] + gamma * u[i+Nx+1];
                        
                    }
                    else{
                        produit[i] = alpha * u[i] + beta * u[i+1] + gamma * u[i+Nx+1];

                    }
                }
                else{
                    if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = alpha * u[i] + beta * u[i-Nx] + gamma * u[i+Nx+1-num_elements];
                        
                    }
                    else{
                        produit[i] = alpha * u[i] + beta * u[i+1] + gamma * u[i+Nx+1-num_elements];

                    }


                
                
                
            }
        }

        if (vy >= 0 && vx >= 0) {
            double gamma, beta, alpha;
            alpha = -dt * (vx / dx);
            beta = -dt * (vy / dy);
            gamma = 1 + alpha + beta;


                if (i < Nx+1)
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        

                        // Calcul du produit matrice-vecteur
                        produit[i] = gamma * u[i] -alpha * u[i+Nx] - beta * u[i+num_elements-Nx-1];
                        
                    }
                    if ((i+1+Nx) % (Nx+1) != 0){
                        produit[i] = gamma * u[i] -alpha * u[i-1] - beta * u[i+num_elements-Nx-1];

                    }
                }
                else{
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        

                        // Calcul du produit matrice-vecteur
                        produit[i] = gamma * u[i] - alpha * u[i+Nx] - beta * u[i-Nx-1];
                        
                    }
                    if ((i+1+Nx) % (Nx+1) != 0){
                        produit[i] = gamma * u[i] - alpha * u[i-1] - beta * u[i-Nx-1];

                    }


                }
                
                
            
        }

            
            
            if (vy >= 0 && vx <= 0) {
                double gamma, beta, alpha;
                beta = -dt * (vy / dy);
                alpha = -dt * (vx / dx);
                gamma = 1 - alpha +  beta;


                if ( i < Nx+1)
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = gamma * u[i] +alpha * u[i-Nx] - beta * u[i+num_elements-Nx-1] ;

                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = gamma * u[i] + alpha * u[i+1] - beta * u[i+num_elements-Nx-1];

                    }
                    else{
                        produit[i] = gamma * u[i] + alpha * u[i+1] - beta * u[i+num_elements-Nx-1] ; 

                    }
                }

                if (i >= (Nx+1)*Ny )
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = - beta * u[i-Nx-1] + alpha * u[i-Nx] + gamma * u[i] ;
                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = - beta * u[i-Nx-1] + alpha * u[i+1] + gamma * u[i] ;
                    }
                    else{
                        produit[i] = - beta * u[i-Nx-1] + alpha * u[i+1] + gamma * u[i] ;
                    }
                }


                if (i >= Nx+1 and i < (Nx+1)*Ny ){
                    
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = - beta * u[i-Nx-1] +alpha * u[i-Nx] + gamma * u[i] ;
                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = - beta * u[i-Nx-1] + alpha * u[i+1] + gamma * u[i] ;
                    }
                    else{
                        produit[i] = - beta * u[i-Nx-1] + alpha * u[i+1] + gamma * u[i] ; 
                    }


                }
                
                
            
        }

        
    
            if (vy <= 0 && vx >= 0) {
                double gamma, beta, alpha;
                beta = -dt * (vy / dy);
                alpha = -dt * (vx / dx);
                gamma = 1 + alpha -  beta;


                if ( i < Nx+1)
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = beta * u[i+Nx+1] -alpha * u[i+Nx]  + gamma * u[i] ;

                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = beta * u[i+Nx+1]  -alpha * u[i-1]   + gamma * u[i] ;

                    }
                    else{
                        produit[i] = beta * u[i+Nx+1] -alpha * u[i-1]  + gamma * u[i] ;

                    }
                }

                if (i >= (Nx+1)*Ny )
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = beta * u[i+Nx+1-num_elements] -alpha * u[i+Nx] +  gamma * u[i] ;
                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = beta * u[i+Nx+1-num_elements] -alpha * u[i-1]  + gamma * u[i] ;
                    }
                    else{
                        produit[i] = beta * u[i+Nx+1-num_elements] -alpha * u[i-1] +  gamma * u[i] ;
                    }
                }


                if (i >= Nx+1 and i < (Nx+1)*Ny ){
                    
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = beta * u[i+Nx+1] -alpha * u[i+Nx] +  gamma * u[i] ;
                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = beta * u[i+Nx+1] -alpha * u[i-1] +  gamma * u[i] ;
                    }
                    else{
                        produit[i] =  beta * u[i+Nx+1] -alpha * u[i-1] +  gamma * u[i] ;
                    }


                }
                
                
            
        }
        } 
    }



    double x , y ,vx, vy ;

     if (space_scheme == 2 && time_scheme == 1) {
        

            for (int i = 0; i < num_elements; i++) {
                x = -1 + (i % (Nx +1) *dx) ;
                y = -1 + (i / (Nx +1 )* dy);
                vx = vxx(x,y,1.0,cas);
                vy = vyy(x,y,1.0,cas);

                double alpha = dt * (vx / (2*dx));
                double beta = dt * (vy / (2*dy));

                if ( i < Nx+1)
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = u[i] - alpha*u[i+Nx] - beta * u[i+num_elements-Nx-1] + alpha * u[i+1] + beta * u[i+Nx+1];
                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = u[i] + alpha * u[i-Nx] + beta * u[i+Nx+1]- beta * u[i+num_elements-Nx-1] - alpha*u[i-1];

                    }
                    else{
                        produit[i] = u[i] - alpha * u[i-1] - beta * u[i+num_elements-Nx-1] + alpha * u[i+1] + beta * u[i+Nx+1];

                    }
                }

                if (i >= (Nx+1)*Ny )
                {
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = u[i] - alpha * u[i+Nx] - beta * u[i-Nx-1] + alpha * u[i+1] + beta * u[i+Nx+1-num_elements];
                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = u[i] + alpha * u[i-Nx] + beta * u[i+Nx+1-num_elements]- beta * u[i-Nx-1] - alpha*u[i-1];
                    }
                    else{
                        produit[i] = u[i] -alpha * u[i-1] - beta * u[i-Nx-1] + alpha * u[i+1] + beta * u[i+Nx+1-num_elements];
                    }
                }


                if (i >= Nx+1 and i < (Nx+1)*Ny ){
                    
                    if ((i+1+Nx) % (Nx+1) == 0) {
                        
                        produit[i] = u[i] - alpha * u[i+Nx] - beta * u[i-Nx-1] + alpha * u[i+1] + beta * u[i+Nx+1];                        
                    }
                    else if ((i+1) % (Nx+1) == 0) {
                        
                        produit[i] = u[i] + alpha * u[i-Nx] + beta * u[i+Nx+1]- beta * u[i-Nx-1] - alpha*u[i-1];
                    }
                    else{
                        produit[i] = u[i] - alpha * u[i-1] - beta * u[i-Nx-1] + alpha * u[i+1] + beta * u[i+Nx+1];
                    }


                }
                
                
            
        }
    }

    return produit;
    

    
}
}

