
#include <iostream>
#include <fstream>
#include <string>
#include "read.h"
#include <cmath> 



#include "fonctions.h"

//le cas test 4 est pour la solution manufacturée

double vxx(double x, double y, double t,int cas) {
    if (cas==1){
        return 1;
    }
    if (cas==2){
        return 0.5;
    }
    if (cas==3){
        return -y;
    }
    if (cas==4){
        return -1/(2*M_PI);
    }
    
}

double vyy(double x, double y, double t,int cas) {
    if (cas==1){
        return 0;
    }
    if (cas==2){
        return 0.5;
    }
    if (cas==3){
        return x;
    }
    if (cas==4){
        return 0;
    }
    
}

double u1(double x, double y, double t,double xmin, double ymin,int cas) {
    if (cas==1){
        double exponent_x = exp(-(x ) * (x ) / 0.0075);
        double exponent_y = exp(-(y) * (y ) / 0.0075);
        return exponent_x*exponent_y;
    }
    if (cas==2){
        if (sqrt(x*x +y*y)<=0.4){
            return 1;
        }
        else {
            return 0;
        }
        
    }
   if (cas==3){
        if (sqrt((x-0.4)*(x-0.4) +(y-0.4)*(y-0.4))>0.4){
            return 0;
        }
        else {
            return 1;
        }
        
    }
    if (cas==4){
        double exponent_x = sin(2*M_PI*x)*sin(2*M_PI*y);
        return exponent_x;
    }
}


