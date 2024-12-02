#ifndef BiCGstab_H
#define BiCGstab_H

#include "Matrice.h"
using namespace std;
vector<double> BiCGstab(int Nx, int Ny, double dx, double dy, double vx, double vy, double dt, int space_scheme, int time_scheme, const vector<double>& b,  double eProd_scailon, int Nmax, int iBeg, int iEnd, int Np, int Me,int cas) ;

#endif