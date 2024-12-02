#ifndef FONCTIONS_H
#define FONCTIONS_H


using namespace std;

double vxx(double x, double y, double t, int cas);
double vyy(double x, double y, double t,int cas);
double u1(double x, double y, double t,double xmin,double ymin,int cas);
void charge_a(int me, int n, int Np, int &iBeg, int &iEnd);


#endif 
