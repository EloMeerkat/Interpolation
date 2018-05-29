#ifndef _METHODES_
#define _METHODES_

void LU(double* A, double* B, double* X, int taille);
void LU_tridiag(double* A, double* B, double* X, int taille);
double* spline_cubique_naturel (double* x, double* f, int taille);
double* spline_parametre (double* x, double* f, int taille);

#endif
