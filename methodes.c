#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "methodes.h"
#include "fonctions.h"

#define A(i,j) A[(i)*taille+(j)]
#define L(i,j) L[(i)*taille+(j)]
#define U(i,j) U[(i)*taille+(j)]
#define I(i,j) I[(i)*taille+(j)]
#define P(j) P[i*100 + j]

void LU(double* A, double* B, double* X, int taille){
	int i,j,k;
	double val;
	int *O = (int*)malloc((taille)*sizeof(int));
	double *L = (double*)malloc((taille*taille)*sizeof(double));
	double *U = (double*)malloc((taille*taille)*sizeof(double));
	double *Y = (double*)malloc((taille*taille)*sizeof(double));

	for(i=0; i<taille; i++){
		O[i]=i;
	}

	//Premiere colone de L
	for(i=0; i<taille; i++){
		L(i,0) = A(i,0);
	}

	meilleur_pivot(A,L,O, taille, 0);

	//Premiere ligne de U
	for(i=1; i<taille; i++){
		U(0,i) = A(0,i)/A(0,0);
	}

	//Complete L et U (sauf la valeur en bas a droite de L)
	for(i=1; i<taille-1; i++){
		val = 0;
		for(k=0; k<i; k++){
			val += (L(i,k)*U(k,i));
		}
		L(i,i) = A(i,i) - val;

		for(j=i+1; j<taille; j++){
			val = 0;
			for(k=0; k<i; k++){
				val +=L(j,k)*U(k,i);
			}
			L(j,i) = A(j,i) - val;
		}

		meilleur_pivot(A,L,O, taille, i);

		for(j=i+1; j<taille; j++){
			val = 0;
			for(k=0; k<i; k++){
				val +=L(i,k)*U(k,j);
			}
			U(i,j) = (A(i,j)-val)/L(i,i);
		}
	}

	//Derniere valeur de L
	val = 0;
	for(k=0; k<taille-1; k++){
		val +=L(taille-1,k)*U(k,taille-1);
	}
	L(taille-1, taille-1) = A(taille-1, taille-1) - val;

	//Diagonale de U
	for(i=0; i<taille; i++){
		U(i,i) = 1;
	}

	//Descente triangulaire
	Y[0] = B[O[0]]/L(0,0);
	for(i=1; i<taille; i++){
		val = 0;
		for(k=0; k<i; k++){
			val +=L(i,k)*Y[k];
		}
		if(L(i,i == 0)){
			Y[i] = 0;
		}
		else {
			Y[i] = (B[O[i]]-val)/L(i,i);
		}
	}

	//Remontee triangulaire
	X[taille-1] = Y[taille-1];
	for(i=taille-2; i>-1; i--){
		val = 0;
		for(k=i+1; k<taille+1; k++){
			val +=U(i,k)*X[k];
		}
		X[i] = Y[i] - val;
	}

	free(L); free(U); free(Y); free(O);
}

void LU_tridiag(double* A, double* B, double* X, int taille){
	int i;
	double *L = (double*)malloc((taille*taille)*sizeof(double));
	double *U = (double*)malloc((taille*taille)*sizeof(double));
	double *Y = (double*)malloc((taille*taille)*sizeof(double));

	for(i=0; i<2; i++){
		L(i,0)=A(i,0);
	}
	U(0,1) = A(0,1)/A(0,0);

	for(i=1; i<taille; i++){
		L(i,i) = A(i,i) - L(i,i-1)*U(i-1,i);
		L(i+1,i) = A(i+1,i);
		U(i,i+1) = A(i,i+1)/L(i,i);
	}

	L(taille-1,taille-1) = A(taille-1, taille-1) - L(taille-1, taille-2)*U(taille-2,taille-1);

	Y[0] = B[0]/L(0,0);

	for(i=1; i<taille-1; i++){
		Y[i]=(B[i]-L(i,i-1)*Y[i-1])/L(i,i);
	}

	X[taille-1] = Y[taille-1];

	for(i=taille-2; i>-1; i--){
		X[i]=Y[i]-U(i,i+1)*X[i+1];
	}

	free(L); free(U); free(Y);
}

double* spline_cubique_naturel (double* x, double* f, int taille){
	int i, j;
	double k;
	double* h = (double*)malloc((taille-1)*sizeof(double));
	double* A = (double*)malloc((taille*taille)*sizeof(double));
	double* B = (double*)malloc(taille*sizeof(double));
	double* P = (double*)malloc((taille-1)*100*sizeof(double));
	double* f_prim = (double*)malloc(taille*sizeof(double));
	double* f_sec = (double*)malloc(taille*sizeof(double));
	double* f_ter = (double*)malloc(taille*sizeof(double));

	for(i=0; i<taille-1; i++){
		h[i] = x[i+1]-x[i];
	}

	for(i=0; i<taille-1; i++){
		for(j=0; j<taille; j++){
			A(i,j) = 0;
		}
	}

	A(0,0) = 1;
	A(taille-1, taille-1) = 1;

	for(i=1; i<taille-1; i++){
		A(i,i-1) = h[i-1]/(h[i-1]+h[i]);
		A(i,i) = 2;
		A(i, i+1) = h[i]/(h[i-1]+h[i]);
	}

	for(i=0; i<taille-1; i++){
		B[i+1] = 6*(diff_div(f,x,i+1,i+2)-diff_div(f,x,i,i+1))/(x[i+2]-x[i]);
	}

	B[0] = 0;
	B[taille-1] = 0;

	LU_tridiag(A, B, f_sec, taille);

	for(i=0; i<taille-1; i++){
		f_prim[i] = diff_div(f,x,i,i+1) - (f_sec[i]/3)*h[i] - (f_sec[i+1]/6)*h[i];
		f_ter[i] = (f_sec[i+1]-f_sec[i])/h[i];
	}

	for(i=0; i<taille-1; i++){
		for(j=0; j<100; j++){
			k = x[i] + j* ((x[i+1] - x[i])/100);
			P(j) = f[i] + f_prim[i]*(k-x[i]) + (f_sec[i]/fact(2))*((k-x[i])*(k-x[i])) + (f_ter[i]/fact(3))*((k-x[i])*(k-x[i])*(k-x[i]));
		}
	}

	free(A); free(B); free(h); free(f_prim); free(f_sec); free(f_ter);

	return P;
}

double* spline_parametre (double* x, double* f, int taille){
	int i, j;
	double k;
	double* h = (double*)malloc((taille-1)*sizeof(double));
	double* A = (double*)malloc((taille*taille)*sizeof(double));
	double* B = (double*)malloc(taille*sizeof(double));
	double* P = (double*)malloc((taille-1)*100*sizeof(double));
	double* f_prim = (double*)malloc(taille*sizeof(double));
	double* f_sec = (double*)malloc(taille*sizeof(double));
	double* f_ter = (double*)malloc(taille*sizeof(double));

	for(i=0; i<taille-1; i++){
		h[i] = x[i+1]-x[i];
	}

	for(i=0; i<taille-1; i++){
		for(j=0; j<taille; j++){
			A(i,j) = 0;
		}
	}

	A(0,0) = h[0]/3;
	A(0,1) = h[0]/6;
	A(0,taille-2) = h[taille-2]/6;
	A(0, taille-1) = h[taille-2]/3;
	A(taille-1,0) = h[0]/3;
	A(taille-1,1) = h[0]/6;
	A(taille-1,taille-2) = h[taille-2]/6;
	A(taille-1, taille-1) = h[taille-2]/3;

	for(i=1; i<taille-1; i++){
		A(i,i-1) = h[i-1]/(h[i-1]+h[i]);
		A(i,i) = 2;
		A(i, i+1) = h[i]/(h[i-1]+h[i]);
	}

	for(i=0; i<taille-1; i++){
		B[i+1] = 6*(diff_div(f,x,i+1,i+2)-diff_div(f,x,i,i+1))/(x[i+2]-x[i]);
	}

	B[0] = diff_div(f,x,0,1)-diff_div(f,x,taille-2,taille-1);
	B[taille-1] = diff_div(f,x,0,1)-diff_div(f,x,taille-2,taille-1);

	LU(A, B, f_sec, taille);

	for(i=0; i<taille; i++){
		printf("%lf\n", f_sec[i]);
	}

	for(i=0; i<taille-1; i++){
		f_prim[i] = diff_div(f,x,i,i+1) - (f_sec[i]/3)*h[i] - (f_sec[i+1]/6)*h[i];
		f_ter[i] = (f_sec[i+1]-f_sec[i])/h[i];
	}

	for(i=0; i<taille-1; i++){
		for(j=0; j<100; j++){
			k = x[i] + j* ((x[i+1] - x[i])/100);
			P(j) = f[i] + f_prim[i]*(k-x[i]) + (f_sec[i]/fact(2))*((k-x[i])*(k-x[i])) + (f_ter[i]/fact(3))*((k-x[i])*(k-x[i])*(k-x[i]));
		}
	}

	free(A); free(B); free(h); free(f_prim); free(f_sec); free(f_ter);

	return P;
}
