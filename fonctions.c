#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fonctions.h"

#define A(i,j) A[(i)*taille+(j)]
#define L(i,j) L[(i)*taille+(j)]

void meilleur_pivot(double* A, double*L, int* O, int taille, int colone){

	//Recherche d'un meilleur pivot dans la colone
	 int i, idx=0;
	 double val;

	 for(i=1; i<taille; i++){
		 if(L(i,colone)>L(idx,colone)){
			 idx = i;
		 }
	 }

	 //s'il y a un meilleur pivot on echange les lignes
	 if(idx != 0){
		 val = O[idx];
		 O[idx] = O[colone];
		 O[colone] = val;

		 for(i=0; i<taille; i++){
			 val = L(idx,i);
			 L(idx,i) = L(colone,i);
			 L(colone,i) = val;

			 val = A(idx,i);
			 A(idx,i) = A(colone,i);
			 A(colone,i) = val;
		 }
	 }
}

double diff_div(double* f, double* x, int idx1, int idx2){
	return ((f[idx2]-f[idx1])/(x[idx2]-x[idx1]));
}

unsigned int fact(unsigned int x){
	if(x == 0){
		return 0;
	}

	unsigned int result = 1;
	for(unsigned int i = 1; i<x+1; i++){
		result = result * i;
	}

	return result;
}
