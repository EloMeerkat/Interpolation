#define _CRT_SECURE_NO_WARNINGS

#define A(i,j) A[(i)*taille+(j)]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methodes.h"
#include "fonctions.h"

int main(int argc, char **argv){
	printf("Debut de l'algorithme\n");
	printf("*********************\n");
	printf("Lecture des donnees\n");

	FILE *in = fopen("entree.txt", "rt");
	FILE *out = fopen("reponse.csv", "wt");

	int nb_ligne;
	double valeur;
	char chaine[5000];
	fgets(chaine, 5000, in);
	char* p = chaine;
	sscanf(p, "%d", &nb_ligne);

//	Definition des matrices
//	double *A = (double*)malloc((taille*taille)*sizeof(double));
//	double *B = (double*)malloc((taille)*sizeof(double));
//	double *X = (double*)malloc((taille)*sizeof(double));
//	int *O = (int*)malloc((taille)*sizeof(int));

	double *t = (double*)malloc((nb_ligne)*sizeof(double));
	double *x = (double*)malloc((nb_ligne)*sizeof(double));
	double *y = (double*)malloc((nb_ligne)*sizeof(double));
	double *X = (double*)malloc(((nb_ligne-1)*100)*sizeof(double));
	double *Y = (double*)malloc(((nb_ligne-1)*100)*sizeof(double));

	for(int i=0; i<nb_ligne; i++){
		fgets(chaine, 5000, in);
		p = chaine;
		sscanf(p, "%lf", &valeur);
		t[i]=valeur;
	}
	for(int i=0; i<nb_ligne; i++){
		fgets(chaine, 5000, in);
		p = chaine;
		sscanf(p, "%lf", &valeur);
		x[i]=valeur;
	}
	for(int i=0; i<nb_ligne; i++){
		fgets(chaine, 5000, in);
		p = chaine;
		sscanf(p, "%lf", &valeur);
		y[i]=valeur;
	}

	//X = spline_cubique_naturel (t, x, nb_ligne);
	//Y = spline_cubique_naturel (t, y, nb_ligne);
	X = spline_parametre (t, x, nb_ligne);
	Y = spline_parametre (t, y, nb_ligne);

	fprintf(out, "x;y \n");

	for(int i=0; i<(nb_ligne-1)*100; i++){
		fprintf(out, "%lf;%lf\n", X[i], Y[i]);
	}

	fclose(out);
	free(t); free(x); free(y);
	printf("Fin\nAppuyer sur entrer pour quitter.\n");
	getchar();

	return 0;
}
