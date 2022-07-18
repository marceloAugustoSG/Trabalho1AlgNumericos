#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <string.h>


#ifdef FUNCOES_H_INCLUDED
#define FUNCOES_H_INCLUDED

double maiorValorVetor(double* vetor, int dimensao);


void imprimirMatriz(double** matriz, int dimensao);
void imprimirVetoresB(double** matriz, int qtdSistemas, int dimensao);


void imprimirVetor(double* vetor, int dimensao);


double** multiplicacaoMatrizes(double** mat1, double** mat2, int dimensao);

void eliminacaoGauss(double** matriz, int dimensao, double* vetorB);

void fatoracaoLU(double** matriz, int dimensao, double* vetorB);

void jacobiIterativo(double** matriz, double* vetorB, int dimensao, double precisao);


void Gauss_Seidel(double** matriz, double* vetorB, int dimensao, double precisao);

#endif //FUNCOES_H_INCLUDED
