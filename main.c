#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include "funcoes.h"
#include <string.h>
#include <stdint.h>

int main(int argc, char* argv[]) {

    int dimensao, qtdSistemas;
    double precisao;
    char isSimetrico;
    double** matrizSistema;
    double** vetoresB;
    double* vetorB;
    char concatena[30];
    char* dat = ".dat";
    char* nomeArquivo = argv[1];
    clock_t h;

    if (argc > 1) {

        strcat(strcpy(concatena, nomeArquivo), dat);


    }

    FILE* arq = fopen(concatena, "r");

    if (arq == NULL) {

        printf("\nArquivo nao encontrado ou invalido\n");
    }
    else {

        printf("\nArquivo encontrado\n");

        fscanf(arq, "%d %d %lf %c", &qtdSistemas, &dimensao, &precisao, &isSimetrico);
        printf("\n%d %d %lf %c\n", qtdSistemas, dimensao, precisao, isSimetrico);

        matrizSistema = (double**)malloc(dimensao * sizeof(double));
        for (int i = 0;i < dimensao;i++) {
            matrizSistema[i] = malloc(dimensao * sizeof(double));
        }

        for (int i = 0;i < dimensao;i++) {
            for (int j = 0; j < dimensao;j++) {
                fscanf(arq, "%lf", &matrizSistema[i][j]);

            }
        }

        imprimirMatriz(matrizSistema, dimensao);

        vetoresB = (double**)malloc(qtdSistemas * sizeof(double));
        for (int i = 0;i < dimensao;i++) {
            vetoresB[i] = malloc(dimensao * sizeof(double));
        }

        for (int i = 0;i < qtdSistemas;i++) {
            for (int j = 0; j < dimensao;j++) {
                fscanf(arq, "%lf", &vetoresB[i][j]);

            }
        }

        printf("\nVetoresB:\n");
        imprimirVetoresB(vetoresB, qtdSistemas, dimensao);

        vetorB = malloc(dimensao * sizeof(double));
        for (int i = 0;i < qtdSistemas;i++) {
            printf("\n=====================================================================\n");
            printf("\nSistema %d\n", i + 1);
            for (int j = 0; j < dimensao;j++) {
                vetorB[j] = vetoresB[i][j];

            }
            printf("\nVetorB\n");
            imprimirVetor(vetorB, dimensao);
            printf("\n");


            //por alguma razao o jacobi tem que ser o primeiro
            printf("\n-----------------------------------------\n");
            h = clock();
            jacobiIterativo(matrizSistema, vetorB, dimensao, precisao);
            h = clock() - h;
            printf("\nTempo jacobi Iterativo Gauss: %f\n", (float)h / CLOCKS_PER_SEC);
            printf("\n-----------------------------------------\n");
            h = clock();
            eliminacaoGauss(matrizSistema, dimensao, vetorB);
            h = clock() - h;
            printf("\nTempo Eliminacao Gauss: %f\n", (float)h / CLOCKS_PER_SEC);
            printf("\n-----------------------------------------\n");
            h = clock();
            fatoracaoLU(matrizSistema, dimensao, vetorB);
            h = clock() - h;
            printf("\nTempo Fatoracao LU: %f\n", (float)h / CLOCKS_PER_SEC);
            printf("\n-----------------------------------------\n");
            h = clock();
            Gauss_Seidel(matrizSistema, vetorB, dimensao, precisao);
            h = clock() - h;
            printf("\nTempo Gauss-Seidel: %f\n", (double)h / CLOCKS_PER_SEC);
        }



    }
        fclose(arq);


    return 0;
}