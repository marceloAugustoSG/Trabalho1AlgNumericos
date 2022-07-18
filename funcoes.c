#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include "funcoes.h"
#include <string.h>

double maiorValorVetor(double* vetor, int dimensao) {
    double maior = vetor[0];
    for (int i = 1;i < dimensao;i++) {
        if (vetor[i] > maior) {
            maior = vetor[i];
        }
    }

    return maior;

}

void imprimirMatriz(double** matriz, int dimensao) {

    for (int l = 0;l < dimensao;l++) {
        for (int c = 0;c < dimensao;c++) {
            printf("%lf ", matriz[l][c]);

        }
        printf("\n");
    }
}

void imprimirVetoresB(double** matriz, int qtdSistemas, int dimensao) {

    for (int l = 0;l < qtdSistemas;l++) {
        for (int c = 0;c < dimensao;c++) {
            printf("%lf ", matriz[l][c]);

        }
        printf("\n");
    }
}


void imprimirVetor(double* vetor, int dimensao) {
    printf("{");
    for (int i = 0;i < dimensao;i++) {
        printf(" %lf, ", vetor[i]);

    }
    printf("}");
    printf("\n");

}
//funcao que retorna uma nova matriz da multiplicacao de duas matrizes
double** multiplicacaoMatrizes(double** mat1, double** mat2, int dimensao) {

    double soma = 0.0;
    double** matResultado;
    //alocando a matriz resultado
    matResultado = (double**)malloc(dimensao * sizeof(double*));
    for (int i = 0;i < dimensao;i++) {
        matResultado[i] = malloc(dimensao * sizeof(double));
    }

    for (int l = 0;l < dimensao;l++) {
        for (int c = 0;c < dimensao;c++) {
            soma = 0;
            for (int k = 0;k < dimensao;k++) {
                soma = soma + mat1[l][k] * mat2[k][c];
            }
            matResultado[l][c] = soma;
        }
    }

    return matResultado;

}

void eliminacaoGauss(double** matriz, int dimensao, double* vetorB) {

    int i, j, k;
    double* X;
    double M, soma = 0.0;


    //zerando abaixo da diagonal principal
    for (i = 0;i < dimensao - 1;i++) {
        for (j = i + 1;j < dimensao;j++) {
            M = matriz[j][i] / matriz[i][i];
            for (k = 0;k < dimensao;k++) {
                matriz[j][k] = matriz[j][k] - M * (matriz[i][k]);
            }
            vetorB[j] = vetorB[j] - M * vetorB[i];

        }

    }



    //vetor respostas
    X = (double*)malloc(dimensao * sizeof(double));
    for (i = 0;i < dimensao;i++) {
        X[i] = 0.0;
    }

    for (j = dimensao - 1;j >= 0;j--) {
        soma = 0.0;
        for (k = j;k < dimensao - 1;k++) {
            soma = soma + matriz[j][k + 1] * X[k + 1];
        }
        X[j] = (vetorB[j] - soma) / (matriz[j][j]);
    }

    printf("Vetor X eliminacao de Gauss:\n");
    imprimirVetor(X, dimensao);
    printf("\n");

}

void fatoracaoLU(double** matriz, int dimensao, double* vetorB) {

    double** matrizL;
    double** matrizU;
    double** somaLU;
    double* vetorY;
    double* vetorX;

    //alocar as matrizes L e U,VetorX  e VetorY

    matrizU = (double**)malloc(dimensao * sizeof(double*));
    for (int i = 0;i < dimensao;i++) {
        matrizU[i] = malloc(dimensao * sizeof(double));

    }

    matrizL = (double**)malloc(dimensao * sizeof(double*));
    for (int i = 0;i < dimensao;i++) {
        matrizL[i] = malloc(dimensao * sizeof(double));

    }

    somaLU = (double**)malloc(dimensao * sizeof(double*));
    for (int i = 0;i < dimensao;i++) {
        somaLU[i] = malloc(dimensao * sizeof(double));

    }
    //alocando espaco pros vetores que irei precisar
    vetorY = malloc((sizeof(double)) * dimensao);
    vetorX = malloc((sizeof(double)) * dimensao);


    // copiar a matriz a em u
    for (int i = 0;i < dimensao;i++) {
        for (int j = 0;j < dimensao;j++) {
            matrizU[i][j] = matriz[i][j];
        }
    }

    //matriz identidade
    for (int i = 0;i < dimensao;i++) {
        for (int j = 0;j < dimensao;j++) {
            if (i == j) {
                matrizL[i][j] = 1;
            }
            else {
                matrizL[i][j] = 0;
            }
        }
    }

    //transformando a matriz U em matriz superior
    for (int j = 0;j < dimensao;j++) {
        for (int i = j + 1;i < dimensao;i++) {
            if (matrizU[i][j] != 0) {
                matrizL[i][j] = matrizU[i][j] / matrizU[j][j];
                for (int c = j;c < dimensao;c++) {
                    matrizU[i][c] = matrizU[i][c] + matrizU[j][c] * (-1 * (matrizL[i][j]));
                }
            }
        }
    }

    //obtendo a matriz LU para verificar se é a matriz do sistema
    somaLU = multiplicacaoMatrizes(matrizL, matrizU, dimensao);

    /*   printf("\nmatriz L*U:\n");
      imprimirMatriz(somaLU, dimensao); */


      //calcular o vetor Y a a partir de L x Y = B

    //printf("\n-----------------------------------------------------------------------\n");

    double somaLinha;

    //calculando o vetorY 
    for (int i = 0;i < dimensao;i++) {
        somaLinha = 0.0;
        for (int j = i - 1;j >= 0;j--) {
            somaLinha += matrizL[i][j] * vetorY[j];
        }
        vetorY[i] = (vetorB[i] - somaLinha) / matrizL[i][i];
    }

    //imprimindo para teste
    //imprimirVetor(vetorY, dimensao);


    //calculando o vetor X

    for (int i = dimensao - 1;i >= 0;i--) {
        somaLinha = 0.0;
        for (int j = i + 1;j <= dimensao - 1;j++) {
            somaLinha += matrizU[i][j] * vetorX[j];
        }
        vetorX[i] = (vetorY[i] - somaLinha) / matrizU[i][i];
    }

    printf("Resposta Fatoracao LU:\n");
    imprimirVetor(vetorX, dimensao);
    printf("\n");

}

void jacobiIterativo(double** matriz, double* vetorB, int dimensao, double precisao) {

    double vetorX0[dimensao];
    double* next;
    double* vetorResposta;
    double** matrizResposta;
    double vetorX0X1[dimensao];


    matrizResposta = (double**)malloc(dimensao * sizeof(double*));
    for (int i = 0;i < dimensao;i++) {
        matrizResposta[i] = malloc(dimensao * sizeof(double));
    }

    vetorResposta = malloc(dimensao * sizeof(double));
    next = malloc(dimensao * sizeof(double));


    //definindo o 1º vetor X0
    for (int i = 0;i < dimensao;i++) {
        for (int j = 0;j < dimensao;j++) {
            if (i == j) {
                vetorX0[i] = matriz[i][j];
            }

        }

    }

    for (int i = 0;i < dimensao;i++) {
        for (int j = 0;j < dimensao;j++) {
            if (i == j) {
                vetorX0[i] = vetorB[i] / vetorX0[i];
            }

        }

    }

    double maiorsubtracaoX1X0, maiorX1, divX1X0 = precisao + 1;

    while (divX1X0 > precisao)
    {

        // for (int k = 0;k < iteracoes;k++) {
        for (int i = 0;i < dimensao;i++) {
            double bi = vetorB[i];
            for (int j = 0;j < dimensao;j++) {
                if (j != i) {
                    bi -= matriz[i][j] * vetorX0[j];
                }
            }
            bi /= matriz[i][i];
            // printf("X%d ^ (%d) = %.2lf\t", i + 1, k + 1, bi);
            next[i] = bi;


            vetorX0X1[i] = next[i] - vetorX0[i];


            maiorX1 = maiorValorVetor(next, dimensao);
            maiorsubtracaoX1X0 = maiorValorVetor(vetorX0X1, dimensao);
            divX1X0 = fabs(maiorsubtracaoX1X0) / fabs(maiorX1);


        }

        if (divX1X0 < precisao) {

            printf("Vetor X Jacobi:\n");
            imprimirVetor(vetorResposta, dimensao);
            printf("\n");


            break;

        }
        else {
            //atualizar o X0
            for (int i = 0;i < dimensao;i++) {

                vetorX0[i] = next[i];
                vetorResposta[i] = vetorX0[i];
            }
        }
    }




}

void Gauss_Seidel(double** matriz, double* vetorB, int dimensao, double precisao) {

    double vetorX0[dimensao];

    double* vetorResposta;
    double** matrizResposta;
    double vetorX0X1[dimensao];

    matrizResposta = (double**)malloc(dimensao * sizeof(double*));
    for (int i = 0;i < dimensao;i++) {
        matrizResposta[i] = malloc(dimensao * sizeof(double));
    }

    vetorResposta = malloc(dimensao * sizeof(double));


    //definindo o 1º vetor X0
    for (int i = 0;i < dimensao;i++) {
        for (int j = 0;j < dimensao;j++) {
            if (i == j) {
                vetorX0[i] = matriz[i][j];
            }

        }

    }
    for (int i = 0;i < dimensao;i++) {
        for (int j = 0;j < dimensao;j++) {
            if (i == j) {
                vetorX0[i] = vetorB[i] / vetorX0[i];
            }

        }

    }
    double maiorsubtracaoX1X0, maiorX1, divX1X0 = precisao + 1;

    while (divX1X0 > precisao)
    {

        for (int i = 0;i < dimensao;i++) {
            double bi = vetorB[i];
            for (int j = 0;j < dimensao;j++) {
                if (j != i) {
                    bi -= matriz[i][j] * vetorX0[j];
                }
            }
            bi /= matriz[i][i];

            // daqui pra baixo ta quase igual ao de jacobi
            vetorX0[i] = bi; // ja atualizar o X0
            vetorX0X1[i] = vetorX0[i] - vetorX0[i];
            maiorX1 = maiorValorVetor(vetorX0, dimensao);
            maiorsubtracaoX1X0 = maiorValorVetor(vetorX0X1, dimensao);
            divX1X0 = fabs(maiorsubtracaoX1X0) / fabs(maiorX1);
            vetorResposta[i] = vetorX0[i];
        }
        // nao obtive uma precisao igual ao do jacobi em meus testes
        //dependendo do sistema ele tem uma precisao melhor 

    }
    printf("Vetor X Gauss Seidel:\n");
    imprimirVetor(vetorResposta, dimensao);
    printf("\n");

}
