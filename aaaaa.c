#include <stdio.h>
#include <stdlib.h>
#include <string.h>




void alocarELerMatriz(int** matriz, int linhas, int colunas) {

  matriz = (int**) malloc(sizeof(int*) * linhas);

  for(int i = 0; i < linhas; i++){
    matriz[i] = (int*) malloc(sizeof(int) * (colunas));

    for(int j = 0; j < colunas; j++){
      scanf("%d", &matriz[i][j]);
    }
  }

  return;

}


void freeMatriz(int** matriz, int linhas) {

  for(int i = 0; i < linhas; i++){
    free(matriz[i]);
  }
  free(matriz);

}


int main() {

  char operacao[13];
  int linhas;
  int** matriz;

  scanf("%s", &operacao);
  scanf("%d", &linhas);

  if(!strcmp(operacao, "resolve")){
    
    alocarELerMatriz(matriz, linhas, linhas+1);

    // Resolucao do sistema

    freeMatriz(matriz, linhas);

    return 0;
  }

  if(!strcmp(operacao, "inverte")){

    alocarELerMatriz(matriz, linhas, linhas);

    // Inversao de matriz

    freeMatriz(matriz, linhas);

    return 0;
  }

  if(!strcmp(operacao, "determinante")){

    alocarELerMatriz(matriz, linhas, linhas);

    // Calculo de determinante

    freeMatriz(matriz, linhas);

    return 0;
  }

  return 0;

}
