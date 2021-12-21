/*
  ##############################
  EP1 - ACH2033-203-2021
  Professor: Flávio Luiz Coutinho
  Alunos:
    Ana Clara Diamantino de Vasconcelos - 12674398
    Jadno Augusto Barbosa da Silva - 12608618
  ##############################
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// constante para ser usada na comparacao de valores double.
// Se a diferenca absoluta entre dois valores double for menor
// do que o valor definido por esta constante, eles devem ser
// considerados iguais.

#define SMALL 0.000001

// struct que representa uma matriz de valores do tipo double.

typedef struct {

	double ** m;
	int lin, col;

} Matriz;


// cria uma matriz de n linhas por m colunas com todas as entradas iguais a zero.

Matriz * cria_matriz(int n, int m){

	int i;

	Matriz * mat = (Matriz *) malloc (sizeof(Matriz));

	mat->lin = n;
	mat->col = m;
	mat->m = (double **) malloc(mat->lin * sizeof(double *));

	for(i = 0; i < mat->lin; i++){

		mat->m[i] = (double *) malloc(mat->col * sizeof(double));
		memset(mat->m[i], 0, mat->col * sizeof(double));
	}

	return mat;
}

// cria uma matriz identidade de tamanho n x n.

Matriz * cria_identidade(int n){

	int i;

	Matriz * mat = cria_matriz(n, n);

	for(i = 0; i < mat->lin; i++) mat->m[i][i] = 1;

	return mat;
}

// imprime a matriz passada como parametro.

void imprime_matriz(Matriz * mat){

	int i, j;

	if(!mat) return;

	for(i = 0; i < mat->lin; i++){

		for(j = 0; j < mat->col; j++){

			printf("%7.2f ", mat->m[i][j]);
		}

		printf("\n");
	}
}

// imprime a matriz expandida formada pela combinacao das matrizes "mat" e "agregada". 
// Isto eh, cada linha da matriz impressa possui as entradas da linha correspondente 
// em "mat", seguida das entradas da linha correspondente em "agregada".

void imprime_matrizes(Matriz * mat, Matriz * agregada){

	int i, j;

	if(!mat || !agregada) return;

	for(i = 0; i < mat->lin; i++){

		for(j = 0; j < mat->col; j++){

			printf("%7.2f ", mat->m[i][j]);
		}

		printf(" |");

		for(j = 0; j < agregada->col; j++){

			printf("%7.2f ", agregada->m[i][j]);
		}

		printf("\n");
	}
}

// funcao que troca as linhas i1 e i2 de lugar.

void troca_linha(Matriz * mat, int i1, int i2){

  double aux;

  for(int j = 0; j < mat->col; j++){
    aux = mat->m[i1][j];
    mat->m[i1][j] = mat->m[i2][j];
    mat->m[i2][j] = aux;
  }

  return;
}

// funcao que multiplica as entradas da linha i pelo escalar k

void multiplica_linha(Matriz * mat, int i, double k){

  for(int j = 0; j < mat->col; j++)
    mat->m[i][j] = mat->m[i][j] * k;

  return;
}

// funcao que faz a seguinte combinacao de duas linhas da matriz:
//
// 	(linha i1) = (linha i1) + (linha i2 * k)
//

void combina_linhas(Matriz * mat, int i1, int i2, double k){

  for(int j = 0; j < mat->col; j++){
    mat->m[i1][j] += (mat->m[i2][j] * k);
  }

  return;
}

// funcao que procura, a partir da linha ini, a linha com uma entrada nao nula que
// esteja o mais a esquerda possivel dentre todas as linhas. Os indices da linha e da 
// coluna referentes a entrada nao nula encontrada sao devolvidos através dos ponteiros 
// "linha_pivo" e "coluna_pivo". Esta funcao ja esta pronta para voces usarem na 
// implementacao da eliminacao gaussiana e eleminacao de Gauss-Jordan.

void encontra_linha_pivo(Matriz * mat, int ini, int * linha_pivo, int * coluna_pivo){

	int i, j;

	*linha_pivo = mat->lin;
	*coluna_pivo = mat->col;

	for(i = ini; i < mat->lin; i++){

		for(j = 0; j < mat->col; j++) if(fabs(mat->m[i][j]) > SMALL) break;

		if(j < *coluna_pivo) {

			*linha_pivo = i;
			*coluna_pivo = j;
		}
	}
}


// implementa a eliminacao gaussiana que coloca a matriz quadrada "mat" na forma escalonada.
// As operacoes realizadas para colocar a matriz "mat" na forma escalonada tambem devem ser
// aplicadas na matriz "agregada" caso esta seja nao nula. Esta funcao tambem deve calcular
// e devolver o determinante de "mat".

double forma_escalonada(Matriz *  mat, Matriz * agregada){

  // Operacoes elementares para triangularizar a matriz
  int linha_pivo; int coluna_pivo;
  for(int i = 0; i < mat->lin; i++){

    encontra_linha_pivo(mat, i, &linha_pivo, &coluna_pivo);

    if (linha_pivo == mat->lin) return 0.0;
    if(linha_pivo != i){
      troca_linha(mat, i, linha_pivo);
      if(agregada) troca_linha(agregada, i, linha_pivo);
    }

    for(int j = i+1; j < mat->lin; j++){
      if(fabs(mat->m[j][coluna_pivo]) > SMALL){
        
        double valorMult = -(mat->m[j][coluna_pivo]) / (mat->m[i][coluna_pivo]);
        combina_linhas(mat, j, i, valorMult);
        if(agregada) combina_linhas(agregada, j, i, valorMult);
      }
    }
  }

  // Calculo do determinante da matriz triangular
  // Logo apos, operacoes para construir uma diagonal principal com elementos = 1
  double determinante = 1;

  for(int i = 0; i < mat->lin; i++){
    determinante *= mat->m[i][i];
    
    int linha_pivo; int coluna_pivo;
    encontra_linha_pivo(mat, i, &linha_pivo, &coluna_pivo);

    double valorMult = 1 / mat->m[i][coluna_pivo];
    multiplica_linha(mat, i, valorMult);
    if(agregada) multiplica_linha(agregada, i, valorMult);
  }

  return determinante;
}

// implementa a eliminacao de Gauss-Jordan que coloca a matriz quadrada "mat" na forma
// escalonada reduzida. As operacoes realizadas para colocar a matriz "mat" na forma
// escalonada reduzida tambem devem ser aplicadas na matriz "agregada" caso esta seja
// nao nula. Nao se pode assumir que "mat" jah esteja na forma escalonada (mas voce
// pode usar a funcao acima para isso).

void forma_escalonada_reduzida(Matriz * mat, Matriz * agregada){

  forma_escalonada(mat, agregada);
  
  if(!agregada) {
    if(fabs(forma_escalonada(mat, agregada)) < SMALL) {
      if(fabs(mat->m[mat->lin - 1][mat->col - 1]) < SMALL) {
        printf("sistema possui diversas soluções\n"); 
        return;
      }
      else {
        printf("sistema sem solução\n"); 
        return;
      }
    }
  }

  else {
    for(int i = mat->lin-1; i > 0; i--){

      int linha_pivo; int coluna_pivo;
      encontra_linha_pivo(mat, i, &linha_pivo, &coluna_pivo);
      

      for(int j = i-1; j >= 0; j--){

        double valorMult = -(mat->m[j][coluna_pivo]) / (mat->m[i][coluna_pivo]);
        combina_linhas(mat, j, i, valorMult);
        if(agregada) combina_linhas(agregada, j, i, valorMult);
        
      }

    }
  }
}

void lerMatriz(Matriz* mat) {
  for(int i = 0; i < mat->lin; i++){
    for(int j = 0; j < mat->col; j++){
      scanf("%lf", &mat->m[i][j]);
    }
  }
}

// funcao principal.

int main(){

	int n;
	char operacao[32];

	scanf("%s", operacao);	// le, a partir da entrada padrao, a string que determina qual operacao deve ser realizada.
	scanf("%d", &n);	// le a dimensão da matriz a ser manipulada pela operacao escolhida.


	if(strcmp("resolve", operacao) == 0){
    Matriz* mat = cria_matriz(n,n+1);
    lerMatriz(mat);
    forma_escalonada_reduzida(mat, NULL);
    if(forma_escalonada(mat, NULL) != 0){
    for(int i = 0; i < mat->lin; i++)
      printf("%.2lf\n", mat->m[i][mat->col-1]);
    }

    if(mat) free(mat);
	}
	else if(strcmp("inverte", operacao) == 0){
    Matriz* mat = cria_matriz(n,n);
    lerMatriz(mat);

    Matriz* inv = cria_identidade(n);

    forma_escalonada_reduzida(mat, inv);
    if(forma_escalonada(mat, inv) == 0) {
      printf("matriz singular");
    }
    else {
      imprime_matriz(inv);
    }

    if(mat) free(mat);
    if(inv) free(inv);
	}
	else if(strcmp("determinante", operacao) == 0){
    Matriz* mat = cria_matriz(n,n);
    lerMatriz(mat);

    double determinante = forma_escalonada(mat, cria_identidade(n));
    printf("%.2lf\n", determinante);

    if(mat) free(mat);
	} else {
		printf("Operação desconhecida!\n");
		return 1;
  }


	return 0;
}
