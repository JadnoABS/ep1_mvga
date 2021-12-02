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
    mat->m[i1][j] = mat->m[i1][j] + (mat->m[i2][j] * k);
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
	
		for(j = 0; j < mat->col; j++) if(fabs(mat->m[i][j]) > 0) break;

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
  for(int i = 0; i < mat->lin; i++){

    int linha_pivo; int coluna_pivo;
    encontra_linha_pivo(mat, i, &linha_pivo, &coluna_pivo);

    if(linha_pivo != i){
      troca_linha(mat, i, linha_pivo);
      if(agregada) troca_linha(agregada, i, linha_pivo);
    }

    for(int j = i+1; j < mat->lin; j++){
      if(abs(mat->m[j][coluna_pivo]) > SMALL){
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

    multiplica_linha(mat, i, 1/(mat->m[i][coluna_pivo]));
    if(agregada) multiplica_linha(agregada, i, 1/(mat->m[i][coluna_pivo]));
  }

	return determinante;
}

// implementa a eliminacao de Gauss-Jordan que coloca a matriz quadrada "mat" na forma 
// escalonada reduzida. As operacoes realizadas para colocar a matriz "mat" na forma 
// escalonada reduzida tambem devem ser aplicadas na matriz "agregada" caso esta seja 
// nao nula. Nao se pode assumir que "mat" jah esteja na forma escalonada (mas voce 
// pode usar a funcao acima para isso).

void forma_escalonada_reduzida(Matriz * mat, Matriz * agregada){

  if(agregada){ imprime_matrizes(mat, agregada); printf("\n"); }
  else { imprime_matriz(mat); printf("\n"); }

  for(int i = 0; i < mat->lin; i++){

    int linha_pivo; int coluna_pivo;
    encontra_linha_pivo(mat, i, &linha_pivo, &coluna_pivo);

    if(linha_pivo != i){
      troca_linha(mat, i, linha_pivo);
      if(agregada) troca_linha(agregada, i, linha_pivo);
    }
    multiplica_linha(mat, i, 1/(mat->m[i][coluna_pivo]));
    if(agregada) multiplica_linha(agregada, i, 1/(mat->m[i][coluna_pivo]));

    if(agregada){ imprime_matrizes(mat, agregada); printf("\n"); }
    else { imprime_matriz(mat); printf("\n"); }

    for(int j = 0; j < mat->lin; j++){
      if(j == i) continue;
      double valorMult = -(mat->m[j][coluna_pivo]/mat->m[i][coluna_pivo]);
      combina_linhas(mat, j, i, valorMult);
      if(agregada) combina_linhas(agregada, j, i, valorMult);
    }

    if(agregada){ imprime_matrizes(mat, agregada); printf("\n"); }
    else { imprime_matriz(mat); printf("\n"); }
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



  // imprime_matriz(mat);
  // printf("\n");

  // forma_escalonada_reduzida(mat, NULL);

  // imprime_matriz(mat);
  // printf("\n");
	
	if(strcmp("resolve", operacao) == 0){
    Matriz* mat = cria_matriz(n,n+1);
    lerMatriz(mat);
    forma_escalonada_reduzida(mat, NULL);

    for(int i = 0; i < mat->lin; i++)
      printf("%.2lf\n", mat->m[i][mat->col-1]);

	}
	else if(strcmp("inverte", operacao) == 0){
    Matriz* mat = cria_matriz(n,n);
    lerMatriz(mat);

    Matriz* inv = cria_identidade(n);

    forma_escalonada_reduzida(mat, inv);

    imprime_matriz(inv);
	}
	else if(strcmp("determinante", operacao) == 0){
    Matriz* mat = cria_matriz(n,n);
    lerMatriz(mat);

    double determinante = forma_escalonada(mat, NULL);
    printf("%.2lf\n", determinante);

	} else {
		printf("Operação desconhecida!\n");
		return 1;
  }


	return 0;
}

