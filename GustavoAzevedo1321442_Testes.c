/**
Gustavo Azevedo - 1321442
Trabalho Final Análise Numérica (INF-1608)
"Movimento de um Pêndulo"

TESTES
*/
#include <math.h>
#include <stdio.h>
#include "GustavoAzevedo1321442_Pendulo.h"

/** Testes para os metodos numericos*/
int main(void){
	double hRk4 = 0.0001, hRk5 = 0.001;
	double theta = M_PI/20.0;	
	double tolerancia = pow(10, -5);
	double esperado;
	double resRk4, erroRk4, resRk5, erroRk5;
	int passosRk4, passosRk5;
	FILE* output;

	output = fopen("theta_tempo.csv","w");
	if(output==NULL) return 2;

	fprintf(output,"Theta0,PERIODO_RK4,NUM_PASSOS_RK4,PERIODO_RK5,NUM_PASSOS_RK5,PERIODO_REAL\n");
	
	while(theta<=M_PI/2.0){
		printf("Theta: %f\n", theta);
		fprintf(output,"%f,", theta);
		
		esperado = periodoEsperado(theta);
 
		//Fixo
		passosRk4 = 0;
		resRk4 = RungeKuttaPassoFixo(hRk4, theta, &passosRk4, 0);
		erroRk4 = fabs(resRk4 - esperado);
		printf("  RK4 - Esperado:%.4f - Encontrado:%.4f - Erro:%f - NumPassos:%d\n", esperado, resRk4, erroRk4, passosRk4);
		fprintf(output,"%.4f,%d,", resRk4, passosRk4);

		//Variavel
		passosRk5 = 0;
		resRk5 = RungeKuttaPassoVariavel(tolerancia, hRk5, theta, &passosRk5, 0);
		erroRk5 = fabs(resRk5 - esperado);
		printf("  RK5 - Esperado:%.4f - Encontrado:%.4f - Erro:%f - NumPassos:%d\n", esperado, resRk5, erroRk5, passosRk5);
		fprintf(output,"%.4f,%d,%f\n", resRk5, passosRk5, esperado);
		
		theta+=M_PI/40.0;
	}

	

	RungeKuttaPassoFixo(hRk4, 1.5707, &passosRk4,1);
	RungeKuttaPassoVariavel(tolerancia, hRk5, 1.5707, &passosRk5, 1);

	return 42;
}
