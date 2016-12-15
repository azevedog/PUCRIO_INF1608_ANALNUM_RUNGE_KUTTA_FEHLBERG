/**
Gustavo Azevedo - 1321442
Trabalho Final Análise Numérica (INF-1608)
"Movimento de um Pêndulo"

IMPLEMENTACAO
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "GustavoAzevedo1321442_Pendulo.h"


/** Gravidade (m/s^2). */
static double G = 9.81;
/** Comprimento em metros do haste do pendulo. */
static double L = 1.0;

/** funcao de configuracao dos parametros do pendulo.
@PARAM gravidade - Aceleracao da gravidade em m/s^2.
@PARAM comprimento - Comprimento em metros do haste do pendulo.
 */
//@OVERRIDE - Implementacao de pendulo.h
void configuraPendulo(double gravidade, double comprimento){
	G= gravidade;
	L= comprimento;
}
 
/** Calcula o periodo em segundos. */
//@OVERRIDE - Implementacao de pendulo.h
double periodoEsperadoSimplificado(){
	return 2*M_PI*sqrt(L/G);
}

/** Calcula o periodo em segundos. */
//@OVERRIDE - Implementacao de pendulo.h
double periodoEsperado(double theta0){
	return 2*M_PI*sqrt(L/G)*(1+(pow(theta0,2)/16.0)+((11*pow(theta0,4))/3072.0)+((173*pow(theta0,6))/737280.0)+((22931*pow(theta0,8))/1321205760.0));
}

/** Calcula a posicao (angulo) para valores baixos de theta. */
double resultadoAnalitico(double theta0, double t){
	return theta0*cos(sqrt((G*t)/L));
}

/** Funcao aceleracao simplificada para baixos valores de theta.
resposta é a aceleracao do pendulo em segundos^-2.*/
double AceleracaoSimplificada(double theta){
	return -(G/L)*theta;
}

/* Calcula o periodo de um pendulo configurado ou do pendulo padrão, com passo FIXO, usando metodo de Runge-Kutta de ordem 4. Esse método foi usado para estabelecer a implementacao no inicio do trabalho para validacao e para comparacao de performance em numero de passos.

@PARAM h - Sugestao de passo inicial. Usado para definir hMin(h/10.0) e hMax(h*10.0).
@PARAM theta0 - Angulo inicial. Posicao de onde o pendulo e solto.
@PARAM numPassos - Numero de passos dado para chegar ao resultado.
@PARAM print - Define se as posicoes passo a passo devem ser impresas em um arquivo. 0 para nao e 1 para sim.

@RETURN - Periodo em segundos.
*/
//@OVERRIDE - Implementacao de pendulo.h
double RungeKuttaPassoFixo(double h, double theta0,  int* numPassos, int print){
	double period = 0.0;
	double theta=theta0, dTheta1, dTheta2, dTheta3, dTheta4;
	double v = 0.0, dV1, dV2, dV3, dV4;
	FILE* output;

	if(print) {
		output = fopen("RK_FIXO.csv","w");
		if(output==NULL) return 2;
		fprintf(output,"Theta: %.6f\n", theta0);
		fprintf(output, "periodo,theta_Calc,theta_Anali\n");
	}

	(*numPassos) = 0;

	/* O uso do periodo na avaliacao de parada eh necessario para garantir que nao pare no inicio do percorrimeto, quando ainda esta proximo do valor inicial.*/
	while(fabs(theta - theta0)>h || (*numPassos)<1000){
		dTheta1=h*v;
		dV1=h*AceleracaoSimplificada(theta);
		dTheta2=h*(v+(dV1/2.0));
		dV2=h*AceleracaoSimplificada(theta+(dTheta1/2.0));
		dTheta3=h*(v+(dV2/2.0));
		dV3=h*AceleracaoSimplificada(theta+(dTheta2/2.0));
		dTheta4=h*(v+dV3);
		dV4=h*AceleracaoSimplificada(theta+dTheta3);
		theta += ((dTheta1+(2*dTheta2)+(2*dTheta3)+dTheta4)/6.0);
		v += ((dV1+(2*dV2)+(2*dV3)+dV4)/6.0);
		period+=h;
		(*numPassos) = (*numPassos)+1;
		if(print) fprintf(output, "%.6f,%.6f,%.6f\n", period, theta, resultadoAnalitico(theta0, period));
	}
	
	if(print) fclose(output);
	return period;
}

double rkStepsPos(double theta,double* dTheta, double v, double* dV, double hMin, double* hAtual, double hMax, double eMin, double eMax){
	double aux_dTheta, aux_dTheta_1, dTheta1, dTheta2, dTheta3, dTheta4, dTheta5 ,dTheta6;
	double aux_dV_1, dV1, dV2, dV3, dV4, dV5, dV6;
	double thetaError;
	double h = *hAtual;
	
	if(h<hMin){
		h=hMin;
	}else if(h>hMax){
		h=hMax;
	}

	dTheta1=h*v;
	dV1=h*AceleracaoSimplificada(theta);
	dTheta2=h*(v+(dV1/4.0));
	dV2=h*AceleracaoSimplificada(theta+(dTheta1/4.0));
	dTheta3=h*(v+(3*dV1/32.0)+(9*dV2/32.0));
	dV3=h*AceleracaoSimplificada(theta+(3*dTheta1/32.0)+(9*dTheta2/32.0));
	dTheta4=h*(v+(1932*dV1/2197.0)-(7200*dV2/2197.0)+(7296*dV3/2197.0));
	dV4=h*AceleracaoSimplificada(theta+(1932*dTheta1/2197.0)-(7200*dTheta2/2197.0)+(7296*dTheta3/2197.0));
	dTheta5=h*(v+(439*dV1/216.0)-(8*dV2)+(3680*dV3/513.0)-(845*dV4/4104.0));
	dV5=h*AceleracaoSimplificada(theta+(439*dTheta1/216.0)-(8*dTheta2)+(3680*dTheta3/513.0)-(845*dTheta4/4104.0));
	dTheta6=h*(v-(8*dV1/27.0)+(2*dV2)-(3544*dV3/2565.0)+(1859*dV4/4104.0)-(11*dV5/40.0));
	dV6=h*AceleracaoSimplificada(theta-(8*dTheta1/27.0)+(2*dTheta2)-(3544*dTheta3/2565.0)+(1859*dTheta4/4104.0)-(11*dTheta5/40.0));

	aux_dTheta = (25*dTheta1/216.0)+(1408*dTheta3/2565.0)+(2197*dTheta4/4104.0)-(dTheta5/5.0);
	aux_dTheta_1 = (16*dTheta1/135.0)+(6656*dTheta3/12825.0)+(28561*dTheta4/56430.0)-(9*dTheta5/50.0)+(2*dTheta6/55.0);
	
	aux_dV_1 = (16*dV1/135.0)+(6656*dV3/12825.0)+(28561*dV4/56430.0)-(9*dV5/50.0)+(2*dV6/55.0);


	thetaError = fabs(fabs(aux_dTheta_1) - fabs(aux_dTheta))/fabs(theta+aux_dTheta_1);
	
	if(thetaError>eMax && h>hMin){
		*hAtual = h/2.0;
		return rkStepsPos(theta, dTheta, v, dV, hMin, hAtual, hMax, eMin, eMax);
	}else{
		if(thetaError<eMin){
			*hAtual = h*2.0;
		}
		*dTheta	= aux_dTheta_1;
		*dV = aux_dV_1;
		return h;
	}
}

/** Funcao de RungeKutta com passo fixo, usada no inicio do trabalho para validacao da implementacao.*/
/* Calcula o periodo de um pendulo configurado ou do pendulo padrão, com passo VARIAVEL, usando metodo acoplado de Runge-Kutta de ordem 4/5.

@PARAM Tolerancia - erro local do passo adaptativo.
@PARAM h - Sugestao de passo inicial. Usado para definir hMin(h/10.0) e hMax(h*10.0).
@PARAM theta0 - Angulo inicial. Posicao de onde o pendulo e solto.
@PARAM numPassos - Numero de passos dado para chegar ao resultado.
@PARAM print - Define se as posicoes passo a passo devem ser impresas em um arquivo. 0 para nao e 1 para sim.

@RETURN - Periodo em segundos.
*/
//@OVERRIDE - Implementacao de pendulo.h
double RungeKuttaPassoVariavel(double tolerancia, double h, double theta0, int* numPassos, int print){
	
	double period = 0.0;
	double theta = theta0, dTheta;
	double v = 0.0, dV;
	double hMin = h/10.0, hMax = h*10.0;
	double eMin = tolerancia/10.0, eMax=tolerancia;
	double usedH;
	FILE* output;

	if(print){
		output = fopen("RK_VARIAVEL.csv","w");
		if(output==NULL) return 2;

		fprintf(output,"Theta: %.6f\n",theta0);
		fprintf(output, "periodo,theta_Calc,theta_Anali,Passo\n");
	}

	*numPassos = 0;
	
	/* O uso do número de passos na avaliacao de parada eh necessario para garantir que nao pare no inicio do percorrimeto, quando ainda esta proximo do valor inicial.*/
	while(fabs(theta - theta0)>0.001 || (*numPassos)<10){
		usedH = rkStepsPos(theta, &dTheta, v, &dV, hMin, &h, hMax, eMin, eMax);
		theta += dTheta;
		v += dV;
		period+=usedH;
		(*numPassos) = (*numPassos)+1;
		if(print) fprintf(output, "%.6f,%.6f,%.6f,%.6f\n", period, theta, resultadoAnalitico(theta0, period), usedH);
	}

	if(print) fclose(output);
	return period;
}
