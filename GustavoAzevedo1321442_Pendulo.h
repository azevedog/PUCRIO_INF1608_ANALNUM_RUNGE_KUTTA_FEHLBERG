/**
Gustavo Azevedo - 1321442
Trabalho Final Análise Numérica (INF-1608)
"Movimento de um Pêndulo"

API
*/

/** funcao de configuracao dos parametros do pendulo.
@PARAM gravidade - Aceleracao da gravidade em m/s^2.
@PARAM comprimento - Comprimento em metros do haste do pendulo.
 */
void configuraPendulo(double gravidade, double comprimento);

/** Calcula o periodo em segundos. */
double periodoEsperadoSimplificado();

/** Calcula o periodo em segundos. */
double periodoEsperado(double theta0);

/* Calcula o periodo de um pendulo configurado ou do pendulo padrão, com passo VARIAVEL, usando metodo acoplado de Runge-Kutta de ordem 4/5.

@PARAM Tolerancia - erro local do passo adaptativo.
@PARAM h - Sugestao de passo inicial. Usado para definir hMin(h/10.0) e hMax(h*10.0).
@PARAM theta0 - Angulo inicial. Posicao de onde o pendulo e solto.
@PARAM numPassos - Numero de passos dado para chegar ao resultado.
@PARAM print - Define se as posicoes passo a passo devem ser impresas em um arquivo. 0 para nao e 1 para sim.

@RETURN - Periodo em segundos.
*/
double RungeKuttaPassoVariavel(double tolerancia, double h, double theta0, int* numPassos, int print);

/* Calcula o periodo de um pendulo configurado ou do pendulo padrão, com passo FIXO, usando metodo de Runge-Kutta de ordem 4. Esse método foi usado para estabelecer a implementacao no inicio do trabalho para validacao e para comparacao de performance em numero de passos.

@PARAM h - Sugestao de passo inicial. Usado para definir hMin(h/10.0) e hMax(h*10.0).
@PARAM theta0 - Angulo inicial. Posicao de onde o pendulo e solto.
@PARAM numPassos - Numero de passos dado para chegar ao resultado.
@PARAM print - Define se as posicoes passo a passo devem ser impresas em um arquivo. 0 para nao e 1 para sim.

@RETURN - Periodo em segundos.
*/
double RungeKuttaPassoFixo(double h, double theta0, int* numPassos, int print);
