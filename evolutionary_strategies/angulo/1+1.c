#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define  TWO_PI 6.2831853071795864769252866
#define c 0.817
#define Gmax 500000
#define error 0.00001
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)
#define e1(x1,x3) (((-1)*(double)x1)+((double)0.0193*(double)x3)) 
#define e2(x2,x3) (((-1)*(double)x2)+((double)0.00954*(double)x3)) 
#define e3(x3,x4) (1296000-(((double)4/(double)3)*(double)M_PI*pow(x3,3))-(M_PI*pow(x3,2)*(double)x4))
#define e4(x4) ((double)x4-240)


//double aleatorio(double,double);
double aleatorio(double,double);
double  *vector(double *,int);
double inic_padre(double);
double evaluarF(double);
int restriccion(double );
double generateGaussianNoise(const double);
double aleatorioGauss(double, double);
float g= 0.0;

double generateGaussianNoise(const double variance)
{
	static int haveSpare = 0;
	static double rand1, rand2;
 
	if(haveSpare)
	{
		haveSpare = 0;
		return sqrt(variance * rand1) * sin(rand2);
	}
 
	haveSpare = 1;
 
	rand1 = ranf();
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = ranf() * TWO_PI;
 
	return sqrt(variance * rand1) * cos(rand2);
}

/*double aleatorio(double min, double max)
{	
	double scaled = ranf();
       return (max - min)*scaled + min;
}*/

double aleatorioGauss(double media, double desE)
{	
	double x1,x2,w,y1;
	static double y2;
	static int utiliza_ult=0;
	
	if(utiliza_ult)
	{
		y1=y2;
		utiliza_ult =0;
	}
	else
	{
		do{
			x1 = 2.0 * ranf() - 1.0;
            x2 = 2.0 * ranf() - 1.0;
            w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

			w = sqrt( (-2.0 * log( w ) ) / w );
			y1 = x1 * w;
			y2 = x2 * w;
			utiliza_ult =1;
	}
	
	return (media + y1 * desE);
}

/*double aleatorio(const double variance)
{
	static int haveSpare = 0;
	static double rand1, rand2;
 
	if(haveSpare == 1)
	{
		haveSpare = 0;
		return sqrt(variance * rand1) * sin(rand2);
	}
 
	haveSpare = 1;
 
	rand1 = rand() / ((double) RAND_MAX);
	if(rand1 < 1e-100) rand1 = 1e-100;
	rand1 = -2 * log(rand1);
	rand2 = (rand() / ((double) RAND_MAX)) * TWO_PI;
 
	return sqrt(variance * rand1) * cos(rand2);
}*/

double aleatorio(double min,double max){
       double scaled = (double)rand()/(double)RAND_MAX;
       return (max - min +1)*scaled + min;
}


double *vector(double *dato,int tam)
{
	dato = (double*)malloc(tam*sizeof(double));
}


double evaluarF(double alfa)
{
	int l =2, h=1, b1=15;
	float A,B,C,D=0.7,E;
	A = l*sin(b1);
	B = l*cos(b1);
	C = ((h + 0.5*D)*sin(b1)) -(0.5*D*tan(b1));
	E = ((h+0.5*D)*cos(b1))-(0.5*D);
	g = ((A*sin(alfa)*cos(alfa))+(B*pow(sin(alfa),2))-(C*cos(alfa))-(E*sin(alfa)));
	return 1000000 -(fabs(g));
}




double inic_padre(double padre)
{

	padre = aleatorio(0.0,90.0);
	
	return padre;
}



int main()
{
	double delta, padre=0.0, hijo=0.0, padreFunc, hijoFunc,num,aprox,*padres=NULL;
	int t, num_var,padre_sus,i,x,ps, ve=0;
	t =0;
	ps =0;
	num_var =1;
	delta = 10; //explorar mas el espacio de busqueda reproduccion
	srand(time(0));
	aprox =1;

	
    padre = inic_padre(padre); //se evalua la funcion
    padreFunc = evaluarF(padre);
	
	
    while((t<Gmax) && (aprox>error))
	{
		//printf("--- Iteracion Numero %d ---\n",t+1);
		//num = aleatorio(0.0,90.0);
	
		hijo= padre+delta*aleatorioGauss(0,1);
		//hijo= num; // mutando el vector;

		hijoFunc = evaluarF(hijo); //evalua la funcion con los nuevos datos
		/*printf("    Evaluacion padre: %f \n",padreFunc);
        printf("    Evaluacion hijo: %f \n",hijoFunc);*/
		aprox=fabs(padreFunc-hijoFunc);
        //printf("    Error: %f \n",aprox);
		//printf("%f,%f",hijoFunc,padreFunc);
        	
		if(hijoFunc>padreFunc)
		{
			
			padre=hijo;	
			padreFunc=hijoFunc;
			ps++;
			ve++;
		}

		
		
		//usando el 1/5
		if((t%(10*num_var)) == 0){
            if(ps<(2*num_var)){
                delta= delta*c;
            }else if(ps>(2*num_var)){
                delta=delta/c;
            }else{
				delta=delta;              
            }
            ps=0;
        }
		/*else
		{
			delta =
		}*/
		
		t++;
	}
	
    printf("    Alfa: %f\n g(alfa) = %f \n f(alfa) = %f",padre,g,padreFunc);
 
	
	printf("\n---%d",ve);
    return 0;
}
