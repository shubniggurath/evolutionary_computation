#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define c 0.817
#define Gmax 1000000
#define error 0.00001


double aleatorio(int,int);
double  *vector(double *,int);
double inic_padre(double *,int);
double evaluarF(double*);
int restriccion(double *);

extern double ranf()
{
    return rand() / ((double) RAND_MAX);
}
double aleatorio(int min,int max){
       double scaled = (double)rand()/(double)RAND_MAX;
       return (max - min +1)*scaled + min;
}
double gauss(double m, double s)
{
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

double *vector(double *dato,int tam)
{
	dato = (double*)malloc(tam*sizeof(double));
}


double evaluarF(double *x)
{
	return ((0.6224*x[0]*x[2]*x[3])+(1.7781*x[1]*pow(x[2],2))+(3.1661*pow(x[0],2)*x[3])+((19.84*pow(x[0],2)*x[2])));
}

int restriccion(double *x)
{
	if ((((-1)*x[0])+(0.0193*x[2]))<=0 && (((-1)*x[1])+(0.00954*x[2]))<=0 && (x[3]-240)<=0 && (1296000-(((double)4/(double)3)*M_PI*pow(x[2],3))-(M_PI*pow(x[2],2)*x[3]))<=0)
		return 1;
	else
		return 0;
}

int main()
{
	double delta, *padre=NULL, *hijo=NULL, padreFunc, hijoFunc,num,aprox;
	int t, num_var,padre_sus,i,ps, ve,fact;
	t =0;
	ps =0;
	fact=0;
	num_var =4;
	delta =5; //explorar mas el espacio de busqueda reproduccion
	srand(time(0));
	padre =vector(padre,num_var);
	hijo=vector(hijo,num_var);
	aprox =1;
	for(i=0;i<num_var;i++)
	{
	    //num=delta*gauss(0,1);
	    if(i==0) {
                num=aleatorio(1,99)*0.0625;
             } else if (i==1) {
                  num=aleatorio(1,99)*0.0625;
             } else if (i==2) {
                 num=aleatorio(10,200);
            } else if (i==3) {
                num=aleatorio(10,200);
            }
		padre[i]=num;
	}

	/*padre[0]=0.826142;//0.778197751897;
	padre[1]=0.417884;//0.384665697936;
	padre[2]=42.712546;//40.321054550108;
	padre[3]=177.071810;//199.980236777701;*/
	padreFunc = evaluarF(padre);
    while((t<Gmax))
    //while((t<Gmax) && (aprox>error))
	{
		//printf("--- Iteracion Numero %d ---\n",t+1);
		for(i=0;i<num_var;i++){
            //num=padre[i]+delta*gauss(0,1);
            if(i==0) {
                num=aleatorio(1,99)*0.0625;
             } else if (i==1) {
                  num=aleatorio(1,99)*0.0625;
             } else if (i==2) {
                 num=aleatorio(10,200);
            } else if (i==3) {
                num=aleatorio(10,200);
            }
            hijo[i]=num;
        }
		hijoFunc = evaluarF(hijo); //evalua la funcion con los nuevos datos
		//printf("    Evaluacion padre: %f \n",padreFunc);
        //printf("    Evaluacion hijo: %f \n",hijoFunc);
		aprox=fabs(padreFunc-hijoFunc);
        //printf("    Error: %f \n",aprox);
		if (restriccion(hijo)==0 && restriccion(padre)==0)
		{
			//printf("NINGUNA ES FACTIBLE");
            fact++;
			if(hijoFunc<padreFunc)
			{
				//printf("  HIJO SUSTITUYO AL PADRE \n");
				for(i=0;i<num_var;i++)
				{
					padre[i]=hijo[i];
				}
				padreFunc=hijoFunc;
				ps++;
				ve++;
			}
			else
			{
				//printf("  SE MANTIENE EL PADRE\n");
			}
		}
		else if (restriccion(hijo)==1 && restriccion(padre)==0)
		{
			//printf("HIJO FACTIBLE");
			//printf("  HIJO SUSTITUYO AL PADRE \n");
			for(i=0;i<num_var;i++)
			{
				padre[i]=hijo[i];
			}
			padreFunc=hijoFunc;
			ps++;
			ve++;

		}
		else if (restriccion(hijo)==0 && restriccion(padre)==1)
		{
			//printf("PADRE FACTIBLE");
			//printf("  EL PADRE SE MANTIENE \n");

		}
		else if (restriccion(hijo)==1 && restriccion(padre)==1)
		{
			//printf("AMBAS SON FACTIBLES");

			if(hijoFunc<padreFunc)
			{
				//printf("  HIJO SUSTITUYO AL PADRE \n");
				for(i=0;i<num_var;i++)
				{
					padre[i]=hijo[i];
				}
				padreFunc=hijoFunc;
				ps++;
				ve++;
			}
			else
			{
				//printf("  SE MANTIENE EL PADRE\n");

			}
		}
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
		t++;
	}
	printf("\n---%d ---\n",fact);
	printf("\n------------RESULTADO---------------\n");
    for(i=0;i<num_var;i++)
        printf("x%d= %f \n",i+1,padre[i]);
    printf("f(Xi)=%f \n",padreFunc);
    printf("\n------------Condiciones---------------\n");
    printf("g1(x)=%f \n",((-1)*padre[0])+(0.0193*padre[2]));
    printf("g2(x)=%f \n",((-1)*padre[1])+(0.00954*padre[2]));
    printf("g3(x)=%f \n",-(M_PI*pow(padre[2],2)*padre[3])-(((double)4/(double)3)*M_PI*pow(padre[2],3))+1296000);
    printf("g4(x)=%f \n",padre[3]-240);
    return 0;
}
