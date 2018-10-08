#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define  TWO_PI 6.2831853071795864769252866
#define c 0.817
#define Gmax 5000
#define error 0.00001
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)
#define e1(x1,x3) (((-1)*x1)+(0.0193*x3)) 
#define e2(x2,x3) (((-1)*x2)+(0.0095*x3)) 
#define e3(x3,x4) (1296000-((4/3)*M_PI*pow(x3,3))-(M_PI*pow(x3,2)*x4))
#define e4(x4) (x4-240)


double aleatorio(double,double);
double  *vector(double *,int);
void inic_padre(double *,int);
double evaluarF(double*);
int restriccion(double *);
double genGausse(double, double);

double genGauss(double media, double desE)
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

double aleatorio(double min, double max)
{	
	double scaled = ranf();
       return (max - min)*scaled + min;
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
	if (e1(x[0],x[2])<=0 && e2(x[1],x[2])<=0 && e3(x[2],x[3])<=0 && e4(x[3])<=0)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void inic_padre(double * padre,int num_var)
{
	int i;
	for(i=0;i<num_var;i++)
	{	
		if((i==0) || (i==1))
		{
			padre[i]=aleatorio(1.0,99.0)*0.0625;
		}
		else if ((i==2) || (i==3))
		{
			padre[i]=aleatorio(10.0,200.0);
		}
		
	}
	
}



int main()
{
	double delta, *padre=NULL, *hijo=NULL, padreFunc, hijoFunc,num,aprox,*padres=NULL,**matr_Pad=NULL;
	int t, num_var,i,x,ps, ve=0,mu;
	t =0;
	ps =0;
	num_var =4;
	delta = 10; //explorar mas el espacio de busqueda reproduccion
	srand(time(0));
	hijo=vector(hijo,num_var);
	aprox =1;
	mu= 10;
	//formando el conjunto de padres
	padres = vector(padres,mu);
	matr_Pad=(double **) malloc (mu*sizeof(double*));
	
	for (x=0;x<mu;x++)
	{
		int y;
		double datosP;
		matr_Pad[x] = vector(matr_Pad[x],num_var);
		for(y=0;y<num_var;y++)
		{	
			if((y==0) || (y==1))
			{
				datosP=aleatorio(1.0,99.0)*0.0625;
			}
			else if ((y==2) || (y==3))
			{
				datosP=aleatorio(10.0,200.0);
			}
			
			matr_Pad[x][y] = datosP; //se evalua la funcion
		}
		
		padres[x]=evaluarF(matr_Pad[x]);
			
	}
	
	
    while((t<Gmax) && (aprox>error))
	{
		printf("--- Iteracion Numero %d ---\n",t+1);

		//evaluar soluciones
		int y, num_ale=0,num_ale2=0, x=0, pos=0,s,opt;
		for (x=0;x<mu;x++)
		{	
			for (y=0;y<num_var;y++)
			{	
				printf("%f ",matr_Pad[x][y]);
			}
			printf("%f\n",padres[x]);
		}
		
		/*se esscogen 2 sol aleatorias*/	
		do
		{
				
			num_ale = rand()%(10);
			num_ale2 = rand()%(10);
			if (num_ale ==num_ale2)
			{
				x=1;
			}
			else
			{
				x=0;
			}
				
		}
		while(x==1);
		//printf("--%d --%d\n",num_ale,num_ale2);
		
		for(i=0;i<num_var;i++){
		
			
             num=genGauss(0,11);//semilla de aleatorios
             //hijo[i]=padre[i]+delta*num; // mutando el vector
			 
			hijo[i]= ((matr_Pad[num_ale][i]+matr_Pad[num_ale2][i])/2)+delta*num; // mutando el vector
            /* printf("    Delta: %f  Aleatorio: %f \n",delta,num);
             printf("    Padre %d : %f \n",i+1,padre[i]);*/
             printf("    Hijo %d : %f \n",i+1,hijo[i]);
        }
		hijoFunc=evaluarF(hijo);
		
		printf("%f %d",hijoFunc,restriccion(hijo));
		
		for(s=1;s<mu;s++)
		{
			if (padres[0]<padres[s])
			{
				pos=s;
			}

        }
		
		for(s=1;s<mu;s++)
		{
			if (padres[0]>padres[s])
			{
				opt=s;
			}

        }
		//buscando al mejor de los padres
		printf("\n----%d  \n",opt);
		//cambiando por la peor solucion
		for(i=0;i<num_var;i++)
		{
			matr_Pad[pos][i]=hijo[i];
        }
		aprox=fabs(padres[opt]-hijoFunc);
		if(hijoFunc<padres[opt])
		{
			ps++;
		}
		padres[pos]=hijoFunc;
		printf("    Error: %f \n",aprox);
		
		/*printf("    Evaluacion padre: %f \n",padreFunc);
        printf("    Evaluacion hijo: %f \n",hijoFunc);
		aprox=fabs(padreFunc-hijoFunc);
        printf("    Error: %f \n",aprox);
		//printf("%f,%f",hijoFunc,padreFunc);*/
        
		/*if (restriccion(hijo)==0 && restriccion(padre)==0)
		{
			//printf("mejor  actitud");
			
			if(hijoFunc<padreFunc)
			{
				printf("  HIJO SUSTITUYO AL PADRE \n");
				for(i=0;i<num_var;i++)
				{
					matr_Pad[pos][i]=hijo[i];
				}
				padreFunc=hijoFunc;
				ps++;
				ve++;
			}
			else
			{
				printf("  SE MANTIENE EL PADRE\n");
			
			}
		}
		
		else if (restriccion(hijo)==1 && restriccion(padre)==0)
		{
			//printf("mejor  actitud");
			printf("  HIJO SUSTITUYO AL PADRE \n");
			for(i=0;i<num_var;i++)
			{
				matr_Pad[pos][i]=hijo[i];
			}
			padreFunc=hijoFunc;
			ps++;
			ve++;

		}
		
		else if (restriccion(hijo)==0 && restriccion(padre)==1)
		{
			//printf("mejor  actitud");
			printf("  EL PADRE SE MANTIENE \n");
			
		
		}
		
		else if (restriccion(hijo)==1 && restriccion(padre)==1)
		{
			//printf("mejor  actitud");
			
			if(hijoFunc<padreFunc)
			{
				printf("  HIJO SUSTITUYO AL PADRE \n");
				for(i=0;i<num_var;i++)
				{
					matr_Pad[pos][i]=hijo[i];
				}
				padreFunc=hijoFunc;
				ps++;
				ve++;
				
				
			}
			else
			{
				printf("  SE MANTIENE EL PADRE\n");
			
			}
		}*/
		
		
		
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

		
		t++;
	}
	printf("\n---%d",ve);
    return 0;
}
