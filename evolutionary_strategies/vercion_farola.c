#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define  TWO_PI 6.2831853071795864769252866
#define c 0.817
#define Gmax 1000
#define error 0.000001  //le ponemos 4 ceros y genera mas rapido
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)
#define codo 0.52355866484
#define aureo (1+sqrt(5))/2


double aleatorio(double,double);
double  *vector(double *,int);
double evaluarF(double*);
int restriccion(double *);
double genGausse(double, double);
int buscarMpad(double *,int);
void aleatorios(int *);
void print(double **m,double *,int, int);

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
	return pow((10*pow(x[0],3)+3*pow(x[0],2)+5),2);
	//return pow((pow(x[0],2)-1),3) *(pow((2*x[0]-5),4));
	//return 3*pow(x[0],4)+pow(x[0],2)-2*x[0]+1; 
	//return pow((x[0]+10*x[1]),2)+5*pow((x[2]-x[3]),2)+pow(x[1]-2*x[2],4)+10*pow((x[0]-x[3]),4);
	//return x[0]*x[0]+1000;
}


int buscarMpad(double *evaluciones,int mu)
{
	
	double var=0.0;
	int opt=0,k;
	var=evaluciones[0];
	for(k=0;k<mu;k++)
	{
		if (var>evaluciones[k])
		{
			var=evaluciones[k];
			opt=k;
		}
    }
	//printf("\nmejor solucion---- %d  \n",(opt+1));
	return  opt;
}


void aleatorios(int *num_ale)
{
	int x,num_1,num_2;
	do
		{
				
			num_1 = rand()%(10); 
			num_2 = rand()%(10);
			if (num_1 ==num_2)
			{
				x=1;
			}
			else
			{
				x=0;
			}
				
		}
	while(x==1);
	num_ale[0]=num_1;
	num_ale[1]=num_2;
}

void print(double **matr_Pad,double *evaluciones,int n_d, int tam_v)
{
	int x,y;
		for (x=0;x<n_d;x++)
		{	
			for (y=0;y<tam_v;y++)
			{	
				printf("%d   %f ",(x+1),matr_Pad[x][y]);
			}
			printf("    f(x)=%f\n",evaluciones[x]);
		}


}
int main()
{
	double **delta=NULL, *padre=NULL, **hijos=NULL, padreFunc, hijoFunc,num,aprox,*evaluciones=NULL,*eval_H=NULL,**matr_Pad=NULL;
	int t, num_var,i,x,ps,mu,landa,mejorS;
	t =0;
	ps =0;
	num_var =2;	
	srand(time(NULL));
	aprox =1;
	mu= 20;
	landa =30;
	delta = (double **) malloc (landa*sizeof(double*)); //explorar mas el espacio de busqueda reproduccion
	hijos=(double **) malloc (landa*sizeof(double*));
	//formando el conjunto de evaluciones
	evaluciones= vector(evaluciones,mu);
	eval_H = vector(eval_H,landa);
	matr_Pad=(double **) malloc (mu*sizeof(double*));
	
	
	//inicilaizar padre
	for (x=0;x<mu;x++)
	{
		int y;
		double datosP;
		matr_Pad[x] = vector(matr_Pad[x],num_var);
		for(y=0;y<num_var;y++)
		{	
			
			datosP=genGauss(0.0,1.0);
			matr_Pad[x][y] = datosP; //se evalua la funcion
		}
		
		evaluciones[x]=evaluarF(matr_Pad[x]);
			
	}
	
	for (x=0;x<landa;x++)
	{
		int y;
		delta[x]=vector(delta[x],num_var);
		for(y=0;y<num_var;y++)
		{	
			//puedes poner 10 y tarda un poco mas  
			delta[x][y] = 0.0; 
		}
			
	}
	
	
	
	mejorS =buscarMpad(evaluciones,mu);
	printf("\nMejor resultado --- %f    f(x)=%f\n",matr_Pad[mejorS][0],evaluciones[mejorS]);
	
    while((t<Gmax) && (aprox>error))
	{
		int y,num_ale[2],s,m,nu;
		double M_hijo=0,var=0.0,tau,tau_p;
		printf("\n--- Iteracion Numero %d ---",t+1);
		
		//print(matr_Pad,evaluciones,mu,num_var);
		printf("\n");
		nu =genGauss(0.0,1.0);
		
		tau=(1/sqrt(4.0*pow(num_var,0.25)));
		tau_p=(1/sqrt(2.0*(num_var)));
		
		
		for(s=0;s<landa;s++)
		{	
			/*se esscogen 2 sol aleatorias que no se repiten*/
            aleatorios(num_ale);
			//printf("--%d --%d\n",num_ale[0],num_ale[1]);
			hijos[s] = vector(hijos[s],num_var);
			
			for(i=0;i<num_var;i++)
			{		
				double h,d;
				num=genGauss(0.0,1.0);//semilla de aleatorios
				//aplicando recombinacion plana
				h = ((num*matr_Pad[num_ale[0]][i])+((1-num)*matr_Pad[num_ale[1]][i]));
				d =	(delta[num_ale[0]][i]+delta[num_ale[0]][i])/2;
				//h = (matr_Pad[num_ale[0]][i]+matr_Pad[num_ale[1]][i])/2; 			
				//printf("//%f",h);
				//aplicando mutacion
				delta[s][i]=d*(exp((tau*num)+(tau_p*nu)));
				//printf("********** %f\n",delta[s][i]);
				h=h+delta[s][i]*num;
				//printf(" -- %f  \n",h);
				hijos[s][i]=h;
			}

			eval_H[s]=evaluarF(hijos[s]);
			//printf("%f\n",evaluarF(hijos[s]));
		}
		//print(hijos,eval_H,landa,num_var);
		m=mejorS;
		mejorS =buscarMpad(eval_H,landa);
		printf("\nMejor resultado --- %f    f(x)=%f\n",hijos[mejorS][0],eval_H[mejorS]);
		aprox=fabs(evaluciones[m]-eval_H[mejorS]);
		printf("Error: %f \n",aprox);
		
		
		//ordenamiento por seleccion
		for(x=0;x<landa;x++)
		{
			int y,min,z;
			double aux;
			min =x;
			for(y=(x+1);y<landa;y++)
			{
				if (eval_H[y]<eval_H[min])
					min=y;
			}
			aux = eval_H[x];
			eval_H[x] =eval_H[min];
			eval_H[min] =aux;
			
			for(z=0;z<num_var;z++)
			{
				double v_a=hijos[x][z];
				double del=delta[x][z];
				
				hijos[x][z]=hijos[min][z];
				hijos[min][z]=v_a;
				
				delta[x][z]=delta[min][z];
				delta[min][z]=del;
			
			}
        }
		
		//print(hijos,eval_H,landa,num_var);
		
		//sustituyecdo a los padres
		for (x=0;x<mu;x++)
		{
			int y;
			for(y=0;y<num_var;y++)
			{					
				matr_Pad[x][y] = hijos[x][y];
			}
			evaluciones[x] = eval_H[x];
		}
		
		t++;
	}

    return 0;
}