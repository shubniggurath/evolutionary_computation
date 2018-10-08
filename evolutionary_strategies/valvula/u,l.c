#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define Gmax 10000
#define error 0.0000001  //le ponemos 4 ceros y genera mas rapido
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)
#define e1(x1,x3) (((-1)*(double)x1)+((double)0.0193*(double)x3)) 
#define e2(x2,x3) (((-1)*(double)x2)+((double)0.00954*(double)x3)) 
#define e3(x3,x4) (1296000-(((double)4/(double)3)*(double)M_PI*pow(x3,3))-(M_PI*pow(x3,2)*(double)x4))
#define e4(x4) ((double)x4-240)

typedef struct dat
{
	double x;
	double sigma;
} x_delt;


double  *vector(double *,int);
double evaluarF(x_delt*);
int restriccion(x_delt *);
double genGausse(double, double);
int buscarMpad(double *,int);
void aleatorios(int *, int);
void print(x_delt **,double *,int, int);
void printMS(x_delt **,double *,int, int);
void Fuv(x_delt **, x_delt **, x_delt **,double*, double *, double *, int,int,int);


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

double aleatorio(double min,double max){
       double scaled = (double)rand()/(double)RAND_MAX;
       return (max - min +1)*scaled + min;
}

double *vector(double *dato,int tam)
{
	dato = (double*)malloc(tam*sizeof(double));
}


double evaluarF(x_delt *x)
{
	return ((0.6224*x[0].x*x[2].x*x[3].x)+(1.7781*x[1].x*pow(x[2].x,2))+(3.1661*pow(x[0].x,2)*x[3].x)+((19.84*pow(x[0].x,2)*x[2].x)));
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

void printMS(x_delt **matr_Pad,double *evaluciones,int optimo, int tam_v)
{
	int y;
		printf("\nMejor resultado ---  ");
		for (y=0;y<tam_v;y++)
		{	
			printf("%f  ",matr_Pad[optimo][y].x);
		}
		printf(" f(x)=%f\n",evaluciones[optimo]);

}

void Fuv(x_delt **u, x_delt **v, x_delt **uv,double* evalTot, double *evaluciones, double *eval_H, int mu,int landa,int num_var)
{

	int x,s,tamT;
	tamT=mu+landa;
	
	for (x=0;x<mu;x++)
	{
		int y;
		uv[x] = (x_delt*)malloc(num_var*sizeof(x_delt));
		for(y=0;y<num_var;y++)
		{	
			uv[x][y].x = u[x][y].x;
			uv[x][y].sigma = u[x][y].sigma;
		}
		
		 evalTot[x] =  evaluciones[x];
		
	}
	
	for (x=0,s=mu;x<landa,s<tamT;x++,s++)
	{
		int y;
		uv[s] = (x_delt*)malloc(num_var*sizeof(x_delt));
		for(y=0;y<num_var;y++)
		{	
			uv[s][y].x = v[x][y].x;
			uv[s][y].sigma = v[x][y].sigma;
		}
		
		 evalTot[s] =  eval_H[x];
			
	}

}

void aleatorios(int *num_ale,int mu)
{
	int x,num_1,num_2;
	do
		{
				
			num_1 = rand()%(mu); 
			num_2 = rand()%(mu);
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

void print(x_delt **matr_Pad,double *evaluciones,int n_d, int tam_v)
{
	int x,y;
		for (x=0;x<n_d;x++)
		{	
			printf("%d " ,(x+1));
			for (y=0;y<tam_v;y++)
			{	
				printf("%f  ",matr_Pad[x][y].x);
			}
			
			printf("f(x)=%f ",evaluciones[x]);
			printf(" res: %d \n",restriccion(matr_Pad[x]));
		}



}

int restriccion(x_delt *x)
{
	if (e1(x[0].x,x[2].x)<=0 && e2(x[1].x,x[2].x)<=0 && e3(x[2].x,x[3].x)<=0 && e4(x[3].x)<=0)
	{	
		return 1;
	}
	else
	{
		return 0;
	}
}

void ordenar(x_delt **uv, double *evalTot,int tamT,int num_var)
{
	int x;
	double *vecEval;	
	
	for(x=1;x<tamT;x++)
	{
		int min,z;
		double aux;
		min =x;
		while(min >0 && evalTot[min]<evalTot[min-1])
		{
			aux = evalTot[min];
			evalTot[min] =evalTot[min-1];
			evalTot[min-1] =aux;
			
			for(z=0;z<num_var;z++)
			{
				double v_a=uv[min][z].x;
				double del=uv[min][z].sigma;
				
				uv[min][z].x=uv[min-1][z].x;
				uv[min-1][z].x=v_a;
				
				uv[min][z].sigma=uv[min-1][z].sigma;
				uv[min-1][z].sigma=del;
			
			}
			min--;
		}
		
			
	}


}

void escoger_validos(x_delt **mat, double *evalTot,int tamT,int num_var)
{
	x_delt **aux,**aux2=NULL;
	double *evaluciones=NULL, *evaluciones2= NULL;
	int u,tamF=0,tam2=0, tam=tamT, conteo=0;
	
	evaluciones= vector(evaluciones,tamT);
	evaluciones2= vector(evaluciones,tamT);
	aux = (x_delt **) malloc (tamT*sizeof(x_delt*));
	aux2 = (x_delt **) malloc (tamT*sizeof(x_delt*));
	
	//printf("----------- %d ",tamT);
	
	for(u =0;u<tamT;u++)
	{
		if (restriccion(mat[u])== 1)
		{	
			int i;
			aux[tamF] = (x_delt *) malloc (num_var*sizeof(x_delt));
			
			for (i =0; i< num_var; i++)
			{
				aux[tamF][i] =mat[u][i];
				
			}
			evaluciones[tamF] = evalTot[u];
			tamF ++;
		}
		else
		{
			aux2[tam2] = (x_delt *) malloc (num_var*sizeof(x_delt));
			int i;
			
			for (i =0; i< num_var; i++)
			{
				aux2[tam2][i] =mat[u][i];
				
			}
			evaluciones2[tam2] = evalTot[u];
			
			tam2 ++;
		}
	}
	//print(aux, evaluciones, tamF, num_var);
	//printf("\n ORDENADOS\n");
	ordenar(aux, evaluciones,tamF,num_var);
	//print(aux, evaluciones, tamF, num_var);
	//printf("\n -------\n");
	//print(aux2, evaluciones2, tam2, num_var);
	//printf("\n ORDENADOS\n");
	ordenar(aux2, evaluciones2,tam2,num_var);
	//print(aux2, evaluciones2, tam2, num_var);*/
	
	//printf(" numero de valores que cumplen: %d ",tamF);
	Fuv(aux, aux2, mat,evalTot, evaluciones, evaluciones2, tamF, tam2, num_var);


	//print(mat, evalTot, tamT, num_var);
	
	
	/*while(tam > 0)
	{
		for(u =0;u<tamF;u++)
		{		
			int i;
			for (i =0; i< num_var; i++)
			{
				mat[conteo][i] = aux[u][i];
			}
			
			evalTot[conteo]= evaluciones[u];
			conteo++;
		}
		
		tam -= tamF;
		if (tam<=tamF)
		{
			tamF = tam;
		}		
		
	}*/
	//free(evaluciones);
	//free(aux);


}

void imprimirG(x_delt *x)
{
	printf("g1 = %f \n g2 = %f \n g3 = %f \n g4 = %f \n",e1(x[0].x,x[2].x),e2(x[1].x,x[2].x), e3(x[2].x,x[3].x),e4(x[3].x));

}

int main()
{
	double **sigma=NULL, *padre=NULL, padreFunc, hijoFunc,num,aprox,*evaluciones=NULL,*eval_H=NULL,tau,tau_p;

	int t, num_var,i,x,ps,mu,landa,mejorS;
	t =0;
	ps =0;
	num_var =4;	
	srand(time(NULL));
	aprox =1;
	mu= 10;
	landa = 70;
	//sigma = (double **) malloc (landa*sizeof(double*)); //explorar mas el espacio de busqueda reproduccion
	//hijos=(double **) malloc (landa*sizeof(double*));
	//formando el conjunto de evaluciones
	evaluciones= vector(evaluciones,mu);
	eval_H = vector(eval_H,landa);
	//matr_Pad=(double **) malloc (mu*sizeof(double*));
	
	x_delt **matr_Pad, **hijos;
		
	matr_Pad = (x_delt **) malloc (mu*sizeof(x_delt*));
	hijos = (x_delt **) malloc (landa*sizeof(x_delt*));
	
	//inicilaizar padre
	for (x=0;x<mu;x++)
	{
		int y;
		double datosP;
		matr_Pad[x] = (x_delt*)malloc(num_var*sizeof(x_delt));
		
		for(y=0;y<num_var;y++)
		{	
			if((y==0) || (y==1))
			{
				datosP=aleatorio(1.0,90.0)*0.0625;
				matr_Pad[x][y].x = datosP; //se evalua la funcion
				matr_Pad[x][y].sigma = 10;
			}
			else if ((y==2) || (y==3))
			{	
				datosP=aleatorio(10.0,200.0);
				matr_Pad[x][y].x = datosP; //se evalua la funcion
				matr_Pad[x][y].sigma = 10;
				//matr_Pad[x][y].sigma = genGauss(10.0,200.0)*0.9;

			}
	
		}
		//printf(" %d\n",restriccion(matr_Pad[x]));
		evaluciones[x]=evaluarF(matr_Pad[x]);
			
	}
	
	mejorS =buscarMpad(evaluciones,mu);
	//print(matr_Pad,evaluciones,mu,num_var);
	//printMS(matr_Pad,evaluciones,mejorS,num_var);
	tau=(1/(4.0*pow(num_var,0.25)));
	tau_p=(1/(pow(2.0*(num_var),0.5)));
	
	while((t<Gmax) && (aprox>error))
	{
		int y,num_ale[2],s,m,nu;
		double M_hijo=0,var=0.0;
		//printf("\n--- Iteracion Numero %d ---\n",t+1);
		
		//print(matr_Pad,evaluciones,mu,num_var);
		//printf("\n");
			
		
		
		for(s=0;s<landa;s++)
		{	
			/*se esscogen 2 sol aleatorias que no se repiten, LA MU SE PUEDE PONER EN 10 y trabaja mejor*/
            aleatorios(num_ale,mu);
			//printf("--%d --%d\n",num_ale[0],num_ale[1]);
			hijos[s] = (x_delt*)malloc(num_var*sizeof(x_delt));
			nu =genGauss(0.0,1.0);
			
			for(i=0;i<num_var;i++)
			{		
				double x,sigma;
				num=genGauss(0.0,1.0);//semilla de aleatorios
				//aplicando recombinacion plana
				//h = ((num*matr_Pad[num_ale[0]][i].x)+((1-num)*matr_Pad[num_ale[1]][i].x))*0.0817;
				x = matr_Pad[num_ale[1]][i].x;
				sigma =	matr_Pad[num_ale[1]][i].sigma;
				//h = (matr_Pad[num_ale[0]][i]+matr_Pad[num_ale[1]][i])/2; 			
				//printf("//%f",h);
				//aplicando mutacion
				hijos[s][i].sigma=sigma*(exp((tau*num)+(tau_p*nu)));
				//printf("********** %f\n",sigma[s][i]);
				//num =genGauss(0.0,pow(d,2));
				x=x+hijos[s][i].sigma*num;
				//printf(" -- %f  \n",h);
				hijos[s][i].x=x;
			}
			
			eval_H[s]=evaluarF(hijos[s]);
			//printf("%f\n",evaluarF(hijos[s]));
			//printf(" %d\n",restriccion(hijos[s]));
		}

		//print(hijos,eval_H,landa,num_var);
		m=mejorS;
		mejorS =buscarMpad(eval_H,landa);
		//printMS(hijos,eval_H,mejorS,num_var);
		//print(hijos,eval_H,landa,num_var);
		aprox=fabs(evaluciones[m]-eval_H[mejorS]);
		//printf("Error: %f \n",aprox);
		
		escoger_validos(hijos, eval_H,landa,num_var);
		//ordenamiento por seleccion
		//ordenar(hijos, eval_H,landa,num_var);
		
		//print(hijos,eval_H,landa,num_var);
		
		//sustituyecdo a los padres
		for (x=0;x<mu;x++)
		{
			int y;
			for(y=0;y<num_var;y++)
			{					
				matr_Pad[x][y].x = hijos[x][y].x;
				matr_Pad[x][y].sigma = hijos[x][y].sigma;
			}
			evaluciones[x] = eval_H[x];
		}
		
		t++;
		
	}
	
	printf("\n");
	print(matr_Pad, evaluciones, mu, num_var);
	printf("\n--- Iteracion Numero %d ---\n",t+1);
	int y;				
	for (y=0;y<num_var;y++){	
		printf(" x%d:%f ",y+1,matr_Pad[0][y].x);
	}
	printf("\n f(x): %f  \n",evaluciones[0]);
	imprimirG(matr_Pad[0]);
	/*mejorS =buscarMpad(evaluciones,mu);
	printMS(matr_Pad,evaluciones,mejorS,num_var);
	imprimirG(matr_Pad[mejorS]);*/
	printf("Error: %f \n",aprox);
	

    return 0;
}