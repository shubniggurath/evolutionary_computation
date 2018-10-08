#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define Gmax 100000
#define error 0.000001  //le ponemos 4 ceros y genera mas rapido
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)

typedef struct dat
{
	double x;
	double sigma;
} x_delt;



double  *vector(double *,int);
double evaluarF(x_delt*);
int restriccion(double *);
double genGausse(double, double);
int buscarMpad(double *,int);
void aleatorios(int *,int);
void print(x_delt **,double *,int, int);
void Fuv(x_delt **, x_delt **,x_delt **,double*,double *,double *,int,int,int);
void ordenar(x_delt **,double *,int,int);
void printMS(x_delt **,double *,int, int);
double aleatorio(double, double);


double aleatorio(double min, double max)
{	
	double scaled = ranf();
       return (max - min)*scaled + min;
}

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


double *vector(double *dato,int tam)
{
	dato = (double*)malloc(tam*sizeof(double));
}


double evaluarF(x_delt *x)
{
	int i,j,k=500;
	double res=0.0;
	int a[2][25]={{-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32},{-32,-32,-32,-32,-32,-16,-16,-16,-16,-16,0,0,0,0,0,16,16,16,16,16,32,32,32,32,32}};
	//return (-0.75/(1+pow(x[0].x,2))) - (0.65*x[0].x*atan(1/x[0].x))+0.65;
	//return ((-4*pow(x[0].x,2)) -20*x[0].x -100) + pow(1-x[0].x,4);
	//return 3*pow(x[0].x,2)+(12/pow(x[0].x,3))-5; //se tiene que ajustar
	//return 3*pow(x[0].x,4)+pow(x[0].x,2)-2*x[0].x+1;
	//return 10+(pow(x[0].x,3))-(2*x[0].x)+(5*exp(x[0].x));  //no tiene minimo
	//return  pow(x[0].x,2) -(10*exp(0.1*x[0].x));
	//return pow((10*pow(x[0].x,3)+3*pow(x[0].x,2)+5),2);
	//return (0.5/sqrt(1+pow(x[0].x,2))) -(sqrt(1+pow(x[0].x,2)))*(1-(0.5/(1+pow(x[0].x,2))))+x[0].x; //no tiene minimo
	//return exp(x[0].x) -(pow(x[0].x,3));
	//return pow((pow(x[0].x,2)-1),3) - (pow((2*x[0].x-5),4));
	//return (-4*pow(x[0].x,2)-20*x[0].x-100)+pow((1-x[0].x),4);
	//return (pow(x[0].x,2)+pow((x[1].x+1),2))*(pow(x[0].x,2)+pow((x[1].x-1),2));
	//return pow(pow(x[0].x,2)-x[1].x,2)+pow(x[1].x,2);
	//return 50*(x[1].x-pow(x[0].x,2),2)+pow((2-x[0].x),2);
	//return 50*pow((x[1].x-pow(x[0].x,2)),2)+pow((2-x[0].x),2);
	//return pow((x[0].x + 2*x[1].x -7),2) + pow((2*x[0].x+x[1].x -5),2);
	//return pow((1.5-x[0].x*(1-x[1].x)),2)+pow((2.25-x[0].x*(1-pow(x[1].x,2))),2)+ pow((2.625-x[0].x*(1-pow(x[1].x,3))),2);
	//return pow((10*(x[1].x-pow(x[0].x,2))),2) +pow((1-x[0].x),2)+ 90*pow((x[3].x-pow(x[2].x,2)),2)+pow((1-x[2].x),2)+10*pow((x[1].x+x[3].x-2),2)+0.1*(x[1].x-x[3].x);
	//return (4-2.1*pow(x[0].x,2)+(pow(x[0].x,4)/3))*pow(x[0].x,2)+(x[0].x*x[1].x)+(-4+(4*pow(x[1].x,2)))*pow(x[1].x,2);
	//return pow((x[0].x+10*x[1].x),2)+(5*pow((x[2].x-x[3].x),2))+pow((x[1].x-2*x[2].x),4)+(10*pow((x[0].x-x[3].x),4));//no tiene minimo
	/*for(i=0;i<3;i++)
	{
		res+=pow(x[i].x,2);
	}
	return res;*/
	//return 100*(pow((pow(x[0].x,2)-x[1].x),2))+pow((1-x[0].x),2);
	/*for(i=0;i<5;i++)
	{
		res+= floor(x[i].x);
	}
	return res;*/
	/*for(i=0;i<30;i++)
	{
		res+= (i+1)*pow(x[i].x,4);
	}
	return res+genGauss(0,1);*/
	for(j =0;j<25;j++)
	{	
		double var; 
		/*aleatorio(-65.536,65.536)*/
		for(i=0;i<2;i++)
		{
			var +=pow((x[i].x-a[i][j]),6);
		}
		res+=1/((j+1)+var);
	}
	
	return 1/((1/k)+res);
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


void aleatorios(int *num_ale, int mu)
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
			for (y=0;y<tam_v;y++)
			{	
				printf("%d   %f ",(x+1),matr_Pad[x][y].x);
			}
			printf("    f(x)=%f\n",evaluciones[x]);
		}


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

void ordenar(x_delt **uv, double *evalTot,int tamT,int num_var)
{
	int x;
	double *vecEval;	
	
	for(x=0;x<tamT;x++)
	{
		int y,min,z;
		double aux;
		min =x;
		for(y=(x+1);y<tamT;y++)
		{
			if (evalTot[y]<evalTot[min])
				min=y;
		}
		aux = evalTot[x];
		evalTot[x] =evalTot[min];
		evalTot[min] =aux;
			
		for(z=0;z<num_var;z++)
		{
			double v_a=uv[x][z].x;
			double del=uv[x][z].sigma;
				
			uv[x][z].x=uv[min][z].x;
			uv[min][z].x=v_a;
				
			uv[x][z].sigma=uv[min][z].sigma;
			uv[min][z].sigma=del;
			
		}
	}


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

int main()
{
	double **sigma=NULL, *padre=NULL, padreFunc, hijoFunc,num,aprox,*evaluciones=NULL,*eval_H=NULL,*evalTot=NULL;
	int t, num_var,i,x,ps,mu,landa,mejorS,tamT;
	t =0;
	ps =0;
	num_var =2;	
	srand(time(NULL));
	aprox =1;
	mu= 20;
	landa =30;
	tamT = mu+ landa;
	//sigma = (double **) malloc (landa*sizeof(double*)); //explorar mas el espacio de busqueda reproduccion
	//hijos=(double **) malloc (landa*sizeof(double*));
	//formando el conjunto de evaluciones
	evaluciones= vector(evaluciones,mu);
	eval_H = vector(eval_H,landa);
	evalTot = vector(evalTot,tamT);
	//matr_Pad=(double **) malloc (mu*sizeof(double*));
	
	x_delt **matr_Pad, **hijos, **uv;
	
	matr_Pad = (x_delt **) malloc (mu*sizeof(x_delt*));
	hijos = (x_delt **) malloc (landa*sizeof(x_delt*));
	uv = (x_delt **) malloc ((landa+mu)*sizeof(x_delt*));
	
	//inicilaizar padre
	for (x=0;x<mu;x++)
	{
		int y;
		double datosP;
		matr_Pad[x] = (x_delt*)malloc(num_var*sizeof(x_delt));
		for(y=0;y<num_var;y++)
		{	
			
			datosP=genGauss(0.0,1.0);
			matr_Pad[x][y].x = datosP; //se evalua la funcion
			matr_Pad[x][y].sigma = 0.0873;
		}
		
		evaluciones[x]=evaluarF(matr_Pad[x]);
			
	}
	
	
	mejorS =buscarMpad(evaluciones,mu);
	//printf("\nMejor resultado ---   %f    f(x)=%f\n",matr_Pad[mejorS][0].x,evaluciones[mejorS]);
	printf("\nMejor resultado ---   %f   %f    f(x)=%f\n",matr_Pad[mejorS][0].x,matr_Pad[mejorS][1].x,evaluciones[mejorS]);
	
	
	while((t<Gmax) && (aprox>error))
	{
		int y,num_ale[2],s,m,nu;
		double M_hijo=0,var=0.0,tau,tau_p;
		printf("\n--- Iteracion Numero %d ---",t+1);
		
		//print(matr_Pad,evaluciones,mu,num_var);
		printf("\n");
		nu =genGauss(0.0,1.0);
		
		tau=(1/(4.0*pow(num_var,0.25)));
		tau_p=(1/(pow(2.0*(num_var),0.5)));
		
		
		for(s=0;s<landa;s++)
		{	
			/*se esscogen 2 sol aleatorias que no se repiten, LA MU SE PUEDE PONER EN 10 y trabaja mejor*/
            aleatorios(num_ale,mu-10);
			//printf("--%d --%d\n",num_ale[0],num_ale[1]);
			//hijos[s] = (x_delt*)malloc(num_var*sizeof(x_delt));
			hijos[s] = (x_delt*)malloc(num_var*sizeof(x_delt));
			for(i=0;i<num_var;i++)
			{		
				double h,d;
				num=genGauss(0.0,1.0);//semilla de aleatorios
				//aplicando recombinacion plana
				h = ((num*matr_Pad[num_ale[0]][i].x)+((1-num)*matr_Pad[num_ale[1]][i].x));
				d =	(matr_Pad[num_ale[0]][i].sigma+matr_Pad[num_ale[1]][i].sigma)/2;
				//h = (matr_Pad[num_ale[0]][i]+matr_Pad[num_ale[1]][i])/2; 			
				//printf("//%f",h);
				//aplicando mutacion
				hijos[s][i].sigma=d*(exp((tau*num)+(tau_p*nu)));
				//printf("********** %f\n",sigma[s][i]);
				h=h+hijos[s][i].sigma*num;
				//printf(" -- %f  \n",h);
				hijos[s][i].x=h;
			}
			
			eval_H[s]=evaluarF(hijos[s]);
			//printf("%f\n",evaluarF(hijos[s]));
		}

		//print(hijos,eval_H,landa,num_var);
		m=mejorS;
		mejorS =buscarMpad(eval_H,landa);
		//printf("\nMejor resultado ---   %f    f(x)=%f\n",hijos[mejorS][0].x,eval_H[mejorS]);
		//printf("\nMejor resultado ---   %f   %f    f(x)=%f\n",hijos[mejorS][0].x,hijos[mejorS][1].x,evaluciones[mejorS]);
		//printf("\nMejor resultado ---   %f   %f  %f  f(x)=%f\n",hijos[mejorS][0].x,hijos[mejorS][1].x,hijos[mejorS][2].x,evaluciones[mejorS]);
		//printf("\nMejor resultado ---   %f   %f  %f %f  f(x)=%f\n",hijos[mejorS][0].x,hijos[mejorS][1].x,hijos[mejorS][2].x,hijos[mejorS][3].x,evaluciones[mejorS]);
		//--printf("\nMejor resultado ---   %f   %f  %f  %f  %f  f(x)=%f\n",hijos[mejorS][0].x,hijos[mejorS][1].x,hijos[mejorS][2].x,hijos[mejorS][3].x,hijos[mejorS][4].x,evaluciones[mejorS]);
		printMS(hijos,evaluciones,mejorS,num_var);
		
		aprox=fabs(evaluciones[m]-eval_H[mejorS]);
		printf("Error: %f \n",aprox);
		
		Fuv(matr_Pad, hijos, uv,evalTot, evaluciones, eval_H, mu, landa, num_var);
		
		//ordenamiento por seleccion, se puede mejorar
		ordenar(uv, evalTot, tamT, num_var);
		
		//print(uv,evalTot,tamT,num_var);
		
		//sustituyecdo a los padres
		for (x=0;x<mu;x++)
		{
			int y;
			for(y=0;y<num_var;y++)
			{					
				matr_Pad[x][y].x = uv[x][y].x;
				matr_Pad[x][y].sigma = uv[x][y].sigma;
			}
			evaluciones[x] = evalTot[x];
		}
		
		t++;
	}


    return 0;
}