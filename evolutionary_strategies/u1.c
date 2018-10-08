#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define M_PI 3.14159265358979323846
#define  TWO_PI 6.2831853071795864769252866
#define c 0.817
#define Gmax 200000
#define error 0.0000001  //le ponemos 4 ceros y genera mas rapido
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)
#define codo 0.52355866484*10
#define aureo (1+sqrt(5))/2


double aleatorio(double,double);
double  *vector(double *,int);
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
	int i,j,k=500;
	double res=0.0;
	int a[2][25]={{-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,16,32},{-32,-32,-32,-32,-32,-16,-16,-16,-16,-16,0,0,0,0,0,16,16,16,16,16,32,32,32,32,32}};
	//return (-0.75/(1+pow(x[0],2))) - (0.65*x[0]*atan(1/x[0]))+0.65;
	//return ((-4*pow(x[0],2)) -20*x[0] -100) + pow(1-x[0],4);
	//return 3*pow(x[0],4)+pow(x[0],2)-2*x[0]+1;
	//return 10+(pow(x[0],3))-(2*x[0])+(5*exp(x[0]));  //no tiene minimo
	//return  pow(x[0],2) -(10*exp(0.1*x[0]));
	//return pow((10*pow(x[0],3)+3*pow(x[0],2)+5),2);
	//return (0.5/sqrt(1+pow(x[0],2))) -(sqrt(1+pow(x[0],2)))*(1-(0.5/(1+pow(x[0],2))))+x[0];//no tiene minimo
	//return exp(x[0]) -(pow(x[0],3));
	//return pow((pow(x[0],2)-1),3) - (pow((2*x[0]-5),4));
	//return (-4*pow(x[0],2)-20*x[0]-100)+pow((1-x[0]),4);
	//return (pow(x[0],2)+pow((x[1]+1),2))*(pow(x[0],2)+pow((x[1]-1),2));
	//return pow(pow(x[0],2)-x[1],2)+pow(x[1],2);
	//return 50*pow((x[1]-pow(x[0],2)),2)+pow((2-x[0]),2);
	//return pow((x[0] + 2*x[1] -7),2) + pow((2*x[0]+x[1] -5),2);
	//return pow((1.5-x[0]*(1-x[1])),2)+pow((2.25-x[0]*(1-pow(x[1],2))),2)+ pow((2.625-x[0]*(1-pow(x[1],3))),2);
	//return pow((10*(x[1]-pow(x[0],2))),2) +pow((1-x[0]),2)+ 90*pow((x[3]-pow(x[2],2)),2)+pow((1-x[2]),2)+10*pow((x[1]+x[3]-2),2)+0.1*(x[1]-x[3]);
	//return (4-2.1*pow(x[0],2)+(pow(x[0],4)/3))*pow(x[0],2)+(x[0]*x[1])+(-4+(4*pow(x[1],2)))*pow(x[1],2);
	//return pow((x[0]+10*x[1]),2)+(5*pow((x[2]-x[3]),2))+pow((x[1]-2*x[2]),4)+(10*pow((x[0]-x[3]),4));//no tiene minimo
	/*for(i=0;i<3;i++)
	{
		res+=pow(x[i],2);
	}
	return res;*/
	//return 100*(pow((pow(x[0],2)-x[1]),2))+pow((1-x[0]),2);
	/*for(i=0;i<5;i++)
	{
		res+= floor(x[i]);
	}
	return res;*/
	
	/*for(i=0;i<30;i++)
	{
		res+= (i+1)*pow(x[i],4);
	}
	return res+genGauss(0,1);*/
	for(j =0;j<25;j++)
	{	
		double var; 
		for(i=0;i<2;i++)
		{
			var +=pow((x[i]-a[i][j]),6);
		}
		res+=1/((j+1)+var);
	}
	
	return 1/((1/k)+res);
	//return 3*pow(x[0],4)+pow(x[0],2)-2*x[0]+1; 
	//return pow((x[0]+10*x[1]),2)+5*pow((x[2]-x[3]),2)+pow(x[1]-2*x[2],4)+10*pow((x[0]-x[3]),4);
}



int main()
{
	double delta, *padre=NULL, *hijo=NULL, padreFunc, hijoFunc,num,aprox,*padres=NULL,**matr_Pad=NULL;
	int t, num_var,i,x,ps, ve=0,mu;
	t =0;
	ps =0;
	num_var =2;
	delta = codo; //explorar mas el espacio de busqueda reproduccion
	srand(time(NULL));
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
			
			datosP=genGauss(0.0,1.0);
			matr_Pad[x][y] = datosP; //se evalua la funcion
		}
		
		padres[x]=evaluarF(matr_Pad[x]);
			
	}
	
	
    while((t<Gmax) && (aprox>error))
	{
		printf("--- Iteracion Numero %d ---\n",t+1);

		//evaluar soluciones
		double var=0.0;
		int y, num_ale=0,num_ale2=0, x=0, pos=0,s,opt=0,k;
		for (x=0;x<mu;x++)
		{	
			for (y=0;y<num_var;y++)
			{	
				printf("%d   %f ",(x+1),matr_Pad[x][y]);
			}
			printf("    f(x)=%f\n",padres[x]);
		}
		
		/*se esscogen 2 sol aleatorias que no se repiten*/	
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
			
            num=genGauss(0.0,1.0);//semilla de aleatorios
			//2 + η(6-2),4 + η(1-4)
			hijo[i]= ((matr_Pad[num_ale][i]+matr_Pad[num_ale2][i])/2)+delta*num; // mutando el vector
            /* printf("    Delta: %f  Aleatorio: %f \n",delta,num);
             printf("    Padre %d : %f \n",i+1,padre[i]);*/
             printf("Hijo %d : %f ",i+1,hijo[i]);
        }
		hijoFunc=evaluarF(hijo);
		
		printf("%f",hijoFunc);
		//busca al peor de los padres
		var=padres[0];
		for(s=0;s<mu;s++)
		{
			if (var<padres[s])
			{
				var=padres[s];
				pos=s;
			}

        }
		printf("\npeor solucion---- %d ",(pos+1));
		
		//printf("\nmejor solucion---- %d  \n",(opt+1));
		//buscando al mejor de los padres
		var=padres[0];
		for(k=0;k<mu;k++)
		{
			if (var>padres[k])
			{
				var=padres[k];
				opt=k;
			}
        }
		printf("\nmejor solucion---- %d  \n",(opt+1));
		
		if(hijoFunc<padres[opt])
		{
			ps++;
			ve++;
		}
		
		
		//cambiando por la peor solucion
		for(i=0;i<num_var;i++)
		{
			matr_Pad[pos][i]=hijo[i];
        }
		aprox=fabs(padres[opt]-hijoFunc);
		printf("Error: %f \n",aprox);
		padres[pos]=hijoFunc;
		
		//usando el 1/5
		if((t%(5*num_var)) == 0){
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
