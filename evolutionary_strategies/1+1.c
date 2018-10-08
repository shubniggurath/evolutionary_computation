#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define  TWO_PI 6.2831853071795864769252866
#define c 0.817
#define Gmax 200000
#define error 0.000001
#define ranf() ((double)rand()/(1.0+(double)RAND_MAX))  //intervalo uniforme de [0,1)


double aleatorio(double,double);
double  *vector(double *,int);
double inic_padre(double *,int);
double evaluarF(double*);
double generateGaussianNoise(const double);
double ale(double, double);


double ale(double min, double max)
{	
	double scaled = ranf();
       return (max - min)*scaled + min;
}

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

double aleatorio(double media, double desE)
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
	//return pow((10*pow(x[0],3)+3*2*x[0]+5),2);
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
	return res+aleatorio(0,1);*/
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
	
}

double inic_padre(double * padre,int num_var)
{
	int i;
	for(i=0;i<num_var;i++)
	{	
		padre[i]=aleatorio(0,1.0);
		//padre[i]=generateGaussianNoise(1);
	}
	
	return evaluarF(padre);
}



int main()
{
	double delta, *padre=NULL, *hijo=NULL, padreFunc, hijoFunc,num,aprox;
	int t, num_var,padre_sus,i,ps, ve=0;
	t =0;
	ps =0;
	num_var =2;
	delta =30; //explorar mas el espacio de busqueda reproduccion
	srand(time(0));
	padre =vector(padre,num_var);
	hijo=vector(hijo,num_var);
	aprox =1;
	padreFunc = inic_padre(padre,num_var); //se evalua la funcion
	
    while((t<Gmax) && (aprox>error))
	{
		printf("--- Iteracion Numero %d ---\n",t+1);
		for(i=0;i<num_var;i++){
             num=aleatorio(0,1.0);//semilla de aleatorios
			 //num= generateGaussianNoise(1);
             hijo[i]=padre[i]+delta*num; // mutando el vector
             printf("    Delta: %f  Aleatorio: %f \n",delta,num);
             printf("    Padre %d : %f \n",i+1,padre[i]);
             printf("    Hijo %d : %f \n",i+1,hijo[i]);
        }
		
		hijoFunc = evaluarF(hijo); //evalua la funcion con los nuevos datos
		printf("    Evaluacion padre: %f \n",padreFunc);
        printf("    Evaluacion hijo: %f \n",hijoFunc);
		
		aprox=fabs(padreFunc-hijoFunc);
        printf("    Error: %f \n",aprox);
		//printf("%f,%f",hijoFunc,padreFunc);
        if(hijoFunc<padreFunc){
            printf("    HIJO SUSTITUYO AL PADRE");
            for(i=0;i<num_var;i++){
                padre[i]=hijo[i];
            }
            padreFunc=hijoFunc;
            ps++;
			ve++;
        }
		
		//usando el 1/5
		/*if((t%(10*num_var)) == 0){
            if((ps/(10*num_var))<(1/5)){
				delta=delta*c;
            }else if((ps/(10*num_var))>(1/5)){
                delta= (delta/c);
            }else{
				delta=delta;
            }
            ps=0;
        }*/
		
		if((t%(10*num_var)) == 0){
            if(ps<2){
                delta= delta*c;
            }else if(ps>2){
                delta=delta/c;
            }else{
				delta=delta;              
            }
            ps=0;
        }
		
		
		t++;
	}
	printf("\n---%d\n",ve);
    return 0;
}
