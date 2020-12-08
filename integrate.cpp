//Numerical recipes en c++

#include<stdio.h>
#include <math.h>
#include <stddef.h>
#include<stdlib.h>
#define EPS 3.0e-11
#define NPOINT 10/*numero de polinomio de legendre*/
#define X1 0//limite inferior
#define X2 10//limite superior
#define X3 10.0
#define NR_END 1
#define H 70
#define FREE_ARG char*



void nrerror(char error_text[])
{
  fprintf (stderr,"error en la rutina\n");
  fprintf (stderr,"%s\n",error_text);
  fprintf (stderr,"saliendo a sistema\n");
  exit (1);
}


float *vector(long nl, long nh);
void gauleg(float x1, float x2, float x[], float w[], int n);

float *vector(long nl, long nh)
{
  float *v;
  v= (float *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
  if(!v) nrerror ("allocation faillure in vector ()");
  return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void gauleg(float x1, float x2, float x[], float w[], int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1; i<=m; i++)
    {
      z=cos (3.141592654*(i-0.25)/(n+0.5));
      /*printf(" here 1 %12.6f \n",z);*/
      do
	{
	  p1=1.0;
	  p2=0.0;
	  for (j=1; j<=n;j++)
	    {
	      p3=p2;
	      p2=p1;
	      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	    }
	  pp=n*(z*p1-p2)/(z*z-1.0);
	  
	  z1=z;
	  
	  z=z1-p1/pp;
	}
      while (fabs(z-z1)>EPS);
      x[i]=xm-xl*z;
      x[n+1-i]=xm+xl*z;/*parametrización de la funcion*/
      w[i]=2.0*xl/((1.0-z*z)*pp*pp);/*evaluación de la función*/
      w[n+1-i]=w[i];
  
    }
}



/*drive for gauleg*/

float func(float x)
{
  
  float Ol=0.0;
  float Or=0.00000;
  float Ok=0;
  float Om=0.3;
  float Omega = Om + Or + Ol;
//return 1/( (1+x)* sqrt((Om*(1+x)/ + (Ol*pow((1+x),2)) ) )) ;//distancia diametral
 //return sqrt(pow((1+x)*Om,3)+pow((1+x),4)*Or+pow((1+x),2)*Ok+1+Ol);
 return 1/sqrt( (Or/(pow(x,2))) + Om/x + Ol*(pow(x,2)) + (1- Omega) );

}

int main(void)
{
  int i,j;
  float xx=0.0;
  float *x,*w,xf,xin,dx,x0;
  float a=X1;
  float b=X2;
  int NI=1000;/*numero de particiones en el intervalo*/
  float integ[NI];
  float integral;
  float suma=0.0;
  
  
  x=vector(1,NPOINT);
  w=vector(1,NPOINT);
  
  
  x0=a;
  dx=(b-a)/(float)NI;

  for(j=1;j<=NI;j++)
    {
      ;
      xf=x0+j*dx;
      gauleg(x0,xf,x,w,NPOINT);
      suma=0.0;
      for(i=1;i<=NPOINT;i++)  
	{ 
	  suma+=(w[i]*func(x[i]));/*integrador*/
	}
      integ[j]=suma;
    }
  
  FILE *Integral_gauss; /*crea el archivo de datos*/
  Integral_gauss= fopen("modelo5.dat","w");

  for(j=1;j<=NI;j++)
    {
    /*printf("L %12.6f %12.6f\n",xf,suma);*/
    fprintf(Integral_gauss,"%f %f \n",a+j*dx,integ[j]);
    }
    fclose(Integral_gauss);

  /*demostrate the use of gauleg for integration*/
  gauleg(X1,X2,x,w,NPOINT);
  for(i=1;i<=NPOINT;i++)
    xx += (w[i]*func(x[i]));
  printf("\n Integral: %12.6f \n",xx);

  free_vector(w,1,NPOINT); 
  free_vector(x,1,NPOINT);   
}
