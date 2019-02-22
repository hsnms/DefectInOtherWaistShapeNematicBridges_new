#include <math.h>
#include <stdio.h>
#include "Header2new.h"
/*#define EPI 1.0e-7*/ 
#define Pi 3.141592653589793
void **nrfunc2new( double **f, double **u, int n, double *q, int a)
{
	double h,v1,v2,v4,v5,v6,v7,v8,v9,K;
	int i,j,m;
	h=1.0/(n-1); 	
	K=1.0;
	m=65;
	for(j=1;j<=n;j++)
		for(i=1;i<=n;i++)
			f[i][j]=0.0;


	for (j=2;j<=n-1;j++)
	{
		for (i=2;i<=n-1;i++)
		{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	        
			f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);
				/*K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);*/
		}
	}
		
	for (i=2;i<=n-1;i++)
	{
			v1=u[i][2]-u[i][1];
		    v2=u[i][2]-u[i][1];
			v4=u[i+1][1]-u[i][1];
			v5=u[i][1]-u[i-1][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];	

		f[i][1]=pow(sin(u[i][1]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][1])*(v4+v5)/(2.0*h)-pow(sin(u[i][1]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);
			/*K*pow(sin(u[i][1]),2)/(h*(i-1))+((K*pow(cos(u[i][1]),2)+1.0*pow(sin(u[i][1]),2))*(i-1)+(K*pow(sin(u[i][1]),2)+1.0*pow(cos(u[i][1]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][1])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][1]),2)+1.0*pow(cos(u[i][1]),2))*(i-1)*pow(q[1]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][1]),2)+1.0*pow(cos(u[i][1]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[1]+(K-1.0)*sin(2.0*u[i][1])*(i-1)*q[1])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][1])-K*(cos(2.0*u[i][1])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][1])-1.0)*q[1]*(v1+v2)/(2.0*h);*/

	}
	for (i=2;i<=n-1;i++)
	{
			v1=u[i][n]-u[i][n-1];
		    v2=u[i][n]-u[i][n-1];
			v4=u[i+1][n]-u[i][n];
			v5=u[i][n]-u[i-1][n];
			v8=q[n]-q[n-1];
			v9=q[n]-q[n-1];	           
		f[i][n]=pow(sin(u[i][n]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][n])*(v4+v5)/(2.0*h)-pow(sin(u[i][n]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);
			/*K*pow(sin(u[i][n]),2)/(h*(i-1))+((K*pow(cos(u[i][n]),2)+1.0*pow(sin(u[i][n]),2))*(i-1)+(K*pow(sin(u[i][n]),2)+1.0*pow(cos(u[i][n]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][n])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][n]),2)+1.0*pow(cos(u[i][n]),2))*(i-1)*pow(q[n]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][n]),2)+1.0*pow(cos(u[i][n]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[n]+(K-1.0)*sin(2.0*u[i][n])*(i-1)*q[n])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][n])-K*(cos(2.0*u[i][n])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][n])-1.0)*q[n]*(v1+v2)/(2.0*h);*/

	}
	for (j=2;j<=n-1;j++)
	{
		    v1=u[1][j+1]-u[1][j];
		    v2=u[1][j]-u[1][j-1];
			v4=u[2][j]-u[1][j];
			v5=u[2][j]-u[1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		       
			i=1;
		f[1][j]=sin(2.0*u[1][j])*(v4+v5)/(2.0*h)-pow(sin(u[1][j]),2)*q[j]*(v1+v2)/(h);
			/*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);*/
	}
	for (j=2;j<=n-1;j++)
	{
		    v1=u[n][j+1]-u[n][j];
		    v2=u[n][j]-u[n][j-1];
			v4=u[n][j]-u[n-1][j];
			v5=u[n][j]-u[n-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		        
			i=n;
		f[n][j]=pow(sin(u[n][j]),2)/(1.0*(n-1)*h)+sin(2.0*u[n][j])*(v4+v5)/(2.0*h)-pow(sin(u[n][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(n-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(n-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);
			/*K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);*/
	}
	f[1][1]=sin(2.0*u[1][1])*(u[2][1]-u[1][1])/(1.0*h)-2.0*pow(sin(u[1][1]),2)*q[j]*(u[1][2]-u[1][1])/(1.0*h);
	f[1][n]=sin(2.0*u[1][n])*(u[2][n]-u[1][n])/(1.0*h)-2.0*pow(sin(u[1][n]),2)*q[j]*(u[1][n]-u[1][n-1])/(1.0*h);
f[n][1]=pow(sin(u[n][1]),2)/(1.0*(n-1)*h)+sin(2.0*u[n][1])*(u[n][1]-u[n-1][1])/(1.0*h)-2.0*pow(sin(u[n][1]),2)*q[j]*(u[n][2]-u[n][1])/(1.0*h)+(1.0*(n-1)*h)*pow((u[n][1]-u[n-1][1])/(1.0*h),2)+(1.0*(n-1)*h)*pow(q[j],2)*pow((u[n][2]-u[n][1])/(1.0*h),2);
	f[n][n]=pow(sin(u[n][n]),2)/(1.0*(n-1)*h)+sin(2.0*u[n][n])*(u[n][n]-u[n-1][n])/(1.0*h)-2.0*pow(sin(u[n][n]),2)*q[j]*(u[n][n]-u[n][n-1])/(1.0*h)+(1.0*(n-1)*h)*pow((u[n][n]-u[n-1][n])/(1.0*h),2)+(1.0*(n-1)*h)*pow(q[j],2)*pow((u[n][n]-u[n][n-1])/(1.0*h),2);
 

	for (i=n-(int)((1/q[1])*((n-1)/(m-1)));i<=n;i++)
		for (j=n-(n-1)/(m-1);j<=n;j++)
			f[i][j]=0.0;
	if (a>=1&&a<=(int)((1/q[1])*((n-1)/(m-1))))
	{
		for (i=1;i<=1+(int)((1/q[1])*((n-1)/(m-1)));i++)
			for (j=1;j<=1+(n-1)/(m-1);j++)
				f[i][j]=0.0;
				
	}
	else if (a>=1+(int)((1/q[1])*((n-1)/(m-1)))&&a<=n-(int)((1/q[1])*((n-1)/(m-1))))
		{
		for (i=a-(int)((1/q[1])*((n-1)/(m-1)));i<=a+(int)((1/q[1])*((n-1)/(m-1)));i++)
			for (j=1;j<=1+(n-1)/(m-1);j++)
				f[i][j]=0.0;
				
	}
	else
	{
		for (i=n-2*(int)((1/q[1])*((n-1)/(m-1)));i<=n;i++)
			for (j=1;j<=1+(n-1)/(m-1);j++)
				f[i][j]=0.0;

	}
	return 0;
}
