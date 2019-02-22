#include <math.h>
#include <stdio.h>
#define Pi 3.141592653589793
#define S 0.00000000000001
double Stat2new (double **u,int n, int a,int b, double R0, double K, double H)
{
	double zero,v1,v2,v4,v5,v6,v7,s1,s2,h2i,h,v8,v9;
	int i,j,k1,k2;
	zero=0.0;h=1.0/(n-1); 
	h2i=1.0/(h*h);
	for(i=2;i<(int)(n*R0-R0+1)-b;i++)
	{
		s2=0.0;
		for (j=b+2;j<(int)(n*H-H+1);j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			         
								 
				/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}
		for(i=2;i<a-b;i++)
	{
		s2=0.0;
		for (j=2;j<=1+b;j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			          
								 
				/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}
	for(i=a+b+1;i<(int)(n*R0-R0+1);i++)
	{
		s2=0.0;
		for (j=2;j<=1+b;j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			           
								 
				/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}
		for(i=(int)(n*R0-R0+1)-b;i<(int)(n*R0-R0+1);i++)
	{
		s2=0.0;
		for (j=b+2;j<(int)(n*H-H+1)-b;j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];            
								 
			/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}
	
	zero=zero/(((int)(n*R0-R0+1)-2)*((int)(n*H-H+1)-2)-(b-1)*(b-1)-(2*b-1)*(b-1));
	return zero;
}
