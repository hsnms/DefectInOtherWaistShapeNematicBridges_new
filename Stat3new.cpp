#include <math.h>
#include <stdio.h>
#define Pi 3.141592653589793
#define S 0.00000000000001
double Stat3new (double **u,int n,int b, double R0, double K, double H, int *qq)
{
	double zero,v1,v2,v4,v5,v6,v7,s1,s2,h2i,h,v8,v9;
	int i,j,k1,k2,mm;
	zero=0.0;h=1.0/(n-1); 
	h2i=1.0/(h*h);
	mm=0;
	for(j=b+2;j<(int)(n*H-H+1)/*-b*/;j++)
	{
		s2=0.0;
		for (i=2;i<qq[j-1];i++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];

                    	/*if (i==qq[j-1])
			 {
			   v7=2.0*(u[i-1][j-1]-u[i][j-1]);
		         }

	                if (i>=qq[j-1]+1)
		         {
		         v2=v1;
		        v6=2.0*(u[i+1][j+1]-u[i-1][j+1]);
		         v7=2.0*(u[i-1][j]-u[i+1][j]);
			 }*/
						
			
			    
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			mm=mm+1;
			
		}
		zero=zero+s2;
	}
		for (j=2;j<=1+b;j++)
	{
		s2=0.0;
		for (i=b+2;i<qq[j-1];i++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			         
			/*	if (i==qq[j-1])
						  {
						    v7=2.0*(u[i-1][j-1]-u[i][j-1]);
						  }

						if (i>=qq[j-1]+1)
						  {
						    v2=v1;
						    v6=2.0*(u[i+1][j+1]-u[i-1][j+1]);
						    v7=2.0*(u[i-1][j]-u[i+1][j]);
						    }*/
											 
				
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			mm=mm+1;
			
		}
		zero=zero+s2;
	}
		//not calculate the lattices near the upper and curved boundary
		/*	for(j=(int)(n*H-H+1)-b;j<(int)(n*H-H+1);j++)
	{
		s2=0.0;
		for (i=2;i<qq[(int)(n*H-H+1)-b]-b;i++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
		          
								 
				
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2))*0.25*pow(v4+v5,2)+(K-1.0)*(
				cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)))*(v4+v5)+(K-1.0)*sin(2.0*u[i][j])*(i-1)*0.5*(v1+v2)
-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2))*(v4-v5)-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v1-v2)
+((K-1.0)*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			mm=mm+1;
			
		}
		zero=zero+s2;
		}*/



	
	zero=zero/mm;
	return zero;
}
