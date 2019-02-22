#include <stdio.h>
double trapzdnew(double *g, int n1,int n2, int n)
{
		int i;
		double h;
		double sum;
		sum=0.0;
		h=1.0/(n-1); 
		for (i=n1;i<=n2;i++)
			sum=sum+g[i]*h;
		sum=sum-0.5*h*(g[n1]+g[n2]);
		
		return sum;
}
