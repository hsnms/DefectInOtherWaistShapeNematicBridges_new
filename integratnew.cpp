#include "Header2new.h"
#include "Header5new.h"
double integratnew (double **f, int n1, int n2, int n3, int n4, int n, double H, double R0)
{
	int i;
        double *k;
	double  integ;
	k=dvectornew(1,(int)(n*R0-R0+1));
	for (i=n1;i<=n2;i++)
	{
		k[i]=trapzdnew(f[i],n3,n4,n);
	}
	integ=trapzdnew(k,n1,n2,n);
	free_dvectornew(k,1,(int)(n*R0-R0+1));
	return integ;
}
