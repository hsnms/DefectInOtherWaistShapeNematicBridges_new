#include <stdio.h>
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Left 1.570796326794897
#define Right 1.570796326794897
#define Pi 3.141592653589793
void fill02new(double **u, int n, double m, int a) /*problem*/
	 
{ 
	int i,j; 
  
	/*for(i=1;i<=3;i++)*/
		u[2][1]=((Down1-0.5*m*Pi)*a/n)+(Down2*(n-a)/n);
	/*for (i=1;i<=3;i++)*/
	
		u[2][3]=Up-0.5*m*Pi;
	
	for (j=1;j<=3;j++)
	{
		u[1][j]=Left-0.5*m*Pi;
	}
	u[2][2]=Up-0.5*m*Pi;
	u[2][3]=Up-0.5*m*Pi;
	u[3][3]=Right;
	u[3][2]=Right;
	u[3][1]=Right;

}
