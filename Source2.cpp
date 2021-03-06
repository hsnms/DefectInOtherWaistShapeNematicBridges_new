//Source2.cpp helps to find the specific equilibrium nematic defect structures in a waist-shaped nematic bridge given specific parameters including the Frank constant ratio and aspect ratio.

//We consider homeotropic boundary conditions and proper inner boundary conditions (after we apply a cut-off of the defect core), and solve the Euler-Lagrange equation. Because of the cylindrical symmetry, we only consider the rectangular diametrical plane.

//The finite difference method (particularly the over-relaxation method) is used. The rectangular region is not rescaled.

//Input: n: number of lattices horizontally and vertically; R0: aspect ratio; K:Frank constant ratio; b: the number of the lattice spaces the length of the defect core occupies; m:types of defect structures (radial:1.0, hyperbolic:-1.0).

//Output: 1,"Bridge.txt": different energies when the defect core is at different locations, error, all the input values; 2, "BridgeVectorField.txt": the vector field.

//How to compile: 1, make sure the in the "makefile", we have the relevant "Source2"; 2, M-X compile (option + X + compile); 3, Compile command: make -k output; 4, Save file, y; 5, M-X shell; 6, bash-3.2$ ./output

//K=0.3 for radial, any value for hyperbolic

#include <stdio.h>
#include <errno.h>
#include <math.h>
//When the angle is 0, the director points up; when it's Pi, it points down
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Right 1.570796326794897
#define Left 1.570796326794897
#define Pi 3.141592653589793
#define S 0.00000000000001
#include "Header1new.h"
#include "Header2new.h"
#include "Header5new.h"
#include "Header6new.h"
#include "Header7new.h"
int main()
{
    int i,j,k,n,ncycle=5,ii,ii2,ii3,ii4,sos,kk,radius,min,n1,RR;
    double **u,*x,*y,*z,*z2,*q,**v,**f;
    int *qq;
    double result,zero,rev2,rev3,v1,v2,v4,K,t,sum,result22,sum2,fff,dfff,KK,XX,h1,h2;
	double m,point,point2,H,h,R0; int a,b; 
	FILE *fp,*fp2,*fp3;
	h1=0.2;//h2=5/4;
	h2=1.0;
	//	printf("h2=%f\n",h2);
	printf ("Number of lattices vertically, n:");
	scanf("%d", &n);
	//n=33;
	printf ("Aspect ratio, R0:");
	scanf("%lf", &R0);
	//R0=3.0;
	H=1.0;
        printf ("Frank constant ratio, K:");
	scanf("%lf", &K);
	//K=0.74;
	//printf("%d\n",(int)(n*R0-R0+1));
	
	fp = fopen ("Bridge.txt","w");
	  if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	
	fp2 = fopen ("BridgeVectorField.txt","w");
	fp3=fopen ("Boundary.txt","w");
  
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
		if (fp3 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
	}
	  
	   printf ("How many of the lattice spaces does the length of the defect core occupy? b:");
	scanf("%d", &b);
	// b=2;/*radius of the point or small radius of the ring*/
	kk=1+2*b;/*location of the smallest ring*/
        printf ("Types of defect structures (radial:1.0, hyperbolic:-1.0), m:");
	scanf("%lf", &m);
	//m=-1.0;
	a=kk;
	//printf("%d\n",(int)(n*R0-R0+1));

	  

	    
    u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	v=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	qq=vectornew(1,(int)(n*H-H+1));
	x=dvectornew(1,(int)(n*R0-R0+1));
	y=dvectornew(1,(int)(n*H-H+1));
	z=dvectornew(1,(int)(n*R0-R0+1));
	z2=dvectornew(1,(int)(n*H-H+1));
	for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
        
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}
	//	printf("Shape of Boundary\n");
	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			
	 q[j]=R0-h1*pow(H,h2)+h1*pow(1.0*(j-1)/(n-1),h2);/*R0+0.75*pow(H,4)-0.75*pow(1.0*(j-1)/(n-1),4) for barrel; 
	R0-0.75*pow(H,4)+0.75*pow(1.0*(j-1)/(n-1),4) for waist*/

		 
	qq[j]=round(n*q[j]-q[j]+1);//number of lattices in each row inside the bridge
        printf("qq[%d]=%d\n",j,qq[j]);
	if (qq[j]>(int)(n*R0-R0+1))//something new
	  {
	    qq[j]=(int)(n*R0-R0+1);
	  }
		 
	}
 

  
  

 
for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=/*Right-0.5*m*Pi-atan(1.0)*/0.0;//initiallly points up
  }
	
   for(i=1;i<=(int)(n*R0-R0+1);i++)
	{
	  u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;//upper boundary,radial points up, hyperbolic points down
	}
   for(j=2;j<=(int)(n*H-H+1)-1;j++)//inside radial points up, hyperbolic points down
		{
			for(i=2;i<=qq[j]-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
   for (i=1;i<=b;i++)//lower boundary, inside core, radial points up, hyperbolic points down
		u[i][1]=Down1-0.5*m*Pi;
	
   for (i=1+b;i<=qq[1];i++)//lower boundary, outside core, points right
		u[i][1]=Down2;



	 for (j=2;j<=(int)(n*H-H+1)-2;j=j+1)
	//for (j=2;j<=(int)(n*H-H+1)-1;j=j+1) //new
		{
		       if (q[j-1]==q[j])
			   {
			     u[qq[j]][j]=Right;//when director points right on the curved boundary
			   }
			   else
			   {
			     u[qq[j]][j]=Right-0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));//directors on the curved boundary

			    
			     /* u[qq[j]][j]=Right-0.5*m*Pi-atan(2.0/((double)((qq[j+1]-qq[j-1])))) for waist;
  u[qq[j]][j]=Right-0.5*m*Pi-atan(2.0/((double)((-qq[j+1]+qq[j-1])))) for barrel*/
			   }
		       u[1][j]=Left-0.5*m*Pi;//left boundary, radial points up, hyperbolic points down
			}

	 j=(int)(n*H-H+1)-1;//second to the last lattice on the curved boundary
		{
		       if (q[j-1]==q[j])
			   {
			   u[qq[j]][j]=Right;
			   }
			   else
			   {
			     u[qq[j]][j]=Right-0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));
			     /* u[qq[j]][j]=Right-0.5*m*Pi-atan(2.0/((double)((qq[j+1]-qq[j-1])))) for waist;
  u[qq[j]][j]=Right-0.5*m*Pi-atan(2.0/((double)((-qq[j+1]+qq[j-1])))) for barrel*/
			    
			   
			   }
		       u[1][j]=Left-0.5*m*Pi;//left boundary,first lattice

		}
        

	       
	       
		for (j=(int)(n*H-H+1)-1; j>=2;j--)//lattices near the curved boundary
		 {
		   for (i=2;i<=qq[j]-1;i++)
		     {
		       	if (i>=qq[j-1])
			  u[i][j]=/*0.5*(u[i][j+1]+u[i+1][j])*/u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1]+1)/(qq[j]-qq[j-1]+1);
		     }
		 }
        
	     



	  
		//hyperbolic type
	       
	if (m==-1.0)
	{
	  //point defect
	for (j=2;j<=1+b;j++)
	{
	  //	u[1+b][j]=atan((y[j]-y[1])/((x[1+b]-x[1])))+0.5*Pi-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*atan((y[j]-y[1])/((x[1+b]-x[1])))));
	  	  u[1+b][j]=atan((y[j]-y[1])/((x[1+b]-x[1])))+0.5*Pi-((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*atan((y[j]-y[1])/((x[1+b]-x[1])))));
		
	}
	for (i=2;i<=1+b;i++)
	{
	  //	u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[1])))+0.5*Pi-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*atan((y[1+b]-y[1])/((x[i]-x[1])))));
	   u[i][1+b]=atan((y[1+b]-y[1])/((x[i]-x[1])))+0.5*Pi-((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*atan((y[1+b]-y[1])/((x[i]-x[1])))));
	 }
	

	
	mgfas5new(u,n,ncycle,m,b,qq,x,y, R0,K, H);//solve PDE, not calculate the lattices near the upper and curved boundary
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc4new(f,u,n,qq,b,R0,K,H);//integrand, the inner region and lower boundary
	
n1=1;
 z[1]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,n1+b,n1,b+1,n, H, R0)-integratnew(f, qq[(int)(n*H-H+1)-b]-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);//total energy
	zero=Stat3new(u,n,b,R0,K, H,qq);//error, not calculate the lattices near the upper and curved boundary.
	result=z[1];
	for (j=1;j<=(int)(n*H-H+1);j++)
		for (i=1;i<=(int)(n*R0-R0+1);i++)
			v[i][j]=u[i][j];
	radius=1;
	
		printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
		//fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);
		fprintf(fp,"{0, %f},\t",z[1]);

	
		//Directors on the boundary
		for (j=2;j<=(int)(n*H-H+1)-2;j=j+1)
			{
		      
			    
			  //	 printf("Arrow[{{%f,%f},{%f,%f}}],\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1),1.0*(qq[j]-2)/(n-1)+0.07*sin(u[qq[j]-1][j]),1.0*(j-1)/(n-1)+0.07*cos(u[qq[j]-1][j]));
			  //	 fprintf(fp,"Arrow[{{%f,%f},{%f,%f}}],\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1),1.0*(qq[j]-2)/(n-1)+0.07*sin(u[qq[j]-1][j]),1.0*(j-1)/(n-1)+0.07*cos(u[qq[j]-1][j]));
		}
    
		//print the curved boundary
		//for (j=2;j<=(int)(n*H-H+1)-2;j=j+1)
		//	{
		// printf("{%f,%f},\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1));
		//  fprintf(fp,"{%f,%f},\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1));
		//	}


      //Ring defect with smallest radius
		a=1+2*b;kk=a;
      for (j=1;j<=(int)(n*H-H+1);j++)//boundary function
		{
			
		  q[j]=R0-h1*pow(H,h2)+h1*pow(1.0*(j-1)/(n-1),h2);
		  qq[j]=round(n*q[j]-q[j]+1);
if (qq[j]>(int)(n*R0-R0+1))//something new
	  {
	    qq[j]=(int)(n*R0-R0+1);
	  }
		 
		}
      //initial and outer boundary conditions
	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }
	
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
        for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=qq[j]-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

        for (j=2;j<=(int)(n*H-H+1)-1;j=j+1)
		{
		       if (q[j-1]==q[j])
			   {
			   u[qq[j]][j]=Right;
			   }
			   else
			   {  u[qq[j]][j]=Right-0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));
			    
			     //u[qq[j]][j]=Right-0.5*m*Pi-atan(2.0/((double)((qq[j+1]-qq[j-1]))));
			 
		}

                           u[1][j]=Left-0.5*m*Pi;
			}
	
	  for (j=(int)(n*H-H+1)-1; j>=2;j--)
		 {
		   for (i=2;i<=qq[j]-1;i++)
		     {
		       	if (i>=qq[j-1])
			  u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1]+1)/(qq[j]-qq[j-1]+1);
		     }
		 }

	  //inner boundary	       
	for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	    u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
  u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	  
	}
	//	u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      
      //	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
      u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));

       u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[j]-y[1])/((x[a]-x[a-b]))))));
	}

	mgfas4new(u,n,ncycle,m,a,b,qq,x,y, R0,K, H);//solve PDE, not calculate the upper and curved boundary
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc3new(f,u,n,qq,kk,b,R0,K,H);//integrand, inner region and lower boundary
		
n1=1;

		z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, qq[(int)(n*H-H+1)-b]-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);
		zero=Stat2new(u,n,a,b,R0,K,H,qq);//error, not calculate the upper and curved boundary


		if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
			}

		printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
		//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);

		// ring defect with different radii
		for(kk=1+2*b;kk<=qq[1]-b-2;kk++) 
	{	
	 a=kk;
	

	u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
        	qq=vectornew(1,(int)(n*H-H+1));
	
	
	for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
		
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
		
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}
	
	//initial and outer boundary conditions
	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			
		  q[j]=R0-h1*pow(H,h2)+h1*pow(1.0*(j-1)/(n-1),h2);
		  qq[j]=round(n*q[j]-q[j]+1);
	if (qq[j]>(int)(n*R0-R0+1))//something new
	  {
	    qq[j]=(int)(n*R0-R0+1);
	  }
		 
		}

	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }


for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
 
for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=qq[j]-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
 
	for (i=1;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

 for (j=2;j<=(int)(n*H-H+1)-1;j=j+1)
		{
		       if (q[j-1]==q[j])
			   {
			   u[qq[j]][j]=Right;
			   }
			   else
			   {
			    
			   u[qq[j]][j]=Right-0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));
			 
			   }

                           u[1][j]=Left-0.5*m*Pi;
			}

  for (j=(int)(n*H-H+1)-1; j>=2;j--)
		 {
		   for (i=2;i<=qq[j]-1;i++)
		     {
		       	if (i>=qq[j-1])
			  u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1]+1)/(qq[j]-qq[j-1]+1);
		     }
		 }

  //inner boundary conditions
	for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	    u[a+b][j]=0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	   u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan((y[1+b]-y[1])/((x[i]-x[a])))));
	}
	//	u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
		u[a][1+b]=0.5*0.5*Pi+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      // u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));

       u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
          u[a-b][j]=-0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan((y[j]-y[1])/((x[a]-x[a-b]))))));
	}



	mgfas4new(u,n,ncycle,m,a,b,qq,x,y, R0,K, H);//solve PDE
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc3new(f,u,n,qq,kk,b,R0,K,H);//integrand
		
	n1=1;
	z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, qq[(int)(n*H-H+1)-b]-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);//total energy

       

	zero=Stat2new(u,n,a,b,R0,K,H,qq);//error

	if (z[kk]<result)//find energy minimum
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);

		fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);
	}
}


	//radial type
	
	if (m==1.0)
	  {
	    //point defect
	  for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
	       
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
        
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}
	//initial and outer boundary conditions
	for (j=1;j<=(int)(n*H-H+1);j++)
		{
			
		  q[j]=R0-h1*pow(H,h2)+h1*pow(1.0*(j-1)/(n-1),h2);
		  qq[j]=round(n*q[j]-q[j]+1);
	if (qq[j]>(int)(n*R0-R0+1))//something new
	  {
	    qq[j]=(int)(n*R0-R0+1);
	  }
		 
		}

 
for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }
	
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=qq[j]-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=b;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+b;i<=qq[1];i++)
		u[i][1]=Down2;

	for (j=2;j<=(int)(n*H-H+1)-1;j=j+1)
		{
		       if (q[j-1]==q[j])
			   {
			   u[qq[j]][j]=Right;
			   }
			   else
			   {
			    
		          u[qq[j]][j]=Right+0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));
			     
			   }
			 u[1][j]=Left-0.5*m*Pi;
			}
	//inner boundary
	       for (j=(int)(n*H-H+1)-1; j>=2;j--)
		 {
		   for (i=2;i<=qq[j]-1;i++)
		     {
		       	if (i>=qq[j-1])
			  u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1]+1)/(qq[j]-qq[j-1]+1);
		     }
		 }
	for (j=2;j<=1+b;j++)
	{
		u[1+b][j]=-atan((y[j]-y[1])/((x[1+b]-x[1])))+0.5*Pi;
		
	}
	for (i=1;i<=1+b;i++)
	{
		u[i][1+b]=-atan((y[1+b]-y[1])/((x[i]-x[1])))+0.5*Pi;
		
	}


        mgfas5new(u,n,ncycle,m,b,qq,x,y, R0,K, H);//solve PDE
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc4new(f,u,n,qq,b,R0,K,H);//integrand
		
n1=1;


 z[1]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, n1,n1+b,n1,b+1,n, H, R0)-integratnew(f, qq[(int)(n*H-H+1)-b]-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);//total energy

 zero=Stat3new(u,n,b,R0,K, H,qq);//error

	result=z[1];
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	radius=1;
	
		printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
		//	fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);
	fprintf(fp,"{0, %f},\t",z[1]);

//Directors on the boundary
		for (j=2;j<=(int)(n*H-H+1)-2;j=j+1)
			{
		      
			    
			  // printf("Arrow[{{%f,%f},{%f,%f}}],\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1),1.0*(qq[j]-2)/(n-1)+0.07*sin(u[qq[j]-1][j]),1.0*(j-1)/(n-1)+0.07*cos(u[qq[j]-1][j]));
			  //	 fprintf(fp,"Arrow[{{%f,%f},{%f,%f}}],\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1),1.0*(qq[j]-2)/(n-1)+0.07*sin(u[qq[j]-1][j]),1.0*(j-1)/(n-1)+0.07*cos(u[qq[j]-1][j]));
		}
		//print the curved boundary
		// for (j=2;j<=(int)(n*H-H+1)-2;j=j+1)
		//	{
		//  printf("{%f,%f},\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1));
		// fprintf(fp,"{%f,%f},\n",1.0*(qq[j]-2)/(n-1),1.0*(j-1)/(n-1));
		//	}


		a=1+2*b;//ring defect with smallest radius
		kk=a;
		//initial and outer boundary conditions
for (j=1;j<=(int)(n*H-H+1);j++)
		{
			
		  q[j]=R0-h1*pow(H,h2)+h1*pow(1.0*(j-1)/(n-1),h2);
		  qq[j]=round(n*q[j]-q[j]+1);
 	if (qq[j]>(int)(n*R0-R0+1))//something new
	  {
	    qq[j]=(int)(n*R0-R0+1);
	  }
		 
		}

	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }
	


	for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
	
for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=qq[j]-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
 
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

	



	       for (j=2;j<=(int)(n*H-H+1)-1;j=j+1)
		{
		       if (q[j-1]==q[j])
			   {
			   u[qq[j]][j]=Right;
			   }
			   else
			   {
			    
		          u[qq[j]][j]=Right+0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));
			   
			 
			   }

                           u[1][j]=Left-0.5*m*Pi;
			}
	
	  for (j=(int)(n*H-H+1)-1; j>=2;j--)
		 {
		   for (i=2;i<=qq[j]-1;i++)
		     {
		       	if (i>=qq[j-1])
			  u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1]+1)/(qq[j]-qq[j-1]+1);
		     }
		 }

	  //inner boundary
for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	  	  u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));

	  	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));
	}
	//	u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(0.5*Pi));
		u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
       u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
      	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
	}



	
	mgfas4new(u,n,ncycle,m,a,b,qq,x,y, R0,K, H);//solve PDE
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc3new(f,u,n,qq,kk,b,R0,K,H);//integrand
		
n1=1;
 z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, qq[(int)(n*H-H+1)-b]-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);//total energy
 zero=Stat2new(u,n,a,b,R0,K,H,qq);//error

 if (z[kk]<result)//find energy minimum
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}

		printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
		//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);
		
		//ring defect with different radii
			for(kk=1+2*b;kk<=qq[1]-b-1;kk++)
	{
	 a=kk;
	

	u=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	q=dvectornew(1,(int)(n*H-H+1));
	
	
		for (i=1;i<=(int)(n*R0-R0+1);i++)
	{
		
		x[i]=R0*(i-1)*(1.0/((int)(n*R0-R0+1)-1));
	}
	for (i=1;i<=(int)(n*H-H+1);i++)
	{
		
		y[i]=H*(i-1)*(1.0/((int)(n*H-H+1)-1));
	}
	//initial and outer boundary conditions
for (j=1;j<=(int)(n*H-H+1);j++)
		{
			
		  q[j]=R0-h1*pow(H,h2)+h1*pow(1.0*(j-1)/(n-1),h2);
		  qq[j]=round(n*q[j]-q[j]+1);
	if (qq[j]>(int)(n*R0-R0+1))//something new
	  {
	    qq[j]=(int)(n*R0-R0+1);
	  }
		 
		}


	for(i=1;i<=(int)(n*R0-R0+1);i++)
 for(j=1;j<=(int)(n*H-H+1);j++)
  {
    u[i][j]=0.0;
  }



for(i=1;i<=(int)(n*R0-R0+1);i++)
		{
			    u[i][(int)(n*H-H+1)]=Up-0.5*m*Pi;
	    }
for(j=2;j<=(int)(n*H-H+1)-1;j++)
		{
			for(i=2;i<=qq[j]-1;i++)
			    u[i][j]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=(int)(n*R0-R0+1);i++)
		u[i][1]=Down2;

 for (j=2;j<=(int)(n*H-H+1)-1;j=j+1)
		{
		       if (q[j-1]==q[j])
			   {
			   u[qq[j]][j]=Right;
			   }
			   else
			   {
			    
			      u[qq[j]][j]=Right+0.5*m*Pi-atan(1.0/((double)((n-1)*(q[j]-q[j-1]))));
			 
			   }

                           u[1][j]=Left-0.5*m*Pi;
			}

  for (j=(int)(n*H-H+1)-1; j>=2;j--)
		 {
		   for (i=2;i<=qq[j]-1;i++)
		     {
		       	if (i>=qq[j-1])
			  u[i][j]=u[qq[j-1]][j-1]+(u[qq[j]][j]-u[qq[j-1]][j-1])*(i-qq[j-1]+1)/(qq[j]-qq[j-1]+1);
		     }
		 }
  //inner boundary
	for (j=2;j<=1+b;j++)
	{
	  //	u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	  u[a+b][j]=-0.5*atan((y[j]-y[1])/((x[a+b]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[j]-y[1])/((x[a+b]-x[a])))));
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));
	  	u[i][1+b]=-0.5*atan((y[1+b]-y[1])/((x[i]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan((y[1+b]-y[1])/((x[i]-x[a])))));
	}
	//	u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(0.5*Pi));
		u[a][1+b]=-0.5*0.5*Pi+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
        u[i][1+b]=0.5*atan((y[1+b]-y[1])/((x[a]-x[i])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin((Pi-atan((y[1+b]-y[1])/((x[a]-x[i]))))));
	}
	for (j=2;j<=1+b;j++)
    {
      //	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
	u[a-b][j]=0.5*atan((y[j]-y[1])/((x[a]-x[a-b])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan((y[1+b]-y[1])/((x[a]-x[a-b]))))));
      
	}

        mgfas4new(u,n,ncycle,m,a,b,qq,x,y, R0,K, H);//solve PDE
		f=dmatrix2new(1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
		
		nrfunc3new(f,u,n,qq,kk,b,R0,K, H);//integrand
		
	n1=1;


	z[kk]=integratnew(f, n1,(int)(n*R0-R0+1),n1,(int)(n*H-H+1),n, H, R0)-integratnew(f, a-b,a+b,n1,b+1,n, H, R0)-integratnew(f, qq[(int)(n*H-H+1)-b]-b,(int)(n*R0-R0+1),(int)(n*H-H+1)-b,(int)(n*H-H+1),n, H, R0);//total energy

	zero=Stat2new(u,n,a,b,R0,K,H,qq);//error

	if (z[kk]<result)//find energy minimum
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=(int)(n*R0-R0+1);i++)
		for(j=1;j<=(int)(n*H-H+1);j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	//	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
		fprintf(fp,"{%f, %f},\t",1.0*(kk-1)/(n-1),z[kk]);
	}
	}
	


	printf("Total energy=%f\t\nRing radius=%d\n",result,radius);
	fprintf(fp,"Number of lattices vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\nDefect type m=%f\nAspect ratio R0=%f\nFrank constant ratio K=%f\nTotal energy=%f\nRing radius=%d\n",n,b,m,R0,K,result,radius);		
	

	

	fclose (fp);
                   
	   

	
		
		

		fprintf(fp2,"{");fprintf(fp3,"{");
	for(i=1;i<=(int)(n*R0-R0+1)-1;i++)
		{
			fprintf(fp2,"{");fprintf(fp3,"{");
			for(k=1;k<=(int)(n*H-H+1)-1;k++)
			  {
			    fprintf(fp2,"{%f,\t%f},\t",sin(v[i][k]),cos(v[i][k]));
			    if (i==qq[k])
			      {
				fprintf(fp3,"{1.0,0},\t");
			      }
			    else
			      {
				fprintf(fp3,"{0,0},\t");
			      }
			  }
			fprintf(fp2,"{%f,\t%f}\t",sin(v[i][(int)(n*H-H+1)]),cos(v[i][(int)(n*H-H+1)]));
		        
			    if (i==qq[k])
			      {
				fprintf(fp3,"{1.0,0}\t");
			      }
			    else
			      {
				fprintf(fp3,"{0,0}\t");
			      }
			    fprintf(fp2,"},");fprintf(fp3,"},");
		}
	fprintf(fp2,"{");fprintf(fp3,"{");
			for(k=1;k<=(int)(n*H-H+1)-1;k++)
			  {
				fprintf(fp2,"{%f,\t%f},\t",sin(v[(int)(n*R0-R0+1)][k]),cos(v[(int)(n*R0-R0+1)][k]));
					fprintf(fp3,"{0,0},\t");
			  }
			fprintf(fp2,"{%f,\t%f}\t",sin(v[(int)(n*R0-R0+1)][(int)(n*H-H+1)]),cos(v[(int)(n*R0-R0+1)][(int)(n*H-H+1)]));
				fprintf(fp3,"{0,0}\t");
			fprintf(fp2,"}");
			fprintf(fp3,"}");
		fprintf(fp2,"}");
			fprintf(fp3,"}");
			fclose (fp2);fclose (fp3);
		
    printf ("File created okay\n");

	
	free_dmatrix2new(u,1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	free_dmatrix2new(v,1,(int)(n*R0-R0+1),1,(int)(n*H-H+1));
	free_dvectornew(q,1,(int)(n*H-H+1));
       	free_vectornew(qq,1,(int)(n*H-H+1));
	free_dvectornew(x,1,(int)(n*R0-R0+1));
	free_dvectornew(y,1,(int)(n*H-H+1));
	free_dvectornew(z,1,(int)(n*R0-R0+1));
	free_dvectornew(z2,1,(int)(n*H-H+1));
    return 0;
}
