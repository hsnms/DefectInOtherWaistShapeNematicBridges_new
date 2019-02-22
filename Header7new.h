#ifndef HEADER7NEW_H
#define HEADER7NEW_H
void mgfas4new(double **u, int n,int ncycle, double m, int a,int b,int  *q,double *x, double *y, double R0, double K, double H);
void mgfas5new(double **u, int n,int ncycle, double m,int b, int *qq, double *x, double *y, double R0, double K, double H);
void **nrfunc3new(double **f,double **u, int n,int *qq, int a,int b,double R0, double K, double H);
void **nrfunc4new(double **f,double **u, int n,int *qq,int b,double R0, double K, double H);
double Stat2new (double **u,int n, int a,int b, double R0, double K, double H, int *qq);
double Stat3new (double **u,int n,int b, double R0, double K, double H, int *qq);
#endif
