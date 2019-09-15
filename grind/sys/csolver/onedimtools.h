#ifndef ONEDIMTOOLS_H_INCLUDED
#define ONEDIMTOOLS_H_INCLUDED
#include <cmath>
double rand(int i,int j=0);
double eq(double x,double y);
double ne(double x,double y);
double real(double x); //useless dummy function for compatibility with MATLAB
double sum(double x); // useless dummy function for  compatibility with MATLAB
double iif (bool cond, double res1, double res2);
double randlogn(int n1,int n2,double mean,double sd); //draw from log normal distribution n1 and n2 must be one
double ln(double x);
double randn(int i=1,int j=0);
double eps(double x=1);
double zeros(int row,int col);
double ones(int row,int col);
int sign(double x);
const double pi = 2.0*acos(0.0);

#endif // ONEDIMTOOLS_H_INCLUDED
