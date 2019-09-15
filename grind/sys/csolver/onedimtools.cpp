#include "onedimtools.h"
#include <cmath>
#include <stdlib.h>
#include <math.h>

double ln (double x)
{
   return log (x);
}
double zeros(int row,int col)
{
    return(0.0);
}
double ones(int row,int col)
{
    return(1.0);
}
int sign (double x)
{
   return (x > 0) - (x < 0);
}

double real (double x)
{
   return x;
}

double randlogn (int n1, int n2, double mean, double sd)
{
//mean is the expected mean
//sd is expected standard dev
//randlogn draws from lognormal distribution
//if sd==0 || mean==0
//   A = ones(n1, n2) * mean;
//else
//   mu = ln(mean^2 / sqrt(sd^2 + mean^2));
//   sigma2 = ln((sd / mean)^2 + 1);
//   A = exp(mu + sigma2 * randn(n1, n2));
//end;
   if ( (n1 > 1) || (n2 > 1))
      throw ("Vector notation not allowed in a simple model");
   if ( (sd == 0) || (mean == 0))
   {
      return mean;
   }
   else
   {
      double  mu = ln ( (mean * mean) / sqrt ( (sd * sd) + (mean * mean)));
      double sigma2 = ln (pow (sd / mean, 2) + 1);
      return exp (mu + sigma2 * randn (1));
   }
}

double eps (double x)
//MATLAB function eps that gives the maximum difference between subsequent doubles
{
   x = abs (x);
   if (x < 2.2251E-308)
      return 4.9407e-324;
   else
   {

      x = pow (2, floor (log2 (x)));
      return 2.2204e-016 * x;
   }
}

double sum (double x)
{
   return x;
}

double iif (bool cond, double res1, double res2)
{
   if (cond)
   {
      return res1;
   }
   else
   {
      return res2;
   }
}
//draw from a Normal distribution
double randn (int i, int j)
// not the most efficient implementation of the Box Muller method
{
   double U1, U2;
   do
   {
      U1  = ( (double) rand() / (RAND_MAX)) ;
      U2 = ( (double) rand() / (RAND_MAX)) ;
   }
   while (U1 == 0);
   return sqrt (-2 * ln (U1)) * cos (2 * pi * U2);
}
double eq (double x, double y)
{
   return abs (x - y) < 1E-10;
}
double ne (double x, double y)
{
   return abs (x - y) > 1E-10;
}
double rand (int i, int j)
{
   return (double) rand() / (RAND_MAX);
}
