/*=================================================================
*
*  GRIND MEX ODE FILE
*  include "grindmex.h" for all declarations
*  define here the yprime function
*
*=================================================================*/
#include "grindmex.h"

static void yprime(double g_X2[], double *t, double g_X1[])
{
// define state variables
#define x g_X1[0]

// define primes
#define dxdt g_X2[0]

// define parameters 
double r = mexGetdouble("r");

// define local variables

// differential equations
dxdt=r*x*(1-x);
return;
}

