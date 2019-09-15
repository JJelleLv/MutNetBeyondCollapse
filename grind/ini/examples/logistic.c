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
#define V g_X1[0]

// define primes
#define dVdt g_X2[0]

// define parameters 
double K = mexGetdouble("K");
double h = mexGetdouble("h");
double hv = mexGetdouble("hv");
double r = mexGetdouble("r");

// define local variables
double prod;
double harvest;

// differential equations
//% Logistic growth with Holling type I harvesting;
//% example to demonstrate optimpars command;
//%;
//%n_=rednoise(t,0,lambda,beta);;
prod=V*r*(1-V/K);
harvest=h*(V/(V+hv));
dVdt =prod-harvest;
return;
}

