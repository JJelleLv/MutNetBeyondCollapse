/*=================================================================
 *
 *  GRIND MEX ODE FILE
 *  INCLUDE THIS FILE IN YOUR MEX PROJECT
 * 
 *=================================================================*/
#include <math.h>
#include "mex.h"

//*********************************************
// mexGetdouble gets a global variable assuming that it is a scalar
//
double mexGetdouble(char* name)
{
   double value;
   mxArray *array_ptr;
#ifdef R11
   array_ptr = mexGetArrayPtr(name, "global");
#else
   array_ptr = mexGetVariablePtr("global",name);
#endif
   value = *mxGetPr(array_ptr);
   return(value);
}
typedef struct {
   int nrows;
   int ncols;
   double *mat;
} mymatrix;

void mexGetmymatrix(mymatrix *result, char *name)
{
   mxArray *array_ptr;
#ifdef R11
   array_ptr = mexGetArrayPtr(name, "global");
#else
   array_ptr = mexGetVariablePtr("global", name);
#endif
   result->mat = mxGetPr(array_ptr);
   result->nrows = mxGetN(array_ptr);
   result->ncols = mxGetM(array_ptr);
   return;
}

double *rmat(mymatrix *m, int row, int col)
{
if (((*m).nrows==1) && ((*m).ncols==1)) return ((*m).mat);

if (row<0) 
  row=row+(*m).nrows;
if (col<0) 
   col=col+(*m).ncols;
if (row>(*m).nrows-1) 
   row=row-(*m).nrows;
if (row>(*m).nrows-1) 
   row=row-(*m).nrows;
return ((*m).mat+row*(*m).nrows+col);

}

double *mat(mymatrix *m, int row, int col)
{
if (((*m).nrows==1) && ((*m).ncols==1)) return ((*m).mat);
return ((*m).mat+row*(*m).nrows+col);

}



double diffusion4(mymatrix *m, double d,int row,int col)
{
return (-d*(4*(*rmat(m,row,col))-
(*rmat(m,row-1,col))-(*rmat(m,row+1,col))-
(*rmat(m,row,col+1))-(*rmat(m,row,col-1))));
}


//*********************************************
// basic wrapper for ODE mex function
//
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray  *prhs[] )
     
{ 
    double *dydt; 
    double *t,*y; 
    unsigned int m,n; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
   
    m = mxGetM(prhs[1]); 
    n = mxGetN(prhs[1]);

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    dydt = mxGetPr(plhs[0]);
    
    t = mxGetPr(prhs[0]); 
    y = mxGetPr(prhs[1]);

    /* Do the actual computations in a subroutine */
    yprime(dydt,t,y); 
    return;
    
}

