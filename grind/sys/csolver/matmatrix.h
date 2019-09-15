#ifndef MATMATRIX_H
#define MATMATRIX_H
#include <vector>
#include <cmath>

using namespace std;

//TMatrix can be a scalar/vector or matrix
//it's data array is 0 based (

class TMatrix  //MATLAB style matrix,
{
   vector<double> Fdata;
public:
   int nrows = 0; //never set directly
   int ncols = 0; //never set directly
   int numel = 0;
   int sizedata()
   {
      return numel;
   }; //numel = number of elements

   double * data = Fdata.data(); //pointer to an array with the data (by default this is the vector Fdata)
   //constructors:
   TMatrix();
   TMatrix (int arows, int acols);
   TMatrix (double m);
   TMatrix (int m);
   void setsize (int arows, int acols = 1); //always use this for setting the numel
   int sub2ind (int row, int col)
   {
      return col * nrows + row;
   };
   void ind2sub (int &row, int &col, int ndx);
   void str2data (string &line, size_t p);
   double getround (int row, int col, int border, double avalue = 0);
   void print (int rows = -1, int cols = -1);
   double & operator() (int elem);
   double & operator() (int row, int col);
};
struct boxcar_result
{
   double outflow;
   TMatrix flow;
};
/*
// very limited support for sparse matrices, only to be used in fast matrix multiplications
struc sprec
{
  int row;
  int col;
  double value;
}

class TSparseMatrix  //Not inherited
{
public:
   vector<sprec> nonzeros;
   int nrows = 0; //never set directly
   int ncols = 0; //never set directly
   int numel=0;
   int nnz=0; //number of nonzeros
   int sizedata() {return numel;}; //numel = number of elements
   //constructors:
   TSparseMatrix();
   TSparseMatrix (int arows, int acols, int annz=0);
   void setsize (int arows, int acols = 1,int annz=0); //always use this for setting the numel
   int sub2ind (int row, int col)
   {
      return col * nrows + row;
   };
   void ind2sub (int &row, int &col, int ndx);
   void print (int rows=-1,int cols=-1);
   double & operator() (int elem);
   double & operator() (int row, int col);
};

*/
const double g_emptyvar = -9999999;

double mtimes (int ndx, TMatrix &A, TMatrix &B);
double mtimes (int ndx, double &A, TMatrix &B);
double mtimes (int ndx, double &A, double &B);
double mtimes (int ndx, TMatrix &A, double &B);
double mtimes1d (int ndx, TMatrix &A, TMatrix &B); //slightly more efficient
double mtimes1d (int ndx, double &A, TMatrix &B);
double mtimes1d (int ndx, double &A, double &B);
double mtimes1d (int ndx, TMatrix &A, double &B);
double repmat (int ndx, TMatrix &A, int tilenrows, int tilencols);
double repmat (int ndx, TMatrix &A, int tilenrows, TMatrix &tilencols);
double repmat (int ndx, TMatrix &A, TMatrix &tilenrows, TMatrix &tilencols);
double repmat (int ndx, TMatrix &A, TMatrix &tilenrows, int tilencols);
double transpose (int ndx, TMatrix &A);
double flipud (int ndx, TMatrix &A);
double fliplr (int ndx, TMatrix &A);
double length (int ndx, TMatrix &A);
double zeros(int ndx,double A, double B);
double ones(int ndx,double A, double B);
double sum (int ndx, TMatrix &A, int dim = 1);
double leftcells (int ndx, TMatrix &A, int border = 1, double avalue = 0.0);
double leftcells (int ndx, TMatrix &A, int border, TMatrix &avalue);
double rightcells (int ndx, TMatrix &A, int border = 1, double avalue = 0.0);
double rightcells (int ndx, TMatrix &A, int border, TMatrix &avalue);
double upcells (int ndx, TMatrix &A, int border, TMatrix &avalue);
double upcells (int ndx, TMatrix &A, int border = 1, double avalue = 0.0);
double downcells (int ndx, TMatrix &A, int border = 1, double avalue = 0.0);
double downcells (int ndx, TMatrix &A, int border, TMatrix &avalue);
double numel (int ndx, TMatrix &A);
double size (TMatrix &A, int dim = -1);
double size (int ndx, TMatrix &A, int dim = -1);
double mean (int ndx, TMatrix &A, int dim = 1);
double iif (int ndx, TMatrix &cond, TMatrix &A, TMatrix &B);
double g_sum (int ndx, TMatrix &A, int dim = 1); //MATLAB's sum should be replaced by g_sum
double g_min (int ndx, double A, TMatrix &B, int dim = 1); //MATLAB's min should be replaced by g_min
double g_max (int ndx, double A, TMatrix &B, int dim = 1);
double g_min (int ndx, TMatrix &A, double B = g_emptyvar, int dim = 1);
double g_max (int ndx, TMatrix &A, double B = g_emptyvar, int dim = 1); //MATLAB's max should be replaced by g_max
double g_min (int ndx, TMatrix &A, TMatrix &B); //MATLAB's min should be replaced by g_min
double g_max (int ndx, TMatrix &A, TMatrix &B); //MATLAB's max should be replaced by g_max
double cumsum (int ndx, TMatrix &A, int dim = 1);
double prod (int ndx, TMatrix &A, int dim = 1);
double cumprod (int ndx, TMatrix &A, int dim = 1);
double eye (int ndx, int i, int j = 0);
double linspace (int ndx, double amin, double amax, int n = 100);
double logspace (int ndx, double amin, double amax, int n = 100);
double parlookup (int ndx, TMatrix &tabl, TMatrix &key1, bool extrapolate = false, bool NaNiszero = true);
double boxcarinflow (int ndx, TMatrix &A, double inflow);
double boxcarinflow (int ndx, TMatrix &A, TMatrix &inflow);

#endif // MATMATRIX_H
