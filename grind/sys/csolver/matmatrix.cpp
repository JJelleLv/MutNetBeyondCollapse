#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "matmatrix.h"


using namespace std;

TMatrix::TMatrix() //constructor
{
   setsize (0, 0);
}


TMatrix::TMatrix (double m) //constructor
{
   setsize (1, 1);
   data[0] = m;
}
TMatrix::TMatrix (int m) //constructor
{
   setsize (1, 1);
   data[0] = m;
}
double & TMatrix::operator() (int row, int col)
{
   return data[sub2ind (row, col)];
}


double & TMatrix::operator() (int elem)
{
   if (numel == 1)
      return data[0];
   else
      return data[elem];
}


TMatrix::TMatrix (int arows, int acols) //constructor
{
   setsize (arows, acols);
}
double zeros(int ndx,double A, double B)
{
    return(0.0);
}
double ones(int ndx,double A, double B)
{
    return(1.0);
}
void TMatrix::setsize (int arows, int acols)
{
   nrows = arows;
   ncols = acols;
   numel = (arows * acols);
   if (data == Fdata.data())
   {
      Fdata.resize (numel, nan (""));
      data = Fdata.data();
   }
}
void TMatrix::ind2sub (int &row, int &col, int ndx)
{
   row = fmod (ndx, nrows) ;
   col = (ndx - row) / nrows ;
}

void TMatrix::str2data (string &line, size_t p = 0)
//read data from string, start at position p+1, or start at first position
{
   size_t p1;
   if (p == 0)
      p1 = line.find (" ", p);
   else
      p1 = line.find (" ", p + 1);
   int nx = atoi (line.substr (p + 1, p - p1).c_str());
   p = line.find (" ", p1 + 1);
   int ny = atoi (line.substr (p1 + 1, p1 - p).c_str());
   setsize (nx, ny);
   for (int i = 0; i < numel - 1; i++)
   {
      p1 = line.find (" ", p + 1);
      data[i] = atof (line.substr (p + 1, p - p1).c_str());
      p = p1;
   }
   data[numel - 1] = atof (line.substr (p + 1).c_str());
}


double TMatrix::getround (int row, int col, int border, double avalue)
{
   switch (border)
   {
   case 0:
      if (col < 0) col = ncols - 1;
      if (col >= ncols) col = 0;
      if (row < 0) row = nrows - 1;
      if (row >= nrows) row = 0;
      return (*this) (row, col);
      break;
   case 1:
      if (col < 0) col = 0;
      if (col >= ncols) col = ncols - 1;
      if (row < 0) row = 0;
      if (row >= nrows) row = nrows - 1;
      return (*this) (row, col);
      break;
   case 2:
      if ( (col < 0) || (col >= ncols) || (row < 0) || (row >= nrows))
         return avalue;
      else
         return (*this) (row, col);
      break;
   default:
      return nan ("");
   }
}

void TMatrix::print (int rows, int cols)
{
   if (rows <= 0) rows = nrows;
   if (cols <= 0) cols = ncols;
   for (int i = 0; i < rows; i++)
   {
      for (int j = 0; j < cols; j++)
      {
         cout << (*this) (i, j) /*data[sub2ind (i, j)] */ << " ";
      }
      cout << endl;
   }
}

//void mtimes (TMatrix &result, TMatrix &A, TMatrix &B)
////Full Matrix multiplication A*B
//{
//   result.setsize (A.nrows, B.ncols);
//   double sum;
//   for (int i = 0; i < A.nrows; i++)
//   {
//      int j;
//      for (j = 0; j < B.ncols; j++)
//      {
//         sum = 0;
//         for (int k = 0; k < A.ncols; k++)
//         {
//            sum = sum + A (i, k) * B (k, j);
//         }
//         result (i, j) = sum;
//      }
//   }
//}

//MATLAB matrix functions that can be translated to a function that returns only one element
double mtimes (int ndx, TMatrix &A, TMatrix &B)
//Matrix multiplication A*B calculating only one element
{
   if ( (B.numel == 1) || (A.numel == 1)) // if one of both args is a scalar return the dot product (like MATLAB)
      return A (ndx) * B (ndx);
   else
   {
      int row = fmod (ndx, A.nrows) ;
      int col = (ndx - row) / A.nrows ;
//     cout<<A.nrows<<' '<<B.ncols<<' '<<ndx<<' '<<row<<' '<<col<<endl;
      double result = 0;
      for (int k = 0; k < A.ncols; k++)
      {
         result = result + A (row, k) * B (k, col);
      }
      return result;
   }
}
double mtimes (int ndx, double &A, TMatrix &B)
{
   return A * B (ndx);
}
double mtimes (int ndx, double &A, double &B)
{
   return A * B;
}
double mtimes (int ndx, TMatrix &A, double &B)
{
   return A (ndx) * B;
}
double mtimes1d (int ndx, TMatrix &A, TMatrix &B)
//Matrix multiplication A*B calculating only one element, little bit more optimized
{
   if ( (B.numel == 1) || (A.numel == 1)) // if one of both args is a scalar return the dot product (like MATLAB)
      return A (ndx) * B (ndx);
   else
   {
      double result = 0;
      for (int k = 0; k < A.ncols; k++)
      {
         result = result + A (ndx, k) * B (k);
      }
      return result;
   }
}
double mtimes1d (int ndx, double &A, TMatrix &B)
{
   return A * B (ndx);
}
double mtimes1d (int ndx, double &A, double &B)
{
   return A * B;
}
double mtimes1d (int ndx, TMatrix &A, double &B)
{
   return A (ndx) * B;
}

double g_sum (int ndx, TMatrix &A, int dim)
//sum along one diminsion
{
   double result = 0;
   if ( ( (dim == 1) && (A.nrows > 1)) || ( (dim == 2) && (A.ncols == 1)))
   {
      int col = min (A.ncols - 1, ndx);
      for (int row = 0; row < A.nrows; row++)
         result = result + A (row, col);
   }
   else
   {
      int row = min (A.nrows - 1, ndx);
      for (int col = 0; col < A.ncols; col++)
         result = result + A (row, col);
   }
   return result;
}

double mean (int ndx, TMatrix &A, int dim)
//sum along one diminsion
{
   double result = 0;
   int n;
   if ( ( (dim == 1) && (A.nrows > 1)) || ( (dim == 2) && (A.ncols == 1)))
   {
      int col = min (A.ncols - 1, ndx);
      n = A.nrows;
      for (int row = 0; row < A.nrows; row++)
         result = result + A (row, col);
   }
   else
   {
      int row = min (A.nrows - 1, ndx);
      n = A.ncols;
      for (int col = 0; col < A.ncols; col++)
         result = result + A (row, col);
   }
   return result / double (n);
}

double transpose (int ndx, TMatrix &A)
{
   int row = fmod (ndx, A.ncols) ;
   int col = (ndx - row) / A.ncols ;
//   row = fmod (ndx, nrows) ;
//  col = (ndx - row) / nrows ;
   return A (col, row);
}
double repmat (int ndx, TMatrix &A, int tilenrows, int tilencols)
{
   int nrows = A.nrows * tilenrows;
//  int ncols = A.ncols * tilencols;
   int row = fmod (ndx, nrows) ;
   int col = (ndx - row) / nrows;
   row = fmod (row, A.nrows);
   col = fmod (col, A.ncols);
   return (A (row, col));
}
double repmat (int ndx, TMatrix &A, int tilenrows, TMatrix &tilencols)
{
   return repmat (ndx, A, tilenrows, tilencols (0));
}
double repmat (int ndx, TMatrix &A, TMatrix &tilenrows, TMatrix &tilencols)
{
   return repmat (ndx, A, tilenrows (0), tilencols (0));
}
double repmat (int ndx, TMatrix &A, TMatrix &tilenrows, int tilencols)
{
   return repmat (ndx, A, tilenrows (0), tilencols);
}

double flipud (int ndx, TMatrix &A)
{
   int col, row;
   A.ind2sub (row, col, ndx);
   return A (A.nrows - 1 - row, col);
}
double fliplr (int ndx, TMatrix &A)
{
   int col, row;
   A.ind2sub (row, col, ndx);
   return A (row, A.ncols - 1 - col);
}
double neighborcells (int ndx, TMatrix &A, int nneighbors, int bordered)
{
   int border;
   if (bordered == 0)
      border = 0;
   else
      border = 2;
   int col, row;
   A.ind2sub (row, col, ndx);
   double value = A.getround (row, col - 1, border, 0)
                  + A.getround (row, col + 1, border, 0)
                  + A.getround (row - 1, col, border, 0)
                  + A.getround (row + 1, col, border, 0);
   if (nneighbors == 8)
   {
      value = value + A.getround (row + 1, col + 1, border, 0)
              + A.getround (row - 1, col + 1, border, 0)
              + A.getround (row + 1, col - 1, border, 0)
              + A.getround (row - 1, col - 1, border, 0);
   }
   return value;
}

double leftcells (int ndx, TMatrix &A, int border, double avalue)
{
   int col, row;
   A.ind2sub (row, col, ndx);
   return A.getround (row, col - 1, border, avalue);
}
double rightcells (int ndx, TMatrix &A, int border, double avalue)
{
   int col, row;
   A.ind2sub (row, col, ndx);
   return A.getround (row, col + 1, border, avalue);
}
double upcells (int ndx, TMatrix &A, int border, double avalue)
{
   int col, row;
   A.ind2sub (row, col, ndx);
   return A.getround (row + 1, col, border, avalue);
}
double leftcells (int ndx, TMatrix &A, int border, TMatrix &avalue)
{
   return leftcells (ndx, A, border, avalue (0));
}

double rightcells (int ndx, TMatrix &A, int border, TMatrix &avalue)
{
   return rightcells (ndx, A, border, avalue (0));
}

double downcells (int ndx, TMatrix &A, int border, TMatrix &avalue)
{
   return downcells (ndx, A, border, avalue (0));
}

double upcells (int ndx, TMatrix &A, int border, TMatrix &avalue)
{
   return upcells (ndx, A, border, avalue (0));
}

double downcells (int ndx, TMatrix &A, int border, double avalue)
{
   int col, row;
   A.ind2sub (row, col, ndx);
   return A.getround (row - 1, col, border, avalue);
}

double numel (int ndx, TMatrix &A)
{
   return A.numel;
}
double size (int ndx, TMatrix &A, int dim)
{
   if ( (dim < 1) && (ndx < 2)) dim = ndx;
   if (dim == 1)
      return (A.nrows);
   else if (dim == 2)
      return (A.ncols);
   else
      return nan ("");
}
double size (TMatrix &A, int dim)
{
   return size (0, A, dim);
}
double length (int ndx, TMatrix &A)
{
   return max (A.ncols, A.nrows);
}
double g_min (int ndx, TMatrix &A, TMatrix &B) //MATLAB's min should be replaced by g_min
{
   return min (A (ndx), B (ndx));
}
double g_max (int ndx, TMatrix &A, TMatrix &B) //MATLAB's max should be replaced by g_max
{
   return max (A (ndx), B (ndx));
}
double g_min (int ndx, double A, TMatrix &B, int dim)
{
   return g_min (ndx, B, A);
}
double g_max (int ndx, double A, TMatrix &B, int dim)
{
   return g_max (ndx, B, A);
}
double g_min (int ndx, TMatrix &A, double B, int dim)
//min along one diminsion
{
   if (B != g_emptyvar)
      return min (A (ndx), B);
   else
   {
      double result = 1E100;
      if ( ( (dim == 1) && (A.nrows > 1)) || ( (dim == 2) && (A.ncols == 1)))
//   if (dim == 1)
      {
         int col = ndx;
         for (int row = 0; row < A.nrows; row++)
            if (result > A (row, col)) result = A (row, col);
      }
      else
      {
         int row = ndx;
         for (int col = 0; col < A.ncols; col++)
            if (result > A (row, col)) result = A (row, col);
      }
      return result;
   }
}

double boxcarinflow (int ndx, TMatrix &A, double inflow)
{
   if (ndx == 0)
      return (inflow);
   else
      return (0.0);
}
double boxcarinflow (int ndx, TMatrix &A, TMatrix &inflow)
{
   if (ndx == 0)
      return (inflow (0));
   else
      return (0.0);
}

double g_max (int ndx, TMatrix &A, double B, int dim)
//max along one diminsion
{
   if (B != g_emptyvar)
      return max (A (ndx), B);
   else
   {
      double result = -1E100;
      if ( ( (dim == 1) && (A.nrows > 1)) || ( (dim == 2) && (A.ncols == 1)))
//   if (dim == 1)
      {
         int col = ndx;
         for (int row = 0; row < A.nrows; row++)
            if (result < A (row, col)) result = A (row, col);
      }
      else
      {
         int row = ndx;
         for (int col = 0; col < A.ncols; col++)
            if (result > A (row, col)) result = A (row, col);
      }
      return result;
   }
}
double cumsum (int ndx, TMatrix &A, int dim)
//cumulative sum along one diminsion
{
   int col, row;
   A.ind2sub (row, col, ndx);
   double result = 0;
   if ( ( (dim == 1) && (A.nrows > 1)) || ( (dim == 2) && (A.ncols == 1)))
//  if (dim == 1)
   {
      for (int k = 0; k < A.nrows; k++)
         result = result + A (k, col);
   }
   else
   {
      for (int k = 0; k < A.ncols; k++)
         result = result + A (row, k);
   }
   return result;
}

double prod (int ndx, TMatrix &A, int dim)
//product along one diminsion
{
   double result = 1;
   if ( ( (dim == 1) && (A.nrows > 1)) || ( (dim == 2) && (A.ncols == 1)))
//   if (dim == 1)
   {
      int col = ndx;
      for (int row = 0; row < A.nrows; row++)
         result = result * A (row, col);
   }
   else
   {
      int row = ndx;
      for (int col = 0; col < A.ncols; col++)
         result = result * A (row, col);
   }
   return result;
}

double cumprod (int ndx, TMatrix &A, int dim)
//cumulative prod along one diminsion
{
   int col, row;
   A.ind2sub (row, col, ndx);
   double result = 1;
   if (dim == 1)
   {
      for (int k = 0; k < A.nrows; k++)
         result = result * A (k, col);
   }
   else
   {
      for (int k = 0; k < A.ncols; k++)
         result = result * A (row, k);
   }
   return result;
}

double eye (int ndx, int nrows, int ncols)
{
   if (ncols <= 0) ncols = nrows;
   int row = fmod (ndx, nrows) ;
   int col = (ndx - row) / nrows ;
   if (col == row)
      return 1.0;
   else
      return 0.0;
}
double linspace (int ndx, double amin, double amax, int n)
{
   if (ndx < n)
      return amin + ndx * (amax - amin) / (n - 1);
   else
      return nan ("");
}
double logspace (int ndx, double amin, double amax, int n)
{
   if (ndx < n)
      return pow (10.0, amin + ndx * (amax - amin) / (n - 1));
   else
      return nan ("");
}
double parlookup (int ndx, TMatrix &tabl, TMatrix &key1, bool extrapolate, bool NaNiszero)
//function par = parlookup(tabl, key1, extrapolate, NaNisZero)
//if nargin < 4
//   NaNisZero = 1;
//end;
//[N, C] = size(tabl);
//if key1 < tabl(1, 1) || key1 > tabl(N, 1)
//   if (nargin < 3) || ~extrapolate
//      par = NaN .* zeros(1, size(tabl, 2) - 1);
//   elseif key1 < tabl(1, 1)
//      par = tabl(1, 2:C);
//   else
//      par = tabl(N, 2:C);
//   end;
//else
//   i = 1;
//   while (i<N) && tabl(i, 1) <= key1
//      i = i + 1;
//   end;
//   a = (tabl(i, 2:C) - tabl(i - 1, 2:C)) ./ (tabl(i, 1) - tabl(i - 1, 1));
//   b = tabl(i - 1, 2:C) - a .* tabl(i - 1, 1);
//   par = a .* key1 + b;
//end;
//if NaNisZero
//   par(isnan(par)) = 0;
//end;
{
   double key = key1 (0);
   double par = 0;
   ndx = ndx + 1;
   if (key < tabl (0, 0) || key > tabl (tabl.nrows - 1, 0))
   {
      if (~extrapolate)
         par = nan ("");
      else if (key < tabl (0, 0))
         par = tabl (0, ndx);
      else if (key > tabl (tabl.nrows - 1, 0))
         par = tabl (tabl.nrows - 1, ndx);
   }
   else
   {
      int i = 0;
      while ( (i < tabl.nrows - 1) && (tabl (i, 0) <= key))
      {
         i++;
      }
      double a = (tabl (i, ndx) - tabl (i - 1, ndx)) / (tabl (i, 0) - tabl (i - 1, 0));
      double b = tabl (i - 1, ndx) - a * tabl (i - 1, 0);
      par = (a * key + b);
   }
   if (NaNiszero && isnan (par))
      return (0);
   else
      return (par);
}

