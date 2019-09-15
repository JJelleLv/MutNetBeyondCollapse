#ifndef CEULER1_H_INCLUDED
#define CEULER1_H_INCLUDED
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include "matmatrix.h"
#include "onedimtools.h"


using namespace std;
//enum solverkind {euler, rk4};

class TExternvar
{
public:
   TMatrix data; //must be sorted
   bool cycle = false;
   bool tofloor = false;  //do not interpolate within time steps
   bool active = true;
   double prevt = 0;
   int previ = 0;
   double defaultvar = 0.0;
 //  TMatrix defaultvar;
   double valueat (double t);
};

class TBoxcar
{
  string name="";
  double gcycl =0.0;
  double outflow =0.0;
  TMatrix flow;
};

class TModel;

typedef void (TModel:: *solverfunction) ();

class TModel
{
private:
   vector<double> dydt;  //derivatives
   vector<double> y0;        //current y0 during run
   vector<TExternvar> rednoisesets; //used by rednoise to store the rednoise data set
   ofstream outfile;
   void writeoutput(double *y,int n=1);
   void parseline (string &line, size_t &p, unsigned int &ndx, unsigned int &nx, unsigned int &ny, double &valu);
   void interpstep (double t1, double t2, vector<double>* X1, vector<double>* X2, double tnew, vector<double> g_perm);
   void readparameters (string filename); //read parameters (binary)
   void openoutfile();
   void runode45();
   void ntrp45 (double tinterp, double t, double * y, unsigned int neq, double h, double * f1, double * f3, double * f4, double * f5, double * f6, double * f7);
   void runrk4 ();
   void runeuler ();
   void rundiffer ();
   void singlerun();
public:
   vector<double> param;     //for a non vector model this is used for parameters
   vector<double> permanent; //permanent variables (are stored and may change during the run)
   vector<double> y_start;  //initial conditions
   vector<double> t_out;     //times for output
   vector<double> externvalues;     //current values of the external variables
   vector<TMatrix> v_y_start; //for vector models: initial conditions
   vector<TMatrix> v_y0; //for vector models: current value
   vector<TMatrix> v_dydt; //for vector models: differential
   vector<TMatrix> v_param; //for vector models: parameters
   vector<TMatrix> v_permanent; //for vector models: permanent variables
   vector<TMatrix> v_auxil; //for vector models: auxiliary variables
   vector<TExternvar> externvars; //external variables
   double tstart = 0;
   double nout = -1;
   double tend = 1000;
   double ndays = 1000;
   bool refine = true;  //refine the results for ode45 (if the nout is not defined)
   bool includet =false;
   bool backward = false;   //true for run backwards
   bool stats = false;      //true if show stats
   double h = 0.1;          //time step for fixed stepping
   double atol = 1e-6;      //absolute tolerance ODE45
   double rtol = 1e-3;      //relative tolerance ODE45
   bool isdiffer = false;   //is true if the model is a difference equation
   bool isvector = false;   //vector mode not yet implemented
   vector<int>  nonnegative; //Non-negative values are forced
   int iters = 1;
   vector<double> boxcar_gcycl;
   vector<string> boxcar_name;
   solverfunction  solverfcn = &TModel::runode45; //ode45 is the default solver
   string outfilename = "out1.tmp";
   TModel (string filename); //Constructor reads parameter file
   double rednoise (double t, double T0, double lambda, double beta, int iset = 1, double deltat = 1);
   double dwiener (double gY, double dgY_dY = 0);
   boxcar_result boxcartrain (int boxnr, TMatrix & A, double devrate, double cv_devrate = 0);
   double solver (string par);
   void  runodefile (double t, double * _ystart,  double * _dydt);
   void runsolver ();
   void setdiffer (bool d);
   void setstatevars (int n);
   void setparameters (int n);
   void setpermanent (int n);
   void setauxil (int n);
};



#endif // CEULER1_H_INCLUDED
