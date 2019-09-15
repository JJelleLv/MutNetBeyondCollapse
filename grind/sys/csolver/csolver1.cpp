#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include "matmatrix.h"
#include "onedimtools.h"
#include "csolver1.h"


using namespace std;

#define power pow
#define solver g_model.solver
#define min min<double>
#define max max<double>
#define mod fmod
#define boxcartrain g_model.boxcartrain
#define rednoise g_model.rednoise
#define dwiener g_model.dwiener
#include "model.tmp"
#undef power
#undef solver
#undef min
#undef max
#undef mod
#undef boxcartrain
#undef rednoise
#undef dwiener

#define binaryoutput  //binary reading is much faster in MATLAB


double TExternvar::valueat (double t)
{
   double tstart, tend;
   if (data.ncols == 2)
   {
      tstart = data (0, 0);
      tend = data (data.nrows - 1, 0);
   }
   else
   {
      tstart = 0;
      tend = data.nrows-1;
   }
   if (tofloor)
      t = floor (t);
   if (cycle)
   {
      t = fmod (t, tend - tstart) + tstart;
   }
   if ( ( ( (t < tstart - 1E-10) || (t > tend + 1E-10))  && !cycle) || (data.nrows == 0) || !active)
      return defaultvar;
   if (data.ncols == 2)
   {
      int it = previ;
      if (t < prevt) // save the previous time and it if the next time is lower go back to 0
      {
         if ( (it > 1) && (t > data (it - 2, 0)))
            it--;
         else
            it = 0;
      }
      while ( (it < data.nrows) && (t > data (it, 0)))
         it++;
      int itprev=it -1;
      if (itprev<0)
        itprev=0;
      prevt = data (itprev, 0);
      double nextt = data (it, 0);
      previ = it;
      if (abs(nextt-prevt)>1E-40)
         return  data (itprev, 1) + (data (it, 1)  - data (itprev, 1)) / (nextt - prevt) * (t - prevt);
      else
          return data(it,1);
   }
   else
   {
      int it = floor (t);
      if (it < data.nrows-1 && it >= 0)
         return data (it, 0) + (data (it + 1, 0) - data (it, 0))  * (t - it);
      else if (it<data.nrows && it >= 0)
         return data (it, 0);
      else
         return nan ("");
   }
}


void TModel::setstatevars (int n)
{
#ifdef vectormodel
   v_y_start.resize (n);
   v_y0.resize (n);
   v_dydt.resize (n);
   for (int i = 0; i < n; i++) //v_y0 and v_dydt have no own data
   {
      v_y0[i].data = nullptr;
      v_dydt[i].data = nullptr;
   }
   isvector = true;
#else
   y_start.resize (n, 0);
#endif
}
void TModel::setparameters (int n)
{
#ifdef vectormodel
   v_param.resize (n);
#else
   param.resize (n, 0);
#endif
}
void TModel::setpermanent (int n)
{
#ifdef vectormodel
   v_permanent.resize (n);
#else
   permanent.resize (n, 0);
#endif
}
void TModel::setauxil (int n)
{
#ifdef vectormodel
   v_auxil.resize (n);
#endif
}


double TModel::solver (string par)
//function to mimic the GRIND command solver('step')
{
   if (par == "step")
      return (h);
   else if (par == "ndays")
      return (ndays);
   else
      return nan ("");
}

void TModel::setdiffer (bool d)
//Call this function if the model is a difference equation (changes the possible solvers)
{
   isdiffer = d;
   if (d)
      solverfcn = &TModel::rundiffer;
}

#ifdef binaryoutput
void TModel::readparameters (string filename)
{
   //read parameters from file
   const int signature = 58275485;
   const int c_param =  1;
   const int c_svar =  2;
   const int c_permanent =  3;
   const int c_auxil =  4;
   const int c_step =  5;
   const int c_reltol =  6;
   const int c_abstol =  7;
   const int c_tstart =  8;
   const int c_tend =  9;
   const int c_nout =  10;
   const int c_srand =  11;
   const int c_outf =  12;
   const int c_solver = 13;
   const int c_tout =  14;
   const int c_extern_data =  15;
   const int c_extern_cycle =  16;
   const int c_extern_active = 17;
   const int c_iters =  18;
   const int c_backw =  19;
   const int c_stats =  20;
   const int c_refine =  21;
   const int c_extern_default = 22;
   const int c_ndays = 23;
   const int c_gcycl = 24;
   const int c_nonnegative = 25;
   nonnegative.resize(0);
   ifstream infile;
   infile.open (filename, ios::binary);
   if (infile.is_open())
   {
      int command, signat;
      int lastcomm = 0;
      double d;
      infile.read ( (char*) &signat, sizeof (signat));
      if (signat == signature)
         while (infile.good())
         {
            infile.read ( (char *) (&command), sizeof (command));
            //        cout <<command<<endl;
            int ndx, siz, nx, ny, i;
            switch (command)
            {
            case c_param:
#ifdef vectormodel
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&nx), sizeof (nx));
               infile.read ( (char *) (&ny), sizeof (ny));
               v_param[ndx].setsize (nx, ny);
               infile.read ( (char*) (v_param[ndx].data), nx * ny * sizeof (double));
#else
               infile.read ( (char *) (&siz), sizeof (siz));
               param.resize (siz);
               infile.read ( (char *) param.data(), siz * sizeof (double));
#endif //vectormodel
               break;
            case c_svar:
#ifdef vectormodel
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&nx), sizeof (nx));
               infile.read ( (char *) (&ny), sizeof (ny));
               v_y_start[ndx].setsize (nx, ny);
               infile.read ( (char*) (v_y_start[ndx].data), nx * ny * sizeof (double));
#else
               infile.read ( (char *) (&siz), sizeof (siz));
               y_start.resize (siz);
               infile.read ( (char *) y_start.data(), siz * sizeof (double));

#endif //vectormodel
               break;
            case c_permanent:
#ifdef vectormodel
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&nx), sizeof (nx));
               infile.read ( (char *) (&ny), sizeof (ny));
               v_permanent[ndx].setsize (nx, ny);
               infile.read ( (char*) (v_permanent[ndx].data), nx * ny * sizeof (double));
#else
               infile.read ( (char *) (&siz), sizeof (siz));
               permanent.resize (siz);
               infile.read ( (char *) permanent.data(), siz * sizeof (double));
#endif //vectormodel
               break;
            case c_auxil:
#ifdef vectormodel
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&nx), sizeof (nx));
               infile.read ( (char *) (&ny), sizeof (ny));
               v_auxil[ndx].setsize (nx, ny);
#endif //vectormodel
               break;
            case c_step:
               infile.read ( (char *) (&h), sizeof (h));
               break;
            case c_reltol://else if (comm == "reltol")
               infile.read ( (char *) (&rtol), sizeof (rtol));
               break;
            case c_abstol://else if (comm == "abstol")
               infile.read ( (char *) (&atol), sizeof (atol));
               break;
            case c_tstart://else if (comm == "tstart")
               infile.read ( (char *) (&tstart), sizeof (tstart));
               break;
            case c_ndays://else if (comm == "tend")
               infile.read ( (char *) (&ndays), sizeof (ndays));
               break;
            case c_tend://else if (comm == "tend")
               infile.read ( (char *) (&tend), sizeof (tend));
               break;
            case c_nout://else if (comm == "nout")
               infile.read ( (char *) (&nout), sizeof (nout));
               break;
            case c_srand://else if (comm == "srand")
               infile.read ( (char *) (&d), sizeof (d));
               srand (round (RAND_MAX * d));
               break;
            case c_outf://else if (comm == "outf")
               infile.read ( (char *) (&siz), sizeof (siz));
               outfilename.resize (siz);
               infile.read ( (char*) (outfilename.data()), sizeof (char) *siz);
               //             cout<<outfilename<<endl;
               break;
            case c_solver://else if ( (comm == "solver") && !isdiffer)     //difference equations have one solver
               infile.read ( (char *) (&ndx), sizeof (ndx));
               if (ndx == 0)
               {
                  solverfcn = &TModel::singlerun;
                  nout = 1;
               }
               if (!isdiffer)
                  switch (ndx)
                  {
                  case 1:
                     solverfcn = &TModel::runode45;
                     break;
                  case 2:
                     solverfcn = &TModel::runrk4;
                     break;
                  case 3:
                     solverfcn = &TModel::runeuler;
                     break;
                  }
               break;
            case c_gcycl:
               infile.read ( (char *) (&siz), sizeof (siz));
               boxcar_gcycl.resize (siz);
               infile.read ( (char *) boxcar_gcycl.data(), siz * sizeof (double));
               break;
            case c_nonnegative:
               infile.read ( (char *) (&siz), sizeof (siz));
               nonnegative.resize (siz);
               infile.read ( (char *) nonnegative.data(), siz * sizeof (int));
               break;
            case c_tout:
               infile.read ( (char *) (&siz), sizeof (siz));
               t_out.resize (siz);
               infile.read ( (char *) t_out.data(), siz * sizeof (double));
               nout = siz;
               tstart = t_out[0];
               tend = t_out[siz - 1];
               break;
            case c_extern_data://else if (comm == "extern.data")
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&nx), sizeof (nx));
               infile.read ( (char *) (&ny), sizeof (ny));
               externvars[ndx].data.setsize (nx, ny);
               infile.read ( (char*) (externvars[ndx].data.data), nx * ny * sizeof (double));
               break;
            case c_extern_default://else if (comm == "extern.default")
               infile.read ( (char *) (&ndx), sizeof (ndx));
               //    infile.read ( (char *) (&nx), sizeof (nx));
               //     infile.read ( (char *) (&ny), sizeof (ny));
               //   externvars[ndx].setdim (nx, ny);
               //     infile.read ( (char*) (externvars[ndx].defaultvar.data), nx * ny * sizeof (double));
               infile.read ( (char *) (&d), sizeof (d));
               externvars[ndx].defaultvar = d;
               break;
            case c_extern_cycle://else if (comm == "extern.cycle")
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&i), sizeof (i));
               externvars[ndx].cycle = (i == 1);
               break;
            case c_extern_active://else if (comm == "extern.active")
               infile.read ( (char *) (&ndx), sizeof (ndx));
               infile.read ( (char *) (&i), sizeof (i));
               externvars[ndx].active = (i == 1);
               break;
            case c_iters://else if (comm == "iters")
               infile.read ( (char *) (&iters), sizeof (iters));
               break;
            case c_backw://else if (comm == "backw")
               infile.read ( (char *) (&i), sizeof (i));
               backward = (i == 1);
               break;
            case c_refine://else if (comm == "refine")
               infile.read ( (char *) (&i), sizeof (i));
               refine = (i == 1);
               break;
            case c_stats://else if (comm == "stats")
               infile.read ( (char *) (&i), sizeof (i));
               stats = (i == 1);
               break;
            default:
               cout << "error reading file, last command=" << lastcomm << endl;
               break;
            }
            lastcomm = command;
         }
   }
   else
   {
      cout << "Unable to open file " << filename;
   }
   infile.close();

}
#else
void TModel::readparameters (string filename)
{
   //read parameters from file
   ifstream infile;
   infile.open (filename);
   if (infile.is_open())
   {
      string line, comm;
      while (infile.good())
      {
         getline (infile, line);
         size_t p = line.find (" ");
         comm = line.substr (0, p);
         unsigned int ndx, nx, ny;
         double valu;
         if (comm == "para")
         {
#ifdef vectormodel
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            v_param[ndx].str2data (line, p1);
#else
            parseline (line, p, ndx, nx, ny, valu);
            if (ndx < param.size())
               param[ndx] = valu;
#endif //vectormodel
         }
         else if (comm == "svar")
         {
#ifdef vectormodel
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            v_y_start[ndx].str2data (line, p1);
            v_y0[ndx].setsize (v_y_start[ndx].nrows, v_y_start[ndx].ncols);
            v_dydt[ndx].setsize (v_y_start[ndx].nrows, v_y_start[ndx].ncols);
#else
            parseline (line, p, ndx, nx, ny, valu);
            if (ndx < y_start.size())
               y_start[ndx] = valu;
#endif //vectormodel
         }
         else if (comm == "perm")
         {
#ifdef vectormodel
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            v_permanent[ndx].str2data (line, p1);
#else
            parseline (line, p, ndx, nx, ny, valu);
            if (ndx < permanent.size())
               permanent[ndx] = valu;
#endif //vectormodel
         }
         else if (comm == "auxil")
         {
#ifdef vectormodel
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            v_auxil[ndx].str2data (line, p1);
#endif //vectormodel
         }
         else if (comm == "step")
         {
            h = atof (line.substr (p + 1).c_str());
         }
         else if (comm == "reltol")
         {
            rtol = atof (line.substr (p + 1).c_str());
         }
         else if (comm == "abstol")
         {
            atol = atof (line.substr (p + 1).c_str());
         }
         else if (comm == "tstart")
         {
            tstart = atof (line.substr (p + 1).c_str());
         }
         else if (comm == "tend")
         {
            tend = atof (line.substr (p + 1).c_str());
         }
         else if (comm == "nout")
         {
            nout = atof (line.substr (p + 1).c_str()) + 1;
         }
         else if (comm == "srand")
         {
            double r = atof (line.substr (p + 1).c_str());
            srand (round (RAND_MAX * r));
         }
         else if (comm == "outf")
         {
            size_t p1 = line.find ('"', p + 1);
            size_t p2 = line.find ('"', p1 + 1);
            outfilename = line.substr (p1 + 1, p2 - p1 - 1);
         }
         else if ( (comm == "solver") && !isdiffer)     //difference equations have one solver
         {
            string solvername = line.substr (p + 1);
            if (solvername == "euler") solverfcn = &TModel::runeuler;
            if (solvername == "rk4") solverfcn = &TModel::runrk4;
         }
         else if (comm == "tout")
         {
            unsigned int n = 0;
            size_t p1;
            bool endfound = false;
            while (!endfound)
            {
               if (t_out.size() <= n)
                  t_out.resize (n + 1000);
               p1 = line.find (' ', p + 1);
               if (p1 != string::npos)
               {
                  t_out[n] = atof (line.substr (p + 1, p1 - p).c_str());
                  n++;
                  p = p1;
               }
               else
               {
                  t_out[n] = atof (line.substr (p + 1).c_str());
                  n++;
                  endfound = true;
               }
            }
            if (t_out.size() > n)
               t_out.resize (n);
         }
         else if (comm == "extern.data")
         {
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            externvars[ndx].data.str2data (line, p1);
         }
         else if (comm == "extern.cycle")
         {
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            p = line.find (' ', p1 + 1) ;
            externvars[ndx].cycle = atoi (line.substr (p + 1).c_str());
         }
         else if (comm == "extern.active")
         {
            int p1 = line.find (' ', p + 1);
            int ndx = atoi (line.substr (p + 1, p1 - p).c_str());
            p = line.find (' ', p1 + 1) ;
            externvars[ndx].active = atoi (line.substr (p + 1).c_str());
         }
         else if (comm == "iters")
         {
            iters = atoi (line.substr (p + 1).c_str());
         }
         else if (comm == "refine")
         {
            refine = atoi (line.substr (p + 1).c_str());
         }
         else if (comm == "backw")
         {
            backward = atoi (line.substr (p + 1).c_str());
         }
         else if (comm == "stats")
         {
            stats = atoi (line.substr (p + 1).c_str());
         }
      }
   }
   else
   {
      cout << "Unable to open file " << filename;
   }
   infile.close();

}
#endif //binary
double TModel::rednoise (double t, double T0, double lambda, double beta, int iset /*= 1*/, double deltat /*= 1*/)
{
   if (iset >= int (rednoisesets.size()))
      rednoisesets.resize (iset);
   iset = iset - 1;
   if (rednoisesets[iset].data.numel == 0)
   {
      if ((deltat == 1) && (tstart==0))
      {
         rednoisesets[iset].tofloor = true;
         int tspan=1;
         if (abs(tend-tstart)>1E-40)
            tspan=round (tend-tstart);
         tspan=tspan+5;
         rednoisesets[iset].data.setsize(tspan,1);
         rednoisesets[iset].data (0,0) = T0;
         for (int i = 1; i < tspan; i++)
            rednoisesets[iset].data (i, 0) = (1 - 1 / lambda) * (rednoisesets[iset].data (i - 1, 0) - T0) + T0 + beta * randn();
      }
      else
      {
         rednoisesets[iset].tofloor = true;
         int n = round ( (tend - tstart + 1) / deltat);
         rednoisesets[iset].data.setsize (n, 2);
         rednoisesets[iset].data (0, 0) = tstart;
         rednoisesets[iset].data (0, 1) = T0;
         for (int i = 1; i < n; i++)
         {
            rednoisesets[iset].data (i, 0) = tstart + i * deltat;
            rednoisesets[iset].data (i, 1) = (1 - 1 / lambda) * (rednoisesets[iset].data (i - 1, 1) - T0) + T0 + beta * randn();
         }
      }
   }
   return rednoisesets[iset].valueat (t);
}

double TModel::dwiener (double gY, double dgY_dY)
{
   double dW = randn() * sqrt (h);
   //Milstein scheme (reduces to Euler-Maruyama if d(gY)/dY=0)
   return (gY * dW + gY * dgY_dY * (dW * dW - h)) / h; //% divide by h as the Euler routine multiplies with h
}

void TModel::parseline (string & line, size_t & p, unsigned int & ndx, unsigned int & nx, unsigned int & ny, double & valu)
{
   size_t p1 = line.find (" ", p + 1);
   ndx = atoi (line.substr (p + 1, p1 - p).c_str());
   size_t p2 = line.find (" ", p1 + 1);
   nx = atoi (line.substr (p1 + 1, p1 - p2).c_str());
   size_t p3 = line.find (" ", p2 + 1);
   ny = atoi (line.substr (p2 + 1, p3 - p2).c_str());
   valu = atof (line.substr (p3 + 1).c_str());
}
void TModel::interpstep (double t1, double t2, vector<double>* X1, vector<double>* X2, double tnew, vector<double> permanent)
{
   int neq = X1->size();
   vector<double> yinterp (neq);
   for (int i = 0; i < neq; i++)
   {
      yinterp[i] = (*X1) [i] + ( (*X2) [i] - (*X1) [i]) / (t2 - t1) * (tnew - t1);
   }
   writeoutput (yinterp.data(), neq);
   writeoutput (permanent.data(), permanent.size());
#ifndef binaryoutput
   outfile << endl;
#endif // binaryoutput
}

void TModel::runsolver ()
//Function to run the current solver
{
   (this->*solverfcn) (); //ugly way to call a pointer to a method
}

void TModel::openoutfile ()
{
#ifdef binaryoutput
   outfile.open (outfilename, ios::binary);
   int const signature = 1234567890;
   outfile.write ( (char *) &signature, sizeof signature);
   int ncol = y_start.size() + permanent.size();
   int nrow = t_out.size();
   if (includet)
   {
      nrow = -1;
      ncol++;
   }
   outfile.write ( (char *) &ncol, sizeof ncol);
   outfile.write ( (char *) &nrow, sizeof nrow);
#else
   outfile.open (outfilename);
   outfile.precision (12);
#endif //binaryoutput
   if (includet && (solverfcn != &TModel::singlerun))
   {
      writeoutput (&tstart);
      writeoutput (y_start.data(), y_start.size()) ;
      writeoutput (permanent.data(), permanent.size()) ;
   }
}
const double Eps = 1E-15;

void TModel::singlerun ()
/*Single run, no integration
*/
{

   openoutfile ();
   tend=tstart;
   runodefile (tstart, y0.data(), dydt.data());
   writeoutput (dydt.data(), dydt.size());
   //outfile.write ( (char *) dydt.data(), sizeof(double)*dydt.size());
   outfile.close();
}
void TModel::writeoutput (double *y, int n)
{
#ifdef binaryoutput
   outfile.write ( (char *) y, sizeof (double) *n);
#else
   for (int i = 0; i < n; i++)
      outfile << y[i] << " ";
#endif // binaryoutput
}

void TModel::runeuler ()
/*Euler integration
Simplest integration, usefull for stochastic models
not very precise


% Test that tspan is internally consistent.
tspan = tspan(:);
ntspan = length(tspan);
if ntspan == 1
   t0 = 0;
 %  next = 1;
else
   t0 = tspan(1);
%   next = 2;
end
tfinal = tspan(ntspan);
if t0 == tfinal
   error('GRIND:euler:tpan','The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t0);
if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan-1)) <= 0)
   error('GRIND:euler:tspan','The entries in tspan must strictly increase or decrease.');
end
t = t0;
y = y0(:);
neq = length(y);
%adapt delta if there should be given more output
%delta=min(delta,median(diff(tspan)));

% Set the output flag.

%outflag = ntspan > 2 && abs(delta*(ntspan-1)-abs(tspan(end)-tspan(1)))>1E-10;     % output only at tspan points

% Allocate memory if we're generating output.
delta = delta * tdir;

if ntspan > 2                         % output only at tspan points
    tout = tspan(:);
    yout = zeros(ntspan, neq);
else
    tout = transpose(linspace(t0,tfinal,round(tfinal/delta)));
    yout = zeros(size(tout, 1),neq);
end
nout=size(tout,1);
i_out = 1;
yout(i_out, :) = y.';
nextt=tout(i_out+1);
%
%MAIN LOOP
%evaluate the odefunction for the next time steps
%fold=[];
h=delta;
running=1;
while running
   f=feval(odefunct, t, y);
   %t is the time of the first point
%  simple Euler
   tnew = min(nextt,t + delta); %adapt the time step if approaching an output point
   h=tnew-t;
   ynew = y + f .*  h;
   if anyNonNegative
       ndx = nonNegative( ynew(nonNegative) < 0 );
       ynew(ndx) = max(ynew(ndx),0);
   end;

   if tnew==nextt
      i_out = i_out + 1;
      yout(i_out, :) = ynew.';
      if i_out+1<=nout
         nextt=tout(i_out+1);
      else
         running=false;
      end
      if haveOutputFcn
        if feval(outputFcn,tnew,ynew','')
          running = false;
        end
      end
   end;
   y = ynew;
   t = tnew;
end;
*/
{
   openoutfile();
   double forwrd_h = h; // if runnning backwards h=-h in simulation
   if (backward) forwrd_h = -forwrd_h;
   unsigned int i;
   unsigned int tt = 0;
   vector<double> Yend (dydt.size(), 0);
   vector<double> * ystart = &y0;
   vector<double> * yend = &Yend;
   vector<double> * yhelp = NULL;
   unsigned int nsteps = round((tend-tstart)/h)+1;
   double t=tstart;
   double nextt;
   for (int j = 0; j< nsteps; j++)
   {
      runodefile (t, ystart->data(), dydt.data());
      for (i = 0; i < ystart->size(); i++)
      {
         (*yend) [i] = (*ystart) [i] + forwrd_h * dydt[i];
      }
      if (nonnegative.size() > 0)
      {
          for (i = 0; i < nonnegative.size(); i++)
          {
             if ((*yend) [i] < 0) (*yend) [i]=0;

          }
      }
     nextt = t + h;
    if (nextt > tend)
           nextt = tend;
     if (includet)
      {
         writeoutput (&nextt);
         writeoutput (yend->data(), yend->size());
         writeoutput (permanent.data(), permanent.size());
      }
      else
        while ((tt < t_out.size()) && (sign(forwrd_h) * (nextt - t_out[tt])) >= -h*0.1)
       //  while ( ((t_out[tt] - 1E-10- t  < h))
         {
            interpstep (t, nextt, ystart, yend, t_out[tt], permanent); //interpolate output
            tt++;
         }
      yhelp = yend;
      yend = ystart;
      ystart = yhelp;
      t+=h;
   }
   outfile.close();
   if (stats)
   {
      cout << nsteps << " steps\n";
      cout << nsteps << " function evaluations\n";
   }
}



void TModel::runrk4 ()
/*Standard Runga Kutta 4 integration, fixed step size
much more precise than euler
k1 = h*f(tn,yn);
k2 = h*f(tn+0.5*h ,yn+ 0.5*k1);
k3 = h*f(tn+0.5*h ,yn+ 0.5*k2);
k4 = h*f(tn+h , yn + k3);
yn+1= yn+1/6*(*k1+2*k2+2*k3+k4)
tn+1=tn+h;
*/
{
   openoutfile ();
   double forwrd_h =  h; // if runnning backwards forwrd=-1, for efficiency this is combined
   if (backward) forwrd_h = -forwrd_h;
   unsigned int i;
   unsigned int tt = 0;
   unsigned int nsteps = 0;
   vector<double> Yend (dydt.size(), 0);
   vector<double> Yh (dydt.size(), 0);
   vector<double> k1 (dydt.size(), 0);
   vector<double> k2 (dydt.size(), 0);
   vector<double> k3 (dydt.size(), 0);
   vector<double> k4 (dydt.size(), 0);
   vector<double> * ystart = &y0;
   vector<double> * yend = &Yend;
   vector<double> * yhelp = NULL;
   for (double t = tstart; t < tend; t += h)
   {
      runodefile (t + h, Yh.data(), k1.data());
      for (i = 0; i < k1.size(); i++)
      {
         k1[i] = forwrd_h * k1[i];
         Yh[i] = (*ystart) [i] + 0.5 * k1[i];
      }
      runodefile (t + 0.5 * h, Yh.data(), k2.data());
      for (i = 0; i < k2.size(); i++)
      {
         k2[i] = forwrd_h * k2[i];
         Yh[i] = (*ystart) [i] + 0.5 * k2[i];
      }
      runodefile (t + 0.5 * h, Yh.data(), k3.data());
      for (i = 0; i < k2.size(); i++)
      {
         k3[i] = forwrd_h * k3[i];
         Yh[i] = (*ystart) [i] + k3[i];
      }
      runodefile (t + h, Yh.data(), k4.data());
      for (i = 0; i < k4.size(); i++)
      {
         (*yend) [i] = (*ystart) [i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + forwrd_h * k4[i]) / 6.0;
      }
      if (nonnegative.size() > 0)
      {
          for (i = 0; i < nonnegative.size(); i++)
          {
             if ((*yend) [i] < 0) (*yend) [i]=0;
          }
      }
      nsteps++;
      if (includet)
      {
         double nextt = t + h;
         if (nextt > tend)
            nextt = tend;
         writeoutput (&nextt);
         writeoutput (yend->data(), yend->size());
         writeoutput (permanent.data(), permanent.size());
      }
      else
         while ( (tt < t_out.size()) && (t_out[tt] - t - Eps < h))
         {
            interpstep (t, t + h, ystart, yend, t_out[tt], permanent);
            tt++;
         }
      yhelp = yend; //Exchange pointers to yend and ystart
      yend = ystart;
      ystart = yhelp;
   }
   outfile.close();
   if (stats)
   {
      cout << nsteps << " steps\n";
      cout << nsteps * 4 << " function evaluations\n";
   }
}

void TModel::ntrp45 (double tinterp, double t, double * y, unsigned int neq, double h, double * f1,  double * f3, double * f4, double * f5, double * f6, double * f7)
//NTRP45  Interpolation helper function for ODE45.
//   YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) uses data computed in ODE45
//   to approximate the solution at time TINTERP.  TINTERP may be a scalar
//   or a row vector.
//   The arguments TNEW and YNEW do not affect the computations. They are
//   required for consistency of syntax with other interpolation functions.
//   Any values entered for TNEW and YNEW are ignored.
{

   double const BI[7][4] =
   {
      {1.0 , -183.0 / 64.0     , 37.0 / 12.0   ,    -145.0 / 128.0},
      {0.0  ,        0.0     ,      0.0      ,      0.0},
      {0  ,     1500.0 / 371.0 ,   -1000.0 / 159.0 ,   1000.0 / 371.0},
      {0  ,     -125.0 / 32.0   ,    125.0 / 12.0  ,   -375.0 / 64.0},
      {0  ,     9477.0 / 3392.0 ,  -729.0 / 106.0  ,  25515.0 / 6784.0},
      {0   ,     -11.0 / 7.0   ,     11.0 / 3.0    ,    -55.0 / 28.0},
      {0  ,       3.0 / 2.0   ,      -4.0      ,      5.0 / 2.0}
   };
   double hBI[7][4];
   for (int i = 0; i < 7; i++)
      for (int j = 0; j < 4; j++)
         hBI[i][j] = h * BI[i][j];
   double s = (tinterp - t) / abs (h);
   vector<double> yinterp (neq);
   for (unsigned int i = 0; i < neq; i++)
   {
      double s1 = s;
      yinterp[i] = y[i] + f1[i] * hBI[0][0] * s1;
      for (int j = 1; j < 4; j++)
      {
         s1 = s1 * s;
         yinterp[i] = yinterp[i] + (f1[i] * hBI[0][j] + f3[i] * hBI[2][j] + f4[i] * hBI[3][j] + f5[i] * hBI[4][j] + f6[i] * hBI[5][j] + f7[i] * hBI[6][j]) * s1;

      }
   }
   writeoutput (yinterp.data(), neq);
   writeoutput (permanent.data(), permanent.size());
#ifndef binaryoutput
   outfile << endl;
#endif // binaryoutput
}




void TModel::runode45 ()
/*ODE45 based on the implementation in MATLAB
%   ODE45 is an implementation of the explicit Runge-Kutta (4,5) pair of
%   Dormand and Prince called variously RK5(4)7FM, DOPRI5, DP(4,5) and DP54.
%   It uses a "free" interpolant of order 4 communicated privately by
%   Dormand and Prince.  Local extrapolation is done.

%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

*/
{
   openoutfile ();
   unsigned int tt = 0;
   unsigned int neq = dydt.size();
   unsigned int nperm = permanent.size();
   vector<double> y1 (neq, 0);
   vector<double> yH (neq, 0);
   vector<double> ynew1 (neq, 0);
   vector<double> * y = &y1;
   vector<double> * ynew = &ynew1;

   vector<double> f01 (neq, 0);
   vector<double> f02 (neq, 0);
   vector<double> f03 (neq, 0);
   vector<double> f04 (neq, 0);
   vector<double> f05 (neq, 0);
   vector<double> f06 (neq, 0);
   vector<double> f07 (neq, 0);
   vector<double> oldpermanent (nperm, 0);

   vector<double> * f1 = &f01;
   vector<double> * f2 = &f02;
   vector<double> * f3 = &f03;
   vector<double> * f4 = &f04;
   vector<double> * f5 = &f05;
   vector<double> * f6 = &f06;
   vector<double> * f7 = &f07;

   vector<double> * fhelp = nullptr;
   bool done = false;
   bool hasnonnegative = nonnegative.size() > 0;
   const double powA = 1.0 / 5.0;

   const double A[6] = {1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};
   const double B[7][6] =
   {
      {1.0 / 5.0, 3.0 / 40.0, 44.0 / 45.0, 19372.0 / 6561.0, 9017.0 / 3168.0, 35.0 / 384.0},
      {0.0, 9.0 / 40.0, -56.0 / 15.0, -25360.0 / 2187.0, -355.0 / 33.0, 0.0},
      {0.0, 0.0, 32.0 / 9.0, 64448.0 / 6561.0, 46732.0 / 5247.0, 500.0 / 1113.0},
      {0.0, 0.0, 0.0, -212.0 / 729.0, 49.0 / 176.0, 125.0 / 192.0},
      {0.0, 0.0, 0.0, 0.0, -5103.0 / 18656.0, -2187.0 / 6784.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 11.0 / 84.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   };
   double hB[7][6];
   const double E1 = 71.0 / 57600.0;
   const double E3 = -71.0 / 16695.0;
   const double E4 = 71.0 / 1920.0;
   const double E5 = -17253.0 / 339200.0;
   const double E6 = 22.0 / 525.0;
   const double E7 = -1.0 / 40.0;
   double hmin = 16.0 * eps (tstart);

   //% Compute an initial step size h using y'(t).

   int tdir = sign (tend - tstart);
   int hbackw = 1;
   if (backward) hbackw = -1;
   runodefile (tstart, y0.data(), f1->data());
   double htspan;
   if (t_out.size() > 1)
      htspan = abs (t_out[1] - tstart);
   else
      htspan = abs (tend - tstart);
   //stats
   unsigned int nsteps = 0;
   unsigned int nfailed = 0;
   unsigned int nfevals = 1;

   double hmax = 0.1 * (tend - tstart); //can be made an option
   double absh = min (hmax, htspan);
   double threshold = atol / rtol;
   double maxf = 0;
   for (unsigned int i = 0; i < neq; i++)
   {
      (*y) [i] = y0[i];
      double f0 = abs ( (*f1) [i] / (max (y0[i], threshold)));
      if (f0 > maxf) maxf = f0;
   }
   double rh = maxf / (0.8 * pow (rtol, powA));
   if (absh * rh > 1)
      absh = 1.0 / rh;
   absh = max (absh, hmin);
// % THE MAIN LOOP
   double t = tstart;
   double tnew;
   double err = 0;
   while (!done)
   {
      //% By default, hmin is a small number such that t+hmin is only slightly
      //% different than t.  It might be 0 if t is 0.
      hmin = 16 * eps (t);
      absh = min (hmax, max (hmin, absh)); //  % couldn't limit absh until new hmin
      h = tdir * absh;

      //% Stretch the step if within 10% of tfinal-t.
      if (1.1 * absh >= abs (tend - t))
      {
         h = tend - t;
         absh = abs (h);
         done = true;
      }

      //% LOOP FOR ADVANCING ONE STEP.
      bool nofailed = true;  //% no failed attempts
      for (unsigned int i = 0; i < nperm; i++)
         oldpermanent[i] = permanent[i];
      while (true)
//advance one step
      {
//       hB=h.*B
         for (int i = 0; i < 7; i++)
            for (int j = 0; j < 6; j++)
               hB[i][j] = hbackw * h * B[i][j];
//       f2 = odefile(t+h*A[1],y+f1*hB(:,1))
         for (unsigned int i = 0; i < neq; i++)
            yH[i] = (*y) [i] + hB[0][0] * (*f1) [i];
         runodefile (t + h * A[0], yH.data(), f2->data());
         if (hasnonnegative) {for (unsigned int i = 0; i < nonnegative.size(); i++){if (yH[i] < 0) yH[i]=0;}}
//       f3 = odefile(t+h*A[2],y+f2*hB(:,2));
         for (unsigned int i = 0; i < neq; i++)
            yH[i] = (*y) [i] + hB[0][1] * (*f1) [i] + hB[1][1] * (*f2) [i];
         runodefile (t + h * A[1], yH.data(), f3->data());
        if (hasnonnegative) {for (unsigned int i = 0; i < nonnegative.size(); i++){if (yH[i] < 0) yH[i]=0;}}
         //       f4 = odefile(t+h*A[3],y+f3*hB(:,3));
         for (unsigned int i = 0; i < neq; i++)
            yH[i] = (*y) [i] + hB[0][2] * (*f1) [i] + hB[1][2] * (*f2) [i] + hB[2][2] * (*f3) [i];
         runodefile (t + h * A[2], yH.data(), f4->data());
        if (hasnonnegative) {for (unsigned int i = 0; i < nonnegative.size(); i++){if (yH[i] < 0) yH[i]=0;}}
//       f5 = odefile(t+h*A[4],y+f4*hB(:,4));
         for (unsigned int i = 0; i < neq; i++)
            yH[i] = (*y) [i] + hB[0][3] * (*f1) [i] + hB[1][3] * (*f2) [i] + hB[2][3] * (*f3) [i] + hB[3][3] * (*f4) [i];
         runodefile (t + h * A[3], yH.data(), f5->data());
        if (hasnonnegative) {for (unsigned int i = 0; i < nonnegative.size(); i++){if (yH[i] < 0) yH[i]=0;}}
//       f6 = odefile(t+h*A[5],y+f5*hB(:,5));
         for (unsigned int i = 0; i < neq; i++)
            yH[i] = (*y) [i] + hB[0][4] * (*f1) [i] + hB[1][4] * (*f2) [i] + hB[2][4] * (*f3) [i] + hB[3][4] * (*f4) [i] + hB[4][4] * (*f5) [i];
         runodefile (t + h * A[4], yH.data(), f6->data());
        if (hasnonnegative) {for (unsigned int i = 0; i < nonnegative.size(); i++){if (yH[i] < 0) yH[i]=0;}}

         tnew = t + h * A[5];
         if (done)
            tnew = tend + eps (1); // Hit end point exactly.
         h = tnew - t;      // Purify h.
         for (unsigned int i = 0; i < neq; i++)
            (*ynew) [i] = (*y) [i] + hB[0][5] * (*f1) [i] + hB[1][5] * (*f2) [i] + hB[2][5] * (*f3) [i] + hB[3][5] * (*f4) [i] + hB[4][5] * (*f5) [i] + hB[5][5] * (*f6) [i];
         runodefile (tnew, ynew->data(), f7->data());
         if (hasnonnegative) {for (unsigned int i = 0; i < nonnegative.size(); i++){if ((*ynew) [i] < 0) (*ynew) [i]=0;}}
         nfevals += 6;
         //calculate the error
         err = 0;
         for (unsigned int i = 0; i < neq; i++)
         {
            double herr = ( (*f1) [i] * E1 + (*f3) [i] * E3 + (*f4) [i] * E4 + (*f5) [i] * E5 + (*f6) [i] * E6 + (*f7) [i] * E7) / max (max (abs ( (*y) [i]), abs ( (*ynew) [i])), threshold);
            if (abs (herr) > err) err = abs (herr);
         }
         err = absh * err;
         if (err > rtol)
         {
            nfailed++;
            for (unsigned int i = 0; i < nperm; i++)
               permanent[i] = oldpermanent[i];
            if (absh <= hmin)
            {
               cout << "code45:IntegrationTolNotMet, t=" << t << "  hmin=" << hmin ;
               return;
            }
            if (nofailed)
            {
               nofailed = false;
               absh = max (hmin, absh * max (0.1, 0.8 * pow ( (rtol / err) , powA)));
            }
            else
               absh = max (hmin, 0.5 * absh);
            h = tdir * absh;
            done = false;
         }
         else
            break;
      }
      nsteps++;
      if (nofailed)
      {
         // % Note that absh may shrink by 0.8, and that err may be 0.
         double temp = 1.25 * pow ( (err / rtol), powA);
         if (temp > 0.2)
            absh = absh / temp;
         else
            absh = 5.0 * absh;
      }
      if (includet)
      {
         if (refine)
         {
            double t1 = t + 0.25 * h;
            writeoutput (&t1);
            ntrp45 (t1, t, y->data(), neq, hbackw * h, f1->data(), f3->data(), f4->data(), f5->data(), f6->data(), f7->data());
            t1 += 0.25 * h;
            writeoutput (&t1);
            ntrp45 (t1, t, y->data(), neq, hbackw * h, f1->data(), f3->data(), f4->data(), f5->data(), f6->data(), f7->data());
            t1 += 0.25 * h;
            writeoutput (&t1);
            ntrp45 (t1, t, y->data(), neq, hbackw * h, f1->data(), f3->data(), f4->data(), f5->data(), f6->data(), f7->data());
         }
         writeoutput (&tnew);
         writeoutput (ynew->data(), neq);
         writeoutput (permanent.data(), permanent.size());
      }
      else
         while ( (tt < t_out.size()) && (t_out[tt] - t - Eps < h))
         {
            if (tnew == t_out[tt])
            {
               writeoutput (ynew->data(), neq);
               writeoutput (permanent.data(), permanent.size());
            }
            else
               ntrp45 (t_out[tt], t, y->data(), neq, hbackw * h, f1->data(), f3->data(), f4->data(), f5->data(), f6->data(), f7->data());
            tt++;
         }
      //% Advance the integration one step.
      t = tnew;
      fhelp = y; //exchange pointers y and ynew
      y = ynew;
      ynew = fhelp;

      fhelp = f1; //exchange pointers f1 and f7 as f7 is the new f1
      f1 = f7;
      f7 = fhelp;
   }
   outfile.close();
   if (stats)
   {
      cout << nsteps << " successful steps\n";
      cout << nfailed << " failed attempts\n";
      cout << nfevals << " function evaluations\n";
   }
}

void TModel::rundiffer ()
/*runnning a difference equation
*/
{
   openoutfile ();
   if (backward)
      throw ("Backward is not implemented for difference equations");
   //BACKWARD NOT YET IMPLEMENTED, NEEDS OPTIMIZATION ITERATION TO FIND THE PREVIOUS VALUE
   unsigned int tt = 0;
   vector<double> * ystart = &y0;
   vector<double> * yend = &dydt;
   vector<double> * yhelp = nullptr;
   for (double t = tstart; t < tend + 1.0; t += 1.0)
   {
      for (int i = 1; i <= iters; i++)
      {
         runodefile (t, ystart->data(), yend->data());
         if (nonnegative.size() > 0)
         {
           for (i = 0; i < nonnegative.size(); i++)
           {
              if ((*yend) [i] < 0) (*yend) [i]=0;
           }
         }
         yhelp = yend;
         yend = ystart;
         ystart = yhelp;
      }
      if (includet)
      {
         double nextt = t + 1;
         if (nextt > tend)
            nextt = tend;
         writeoutput (&nextt);
         writeoutput (yend->data(), yend->size());
         writeoutput (permanent.data(), permanent.size());
      }
      else
         while ( (tt < t_out.size()) && (t_out[tt] - Eps - t  < 1))
         {
            interpstep (t, t + 1, ystart, yend, t_out[tt], permanent); //interpolate output
            tt++;
         }
   }
   outfile.close();
}


//we include the odefile function here, so other methods cannot use it directly.

void TModel::runodefile (double t, double * _ystart,  double * _dydt)
//Wrapper around odefile for later use (to initialize Matrices etc).
{
#ifdef vectormodel
   //set the data pointers of the vectors to _ystart and _dydt
   int dim = 0;
   for (unsigned int i = 0; i < v_y0.size(); i++)
   {
      v_y0[i].data = _ystart + dim;
      v_dydt[i].data = _dydt + dim;
      dim += v_y0[i].numel;
   }
#endif // vectormodel
   for (unsigned int i = 0; i < externvars.size(); i++) //update the external variables
      externvalues[i] = externvars[i].valueat (t);
   odefile (t, _ystart, _dydt, *this);
}

TModel::TModel (string filename)
//Constructor of TModel
{
   boxcar_name.resize(0);
   initialize (*this);
   externvalues.resize (externvars.size(), 0);
   readparameters (filename);
   includet = nout < 0;
   if ( (t_out.size() == 0) && !includet) //tout values are not given
   {
      t_out.resize (nout);
      for (int i = 0; i < nout; i++)
         t_out[i] = tstart + i * (tend - tstart) / (nout - 1);
   }
   else if (!includet)    //tout values are given
   {
      nout = t_out.size();
      tstart = t_out[0];
      tend = t_out[t_out.size() - 1];
   }
#ifdef vectormodel
   unsigned int totdim = 0;
   for (unsigned int i = 0; i < v_y_start.size(); i++)
   {
      totdim += v_y_start[i].numel;
      v_y0[i].setsize (v_y_start[i].nrows, v_y_start[i].ncols);
      v_dydt[i].setsize (v_y_start[i].nrows, v_y_start[i].ncols);
   }
   y_start.resize (totdim);
   unsigned int dim = 0;
   for (unsigned int i = 0; i < v_y_start.size(); i++)
   {
      for (int j = 0; j < v_y_start[i].numel; j++)
         y_start[j + dim] = v_y_start[i] (j);
      v_y_start[i].data = y_start.data() + dim;
      dim += v_y_start[i].numel;
   }
   totdim = 0;
   for (unsigned int i = 0; i < v_permanent.size(); i++)
   {
      totdim += v_permanent[i].numel;
   }
   permanent.resize (totdim);
   dim = 0;
   for (unsigned int i = 0; i < v_permanent.size(); i++)
   {
      for (int j = 0; j < v_permanent[i].numel; j++)
         permanent[j + dim] = v_permanent[i] (j);
      v_permanent[i].data = permanent.data() + dim;
      dim += v_permanent[i].numel;
   }
#endif //vectormodel
   dydt.resize (y_start.size(), 0);
   y0.resize (y_start.size(), 0);
   y0.assign (y_start.begin(), y_start.end());
   boxcar_gcycl.resize(boxcar_name.size(),0);
}

boxcar_result TModel::boxcartrain (int boxnr, TMatrix &A, double devrate, double cv_devrate)
{
   double deltat = h;
   int N = A.nrows;
   boxcar_result res;
   double gamma = 1.0 / N; //%%%% origial code of Goudriaan!!
   if ((solverfcn != &TModel::runeuler) && (solverfcn != &TModel::singlerun))
      cout << "Boxcartrain needs Euler integration";
   //flow calculation

   if (isnan (cv_devrate) || (cv_devrate < 0))
   {
      //
      //fixed boxcar train
      //
      //flows from each boxcar
      res.flow.setsize (A.nrows, A.ncols);
      for (int i = 0; i < N; i++)
      {
         res.flow (i) = A (i) / gamma * devrate;
      }
      // %shift the flows to get the inflow to the next boxcar and keep outflow
      res.outflow = res.flow (N-1);
      res.flow (0) = -res.flow (0);
      for (int i = 1; i <N; i++)
      {
         res.flow (i) = A (i-1) / gamma * devrate - res.flow (i);
      }
    }
   else if (cv_devrate == 0)
   {

      //
      //escalator boxcar train
      //

      res.flow.setsize (A.nrows, A.ncols);
      for (int i = 0; i < N; i++)
        res.flow(i)=0;
      res.outflow = 0;
      if (isnan (boxcar_gcycl[boxnr]))
         boxcar_gcycl[boxnr] = 0.5 * gamma; // !! boxcar.gcycl should be reset to nan
      else
         boxcar_gcycl[boxnr] = boxcar_gcycl[boxnr] + devrate * deltat;

      if (boxcar_gcycl[boxnr] >= gamma)
      {
         //flow is the difference between next and current
         res.outflow = res.outflow + A (N-1) / deltat; //all biomass is moved in one moment to the next stage
         res.flow (0) = -A (0)/ deltat;
         for (int i = 1; i < N; i++)
         {
            res.flow (i) = (A (i-1) - A (i)) / deltat; //divide by deltat such that the whole A will be shifted
         }
         boxcar_gcycl[boxnr] = boxcar_gcycl[boxnr] - gamma;
      }
   }
   else
   {

      //
      //fractional boxcar train
      //
      res.flow.setsize (A.nrows, A.ncols);
      for (int i = 0; i < N; i++)
        res.flow(i)=0;
      res.outflow = 0;

      if (isnan (boxcar_gcycl[boxnr]))
         boxcar_gcycl[boxnr] = 0.5 * gamma; // !! g_grind.boxcar.gcycl should be reset to NaN
      else
         boxcar_gcycl[boxnr] = boxcar_gcycl[boxnr] + devrate * deltat;
      //
      //Fraction F to be moved (Goudriaan)
      double F = 1 - (N * pow (cv_devrate, 2));
      if (F < 0)
         throw std::string ("Boxcartrain: coefficent of variation too big or number of stages too big");
//     sprintf("Boxcartrain: coefficent of variation (%g) too big or number of stages (%d) too big\nboxcar %s CV should at least be < %g (<sqrt(1/N))",cv_devrate,N,g_grind.boxcar.names{boxnr},sqrt(1/N));

      if (deltat > F * gamma / (devrate + 1e-30))
         throw std::string ("Boxcartrain: integration step too large");
//    error('GRIND:boxcartrain:IntSteptoobig','Boxcartrain: Integration step too large: deltat should be < %0.5g\nAlternatively decrease the number of classes of %s to %d\ncv = %g devr = %g',...
      //       F * gamma / (devrate + 1e-30),g_grind.boxcar.names{boxnr},floor(1/(deltat*devrate+cv_devrate^2)),cv_devrate,devrate);

      if (gamma - boxcar_gcycl[boxnr]  > 1E-30)
      {

         res.outflow = min (A (N-1) / deltat, A (N-1) / (gamma - boxcar_gcycl[boxnr]) * devrate);
         res.flow (N-1) = -res.outflow;
      }
      else
         res.outflow = 0;
      if (boxcar_gcycl[boxnr] >= F * gamma-1E-10)
      {
         //flow is the difference between next and current
         res.flow (0) =  -A (0)* F / deltat;
//         cout<<res.flow(0)<<endl;
         for (int i = 1; i < N; i++)
         {
            res.flow (i) = (A (i-1) - A (i)) * F / deltat;
//            cout<<res.flow(i)<<endl;
         }
         //    res.flow = ([zeros(1,N2); A(1:N - 1,:)]-A) .* F ./ deltat; %divide by deltat such that the whole A will be shifted
         boxcar_gcycl[boxnr] = boxcar_gcycl[boxnr] - F * gamma;
      }
   }
//update boxcar train (is done directly in the state variable!!)
   for (int i = 0; i < N; i++)
   {
      A (i) = A (i) + res.flow (i) * deltat;
     if (A(i)<0) A(i)=0;
    res.flow(i)=0;
    }
  return(res);
}
/*
function [res,Anew] = boxcartrain(boxnr,A, devrate, cv_devrate)
global g_grind;
if nargin < 4
   %if no sd at default a escalator boxtrain;
   cv_devrate = 0;
end;
%initiations;
%
deltat = solver('step');
N = size(A,1);
N2= size(A,2);
%gamma = 1 / (N-0.5); %%%% should last a bit longer as the last class transfers gradually into the next stage
 gamma = 1 / N; %%%% origial code of Goudriaan!!
if ~any(strcmpi(g_grind.solver.name, {'euler','c.euler'}))
   warning('GRIND:boxcartrain:euler','Boxcartrain needs Euler integration');
end;
%flow calculation
%
if isempty(cv_devrate) || (cv_devrate < 0)
   %
   %fixed boxcar train
   %
   %flows from each boxcar
   res.flow = A ./ gamma .* devrate;
   %shift the flows to get the inflow to the next boxcar and keep outflow
   res.outflow = res.flow(N,:);
   res.flow = [0; res.flow(1:N - 1,:)]-res.flow; %nett flow: what goes in what out?
elseif cv_devrate == 0
   %
   %escalator boxcar train
   %
   res.flow = zeros(N, N2);
   res.outflow = zeros(1,N2);
   if isnan(g_grind.boxcar.gcycl(boxnr))
      g_grind.boxcar.gcycl(boxnr) = 0.5 * gamma; % !! boxcar.gcycl should be reset to []
   else
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) + devrate * deltat;
   end;
%the last class empties gradually instead of all at once
% Goudriaan code
% if (gamma-boxcar.gcycl > 1E-30)
 %    res.outflow =  min(A(N)/deltat, A(N) / (gamma - boxcar.gcycl) * devrate);
%      res.flow(N) = -res.outflow;
% end;

   if g_grind.boxcar.gcycl(boxnr) >= gamma
      %flow is the difference between next and current
      res.outflow = res.outflow+A(N,:)/deltat; % all biomass is moved in one moment to the next stage
      res.flow = ([zeros(1,N2); A(1:N - 1,:)]-A) / deltat; %divide by deltat such that the whole A will be shifted
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) - gamma;
   end;
else
   %
   %fractional boxcar train
   %
   res.flow = zeros(N, N2);
   if isnan(g_grind.boxcar.gcycl(boxnr))
      g_grind.boxcar.gcycl(boxnr) = 0.5 * gamma; % !! g_grind.boxcar.gcycl should be reset to NaN
   else
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) + devrate * deltat;
   end;
   %
   %Fraction F to be moved (Goudriaan)
   F = 1 - (N * cv_devrate^2);
   if F<0
      error('GRIND:boxcartrain:CVtoobig','Boxcartrain: coefficent of variation (%g) too big or number of stages (%d) too big\nboxcar %s CV should at least be < %g (<sqrt(1/N))',cv_devrate,N,g_grind.boxcar.names{boxnr},sqrt(1/N));
   end;
   if deltat > F * gamma / (devrate + 1e-30)
      error('GRIND:boxcartrain:IntSteptoobig','Boxcartrain: Integration step too large: deltat should be < %0.5g\nAlternatively decrease the number of classes of %s to %d\ncv = %g devr = %g',...
         F * gamma / (devrate + 1e-30),g_grind.boxcar.names{boxnr},floor(1/(deltat*devrate+cv_devrate^2)),cv_devrate,devrate);
   end;
   if (gamma -g_grind.boxcar.gcycl(boxnr) > 1E-30)
      res.outflow = min(A(N,:)/deltat,A(N,:) / (gamma - g_grind.boxcar.gcycl(boxnr)) * devrate);
      res.flow(N,:) = -res.outflow;
   else
      res.outflow = zeros(1,N2);
   end;
   if g_grind.boxcar.gcycl(boxnr) >= F * gamma
      %flow is the difference between next and current
      res.flow = ([zeros(1,N2); A(1:N - 1,:)]-A) .* F ./ deltat; %divide by deltat such that the whole A will be shifted
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) - F * gamma;
   end;
end;
%update boxcar train
Anew = A + res.flow * deltat;
%g_grind.boxcar.trains{boxnr} = boxcar;

*/
