#include "csolver1.h"
#include <string>
#include <string.h>
#include <iostream>
#include <ctime>

using namespace std;
string helpmodel();
int main (int argc, char* argv[])
//Main routine, just reads the arguments and opens/runs the model

{
// time_t tstart = clock();
   try
   {
      string filename, psthname;
      string pathname = argv[0];
      int p = pathname.rfind ("\\");
      pathname = pathname.substr (0, p + 1);
      if (argc == 2)
      {
         if (strcmp (argv[1], "?") == 0)
         {
            cout << helpmodel();
            return 0;
         }
         filename = argv[1];
      }
      else
      {
         filename = pathname + "param1.tmp";
      }
      TModel TheModel (filename);

//cout << "Reading pars took "<< difftime(clock(), tstart)/1000 <<" second(s)."<< endl;
      TheModel.runsolver();

//cout << "Total "<< difftime(clock(), tstart)/1000 <<" second(s)."<< endl;
      return 0;
   }
   catch ( string err )
   {
      cout << "Exception raised: " << err << '\n';
      return 1;
   }
}
