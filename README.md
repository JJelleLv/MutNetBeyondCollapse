MutNetBeyondCollapse - A MATLAB(R) Package
===================================

Copying and distribution of this file, with or without modification, are permitted in any medium without royalty provided the copyright notice and this notice are preserved. This file is offered as-is, without any warranty.

-------------------------------------------------------------------

When using this code for further research, please cite the original publication: Lever, J.J., Van de Leemput, I.A., Weinans, E., Quax, R., Van Nes, E. H., Bascompte, J., & Scheffer, M. (2019). Foreseeing the future of mutualistic communities beyond collapse. *Ecology Letters*, XX(X), XXX-XXX.

-------------------------------------------------------------------

#### PURPOSE AND SETUP

The main purpose of this code is to provide insight in the analysis done for the above mentioned publication. Adjustments to the code will likely be nescessary when applying the method proposed in this publication to your own data or when making additional simulations. Please contact the authors for more information.

The folder '+code' contains the Matlab functions supporting the main function (MAIN_MutNetBeyondCollapse).  
The folder 'data' contains data needed to make simulations, the output from simulations and time-series analysis.  
The folder 'grind' contains Grind for MATLAB, a package needed to generate time series (https://www.sparcs-center.org/resources/, downloaded d.d. 15/09/2019).  
The file 'MAIN_MutNetBeyondCollapse.m' is the main function executing six steps from making a simulation to plotting figures for initial network 'NETnr' and final network 'SETnr'.  

-------------------------------------------------------------------

#### NOTES ON MAIN_MutNetBeyondCollapse.m

To run the main function: open Matlab, go to the folder where the main function is stored, and type 'MAIN_MutNetBeyondCollapse(NETnr,SETnr)' in the command line, i.e. to make simulations, analysis, or plot figures for NETnr=1 and SETnr=1, type: 'MAIN_MutNetBeyondCollapse(1,1)'. 

In total examples of three initial networks (max NETnr=3) and five final networks per initial network (max SETnr=5) are provided.

Simulations or analysis are loaded from a file when replace\*File = false. 
New simulations or analysis are loaded from a file when replace\*File = true.

Step 1: Simulations are made with a length of 20.000 time steps. This step assumes that Grind for Matlab is installed (see above). When this is not the case, please go to the grind folder and run 'setupgrind' (see 'readme.txt' in this folder for additional information).  
Step 2: The mean abundance and a CI are determined in a rolling window with length 'meanRange_window' (in number of time steps). This window is moved forward with a number of steps equal to 'meanRange_stepsize'. Windows (partially) overlap when meanRange_stepsize &lt; meanRange_window.  
Step 3: The point at which a major shift in abundance occurs is determined (a Tipping Point), the observed change in abundance, and the direction and magnitude of the indicator prior to this tipping point. 'analyseTPP_window_dAb' specifies the size of the window within which the mean abundance is determined. 'analyseTPP_window_PCA' specifies the size of the window within which direction and magnitude of the indicator is determined prior to a tipping point. Both window are moved forward with a number of steps equal to 'analyseTPP_stepsize' (see Methods and Appendix S4).  
Step 4: The direction and magnitude of the indicator are determined for the full time series. 'makePCAskew_window' specifies the size of the window within which direction and magnitude of the indicator is determined. This window is moved forward with a number of steps equal to 'makePCAskew_stepsize'. Windows (partially) overlap when meanRange_stepsize &lt; meanRange_window (see Methods and Appendix S4).  
Step 5: Analysis of the trend in explained variance. A 'critical range' is found when there is a significant increase in explained variance. This trend is determined with a moving window with initial length 'ExplVartrend_window'. If an upward trend is detected within this window the length of the window is increased (see Methods and Appendix S4).  
Step 6: Data is plotted as in Fig. 2 of the above mentioned publication.