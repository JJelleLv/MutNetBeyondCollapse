 +++++++++++++++++++++ GRIND for MATLAB ++++++++++++++++++++++++
 Version 2.00
 November 2018
 MATLAB versions: R2008a to R2018b

 grind for MATLAB
    grind for MATLAB is a toolbox for analysing sets of differential equation, 
    and it also supports difference equations, delay differential equations (DDE), vectors, 
    matrices (cellular automata). It was originally based on grind, which was a convenient DOS program for analysing sets of 
    differential equations made by Rob de Boer.
        
    grind for MATLAB is a command based system, i.e. the user types commands
    in the MATLAB command Window to do most analyses and make figures. 
    These figures can be edited using standard MATLAB commands and menus. 
    Additionally, a user can click in some figures for instance to show 
    trajectories in a phase plane. 
    The system includes commands for:
    - simulating sets of n ordinary differential equations, including time 
      delays (see lag) and vector- and matrix notations. 
    - simulating sets of difference equations.
    - creating null-isoclines (nullclines) in phase spaces of 2 or 3 dimensions. (see null, null3)
    - creating nullclines in spaces spanned by state variables and parameters. 
    - simple one-dimensional bifurcation analysis by simulation. (see paranal)
    - continuation of equilibria using the engine of MATCONT or COCO. (see conteq) 
    - stability analysis of equilibria, using eigenvalues of the Jacobian matrix. (see eigen)
    - creating plots of user-defined functions (including parameters from 
      the current model)(see funplot)
    - using (red) noise and interpolation of external variables (for example 
      real temperature data) (see rednoise)
    - automatic calibration of parameters by optimizing the sum of squares 
      between observed and predicted values(see optimpars).
    - various special analyses such as determining the Lyapunov coefficient
      to detect chaos, Poincaré sections etc. (see lyapunov)
    - discrete events within a continuous model (see setevent)
    - stochastic differential equations (see dwiener)
    - three kinds of "boxcartrains" (Goudriaan, 1989) to model stage structured populations (see boxcartrain)
 
 
    Requirements: MATLAB R2008a or later, no MATLAB toolboxes are required, but for some commands/systems
    the Statistics and the Symbolic toolboxes are recommended, and optimization tools from the Optimization toolbox can optionally be used
    (see toolboxes)    
    
    
    Use the command model to create a model (which is saved to an ini
    file). Use the upper panel in the window to enter differential 
    equation or difference equations. You can also define the model as a Forrester diagram 
    (see vismod).
 
    Example:
    Logistic differential equation:
      N'=N*r*(1-N/K)
    Logistic difference equation (or recurrence relation):
       N(t+1)=N(t)*r*(1-N(t)/K)
 
    Note that parameters and state variables are case-sensitive.
    Parameter names can be all alphanumeric names except the
    following reserved words:
       t   time
      pi   pi=3.1416
      inf  Infinity
      Inf  Infinity
      nan  NaN (not-a-number)
      NaN  NaN (not-a-number)
      eps  Floating point relative accuracy eps=2.2204e-016 
 
    In the lower panel default values of the parameters and initial
    values of the state variables are entered.  
    Furthermore default commands can be entered here.
 
    Example:
      N=0.01;
      r=0.5;
      ax x N [0 10];
    (the semicolon is not required, but suppresses unnecessary
    output) 
 
    Author: 
 
    Egbert van Nes (egbert.vanNes@wur.nl)
    Wageningen University
    Aquatic Ecology and Water Quality Management group
    PO Box 47
    6700 AA Wageningen 
    The Netherlands
       
 
    grind for Matlab website:
    http://www.sparcs-center.org/grind/ 
    
+++++++++++++++++++++ INSTALLING GRIND ++++++++++++++++++++++++
 - Extract the files of GRIND.zip to any directory

 - to start GRIND:
 (1) start MATLAB
 (2) go to the directory where GRIND was extracted.
 (3) type CD GRIND
 (4) type SETUPGRIND (only the first time)
 (5) type MODEL and enter the model (an existing model can be opened
     with USE INIFILE)
 (6) type commands

     type HELP COMMANDS to get an overview of the GRIND commands.


 ++++++++++++++++++++++++ COMMANDS +++++++++++++++++++++++++++
 
 General information
   grind        - analysing difference and differential equations
      
 Define model:
   assign(=)    - Assign parameters/state variables 
   enterjac     - Enter the derivatives used to calculate the Jacobian matrix 
   funcs        - Add/edit functions that can be used for axes and funplot 
   initgrind    - Initiate global variables (required before running any model). 
   loaddata     - Load data file for parameter optimizing or external variables 
   model        - Create or open a model 
   savemodel    - Save the model to an inifile 
   savepar      - Save model and current parameter settings to a file 
   setdata      - Enter data for parameter optimizing or external variables 
   setodefile   - Use a specific ODEFILE (g_grind.odefile) 
   setupgrind   - Add the grind\sys directory to the default MATLAB path 
   use          - Load inifile and initialize the model 
   vismod       - Create a model as a Forrester diagram
   
 Simulating:
   backw        - Run the model with negated right hand sides 
   dirfield     - Create a direction field of t vs. state variable 
   e2n          - Combination of the commands ERA and NULL 
   e2p          - Combination of the commands ERA, PHAS 2, RU and NULL 
   e2r          - Combination of the commands ERA, PHAS 2 and RU 
   e3r          - Combination of the commands ERA, PHAS 3 and RU 
   era          - Erase and close all figures 
   gmex         - Compile/create mex ODE file
   itermap      - Plot x(t) versus x(t+n) for a 1D difference equation 
   ke           - Fix the final state of the last simulation as new starting point 
   null         - Create a phase plane with null isoclines 
   null3        - Create a 3D phase space with null isoclines 
   phas         - Create a 2d or 3d phase space 
   plotdiff     - Plot a one dimensional differential equation 
   plotreldiff  - Plot per capita growth of a one dimensional differential equation 
   ru           - Run the model and show the results in the current plot (or the phase plane) 
   rungrid      - Create grid of trajectories and show the results in the current plot (or the phase plane) 
   randinit     - Set random initial values. 
   stabil       - Run the current model without showing results and keep the final state as initial value 
   time         - Run the model and show the results in a time plot 
   timesens     - Sensitivity analysis during a time run 
   val          - Show initial and final values of state variables 
   vector       - Show vectors of change in the phase plane 
   vectplot     - Create special plots of vector state variables 
   viewcells    - Make a movie of matrix state variable (cellular automata) 
   where        - Show the current initial conditions in the phase plane 
   
 Analysing:
   attrbasin    - Find basin of attraction by simulation. 
   autocorrs    - Autocorrelation function of state variables 
   conteq       - Continue an equilibrium and locate some bifurcations 
   contbif      - Continue a bifurcation in two dimensions
   eigen        - Calculate the eigenvalues of the Jacobian in the current initial point 
   findeq       - Find closest equilibrium and set initial values to the equilibrium 
   funcs        - Add/edit functions that can be used for axes and funplot 
   funplot      - Create a plot with a user defined function (may include parameters and state vars) 
   growths      - Make filled contour plots with the growth of each state variable 
   implicitplot - Create a 2D implicit plot of a unsolved function written as f(x,y)=0 (may include pars and state vars) 
   lyapspect    - Determine the whole Lyapunov exponent spectrum  
   lyapunov     - Determine the Lyapunov exponent (?), expressing the sensitivity to initial conditions  
   lyapsimple   - Simple Ellner algorithm for Lyapunov exponent
   locallyap    - Ellner algorithm for local Lyapunov exponent
   lorenzmap    - Create a Lorenz map (the maxima of a state variable versus the next maxima) 
   marbleplot   - Calculate marble-in-a-cup plot (1D model only) 
   manifolds     - Trajectories in the direction of the eigenvectors
   optimpars    - Optimize parameters (use setdata to load data first) 
   paranal      - Change one parameter step-by step and show the attractor by simulation 
   paranal2d    - Change two parameters step-by step and show the mean attractor by simulation 
   perturb      - Perturb a saddle in the direction of the attracting eigenvector (combine with backw to plot separatrix) 
   potential    - Create potential (marble-in-a-cup) plot (works only in 1D models correctly) 
   poincaremap  - Make Poincaré map of previous simulation run 
   poincaresect - Make Poincaré section of previous simulation run. 
   replayall    - Replay the model updating one or more graphs
   returntime   - Estimates the time necessary to reach a stable node 
   returntime2d - Plot of the time necessary to reach a stable node. 
   takens       - Generates a Takens reconstruction of the variable. 
   torus        - Analyse cycles by plotting in polar coordinates. 
   trdet        - Determine trace and determinant of equilibrium    
      
 Settings:
   addmode      - Sets the simulation modus to adding or not 
   ax           - Define an axis for the phase plane 
   enterJac     - Enter the derivatives used to calculate the Jacobian matrix 
   funcs        - Add/edit functions that can be used for axes and funplot 
   nextpen      - Assign the next colour for g_grind.pen (the colour of the next trajectory in the phase plane). 
   out          - Sets variables/functions for time plots (at default all state variables) 
   par          - Show values of parameters 
   savemodel    - Save the model to an inifile 
   savepar      - Save model and current parameter settings to a file 
   setdata      - Enter data for parameter optimizing 
   setdefaults  - Save or get default settings 
   setpen       - Set style and colour of runs in the phase plane (and dirfield) 
   setcolormap  - Set the colour gradient in 3D plots 
   setupgrind   - Add the grind\sys directory to the default MATLAB path (see pathdef) 
   simtime      - Set the simulation parameters (t g_grind.ndays g_grind.tstep) 
   solver       - Select a solver and solver settings (only for ODE's) 
   stat         - Display the complete status of GRIND (global variables, parameters and state variables) 
      
 Special functions:
   arrows       - Add arrows to a figure. 
   boxcartrain  - Update a boxcartrain variable. 
   boxcarinflow  - Inflow into a boxcartrain
   combfig      - Combine two or more figures in one figure. 
   copyfig      - Copy a figure to a new figure on screen. 
   defextern    - Use an external variable as parameter (interpolate values in matrix) 
   figdata      - Extract data from the currently selected figure 
   fftconv       - Convolution using fft (for Integro-Difference Equations)
   grindpath    - Function that returns a full string to the grind system directory 
   iif          - Immediate if, handy way of using if conditions 
   insim        - Phenology and populatIoN SIMulator (Mols)
   lag          - Function for using a time lag on one of the state variables in a differential equation 
   leftcells    - Shift a matrix one position to the left, see also rightcells, upcells, downcells. 
   outf         - various function (e.g. sum, mean, min, max) of (vector) state variable. 
   rednoise     - Function that generates red or white noise (Steele)  
   setevent     - Create a discrete event queue. 
   setmat       - Set the values of a matrix (disturbance/gradient)
   setdimension - Change the dimension of a (state) variable
   sp_corr      - Function to calculate the spatial correlation of a matrix for
   varcopy      - Copy a matrix to the clipboard. 
   varpaste     - Paste a matrix from the clipboard. 
