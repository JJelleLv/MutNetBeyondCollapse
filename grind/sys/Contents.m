%GRIND for MATLAB
%
%General information
%  grind        - analysing difference and differential equations
%     
%Define model:
%  assign(=)    - Assign parameters/state variables
%  defextern    - Use an external variable as parameter (interpolate values in matrix)
%  definepars   - Define extra parameters that are not in any of the equations
%  definespace  - Define lattice to add space
%  defpermanent - Define permanent variables (variables that can be changed during a run)
%  dwiener      - Define a Wiener process (= stochastic differential equation)
%  djump        - Define a jump process in a stochastic differential equation
%  enterjac     - Enter the derivatives used to calculate the Jacobian matrix
%  funcs        - Add/edit functions that can be used for axes and funplot
%  initgrind    - Initiate global variables (required before running any model). 
%  loaddata     - Load data file for parameter optimizing or external variables
%  makemap      - Make a discrete model of a cyclic differential equation
%  model        - Create or open a model
%  onstart      - Define a callback function that is called before each run
%  savemodel    - Save the model to an inifile
%  savepar      - Save model and current parameter settings to a file
%  setdata      - Enter data for parameter optimizing or external variables
%  setodefile   - Use a specific ODEFILE (g_grind.odefile)
%  setupgrind   - Add the grind\sys directory to the default MATLAB path
%  use          - Load inifile and initialize the model
%  vismod       - Create a model as a Forrester diagram
% 
%Simulating:
%  backw        - Run the model with negated right hand sides 
%  dirfield     - Create a direction field of t vs. state variable
%  e2n          - Combination of the commands ERA and NULL
%  e2p          - Combination of the commands ERA, PHAS 2, RU and NULL
%  e2r          - Combination of the commands ERA, PHAS 2 and RU
%  e3r          - Combination of the commands ERA, PHAS 3 and RU
%  era          - Erase and close all figures
%  gmex         - Compile/create mex ODE file
%  itermap      - Plot x(t) versus x(t+n) for a 1D difference equation
%  ke           - Fix the final state of the last simulation as new starting point
%  null         - Create a phase plane with null isoclines
%  null3        - Create a 3D phase space with null isoclines
%  phas         - Create a 2d or 3d phase space
%  plotdiff     - Plot a one dimensional differential equation
%  plotreldiff  - Plot per capita growth of a one dimensional differential equation
%  ru           - Run the model and show the results in the current plot (or the phase plane)
%  rungrid      - Create grid of trajectories and show the results in the current plot (or the phase plane)
%  randinit     - Set random initial values.
%  stabil       - Run the current model without showing results and keep the final state as initial value
%  time         - Run the model and show the results in a time plot
%  timesens     - Sensitivity analysis during a time run
%  val          - Show initial and final values of state variables
%  vector       - Show vectors of change in the phase plane
%  vectplot     - Create special plots of vector state variables
%  viewcells    - Make a movie of matrix state variable (cellular automata)
%  where        - Show the current initial conditions in the phase plane
% 
%Analysing:
%  analunits    - Analyse and add units to the parameters and state variables.
%  attrbasin    - Find basin of attraction by simulation.
%  autocorrs    - Autocorrelation function of state variables
%  conteq       - Continue an equilibrium or bifurcation using MATCONT or COCO
%  eigen        - Calculate the eigenvalues of the Jacobian in the current initial point
%  findeq       - Find closest equilibrium and set initial values to the equilibrium
%  fokkerplanck - Analyse a 1D model as a Fokker-Planck equation
%  funcs        - Add/edit functions that can be used for axes and funplot
%  funplot      - Create a plot with a user defined function (may include parameters and state vars)
%  growths      - Make filled contour plots with the growth of each state variable
%  implicitplot - Create a 2D implicit plot of a unsolved function written as f(x,y)=0 (may include pars and state vars)
%  lyapspect    - Determine the whole Lyapunov exponent spectrum 
%  lyapunov     - Determine the Lyapunov exponent (?), expressing the sensitivity to initial conditions 
%  lorenzmap    - Create a Lorenz map (the maxima of a state variable versus the next maxima)
%  mcarlo       - Draw parameters from distributions for sensitivity or uncertainty analysis
%  marbleplot   - Calculate marble-in-a-cup plot (1D model only)
%  manifolds    - Trajectories in the direction of the eigenvectors
%  optimpars    - Optimize parameters (use setdata to load data first)
%  outfun       - Extract data from last run (any equation)
%  paranal      - Change one parameter step-by step and show the attractor by simulation
%  paranal2d    - Change two parameters step-by step and show the mean attractor by simulation
%  perturb      - Perturb a saddle in the direction of the attracting eigenvector (combine with backw to plot separatrix) 
%  potential    - Create potential (marble-in-a-cup) plot (works only in 1D models correctly)
%  poincaremap  - Make Poincaré map of previous simulation run
%  poincaresect - Make Poincaré section of previous simulation run.
%  replayall    - Replay the model updating one or more graphs
%  returntime   - Estimates the time necessary to reach a stable node
%  returntime2d - Plot of the time necessary to reach a stable node.
%  takens       - Generates a Takens reconstruction of the variable.
%  torus        - Analyse cycles by plotting in polar coordinates.
%  trdet        - Determine trace and determinant of equilibrium   
%    
%Settings:
%  addmode      - Sets the simulation modus to adding or not 
%  ax           - Define an axis for the phase plane
%  enterJac     - Enter the derivatives used to calculate the Jacobian matrix
%  funcs        - Add/edit functions that can be used for axes and funplot
%  gstat        - Display the complete status of GRIND (global variables, parameters and state variables)
%  nextpen      - Assign the next colour for g_grind.pen (the colour of the next trajectory in the phase plane).
%  out          - Sets variables/functions for time plots (at default all state variables)
%  par          - Show values of parameters
%  savemodel    - Save the model to an inifile
%  savepar      - Save model and current parameter settings to a file
%  setdata      - Enter data for parameter optimizing
%  setdefaults  - Save or get default settings
%  setpen       - Set style and colour of runs in the phase plane (and dirfield)
%  setcolormap  - Set the colour gradient in 3D plots
%  setupgrind   - Add the grind\sys directory to the default MATLAB path (see pathdef)
%  simtime      - Set the simulation parameters (t g_grind.ndays g_grind.tstep)
%  solver       - Select a solver and solver settings (only for ODE's)
%    
%Special functions:
%  arrows       - Add arrows to a figure.
%  boxcartrain  - Update a boxcartrain variable.
%  boxcarinflow - Inflow into a boxcartrain
%  connectmat   - Make a matrix for connections between cells of a lattice.
%  combfig      - Combine two or more figures in one figure.
%  copyfig      - Copy a figure to a new figure on screen.
%  figdata      - Extract data from the currently selected figure
%  fftconv      - Convolution using fft (for Integro-Difference Equations)
%  externlag    - Function for using a time lag on one of external variables
%  grindpath    - Function that returns a full string to the grind system directory
%  iif          - Immediate if, handy way of using if conditions
%  insim        - Phenology and populatIoN SIMulator (Mols)
%  lag          - Function for using a time lag on one of the state variables in a differential equation
%  leftcells    - Shift a matrix one position to the left, see also rightcells, upcells, downcells.
%  outf         - various function (e.g. sum, mean, min, max) of (vector) state variable.
%  rednoise     - Function that generates red or white noise (Steele) 
%  setevent     - Create a discrete event queue.
%  setmat       - Set the values of a matrix (disturbance/gradient)
%  setdimension - Change the dimension of a (state) variable
%  sp_corr      - Function to calculate the spatial correlation of a matrix for
%  varcopy      - Copy a matrix to the clipboard.
%  varpaste     - Paste a matrix from the clipboard.
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('Contents')">commands Contents</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
%start


