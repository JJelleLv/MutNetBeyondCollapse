%GRIND for MATLAB
%   GRIND for MATLAB is a toolbox for analysing sets of differential equation, 
%   and it also supports difference equations, delay differential equations (DDE), vectors, 
%   matrices (cellular automata). It was originally based on GRIND, which was a convenient DOS program for analysing sets of 
%   differential equations made by Rob de Boer.
%      
%   GRIND for MATLAB is a command based system, i.e. the user types commands
%   in the MATLAB command Window to do most analyses and make figures. 
%   These figures can be edited using standard MATLAB commands and menus. 
%   Additionally, a user can click in some figures for instance to show 
%   trajectories in a phase plane.
%   The system includes commands for:
%   - simulating sets of n ordinary differential equations, including time 
%     delays (see <a href="matlab:help lag">lag</a>) and vector- and matrix notations. 
%   - simulating sets of difference equations.
%   - creating null-isoclines (nullclines) in phase spaces of 2 or 3 dimensions. (see <a href="matlab:help null">null</a>, <a href="matlab:help null3">null3)</a>
%   - creating nullclines in spaces spanned by state variables and parameters.
%   - simple one-dimensional bifurcation analysis by simulation. (see <a href="matlab:help paranal">paranal</a>)
%   - continuation of equilibria using the engine of MATCONT or COCO. (see <a href="matlab:help conteq">conteq</a>) 
%   - stability analysis of equilibria, using eigenvalues of the Jacobian matrix. (see <a href="matlab:help eigen">eigen</a>)
%   - creating plots of user-defined functions (including parameters from 
%     the current model)(see <a href="matlab:help funplot">funplot</a>)
%   - using (red) noise and interpolation of external variables (for example 
%     real temperature data) (see <a href="matlab:help rednoise">rednoise</a>)
%   - automatic calibration of parameters by optimizing the sum of squares 
%     between observed and predicted values(see <a href="matlab:help optimpars">optimpars</a>).
%   - various special analyses such as determining the Lyapunov coefficient
%     to detect chaos, Poincar� sections etc. (see <a href="matlab:help lyapunov">lyapunov</a>)
%   - discrete events within a continuous model (see <a href="matlab:help setevent">setevent</a>)
%   - stochastic differential equations (see <a href="matlab:help dwiener">dwiener</a>)
%   - three kinds of "boxcartrains" (Goudriaan, 1989) to model stage structured populations (see <a href="matlab:help boxcartrain">boxcartrain</a>)
%
%
%   Requirements: MATLAB R2008a or later, no MATLAB toolboxes are required, but for some commands/systems
%   the Statistics and the Symbolic toolboxes are recommended, and optimization tools from the Optimization toolbox can optionally be used
%   (see <a href="matlab:commands toolboxes">toolboxes</a>)    
% 
% 
%   Use the command <a href="matlab:help model">model</a> to create a model (which is saved to an ini
%   file). Use the upper panel in the window to enter differential 
%   equation or difference equations. You can also define the model as a Forrester diagram 
%   (see <a href="matlab:help vismod">vismod</a>).
%
%   Example:
%   Logistic differential equation:
%     N'=N*r*(1-N/K)
%   Logistic difference equation (or recurrence relation):
%      N(t+1)=N(t)*r*(1-N(t)/K)
%
%   Note that parameters and state variables are case-sensitive.
%   Parameter names can be all alphanumeric names except the
%   following reserved words:
%     <a href="matlab:commands t">t</a>   time
%     pi   pi=3.1416
%     inf  Infinity
%     Inf  Infinity
%     nan  NaN (not-a-number)
%     NaN  NaN (not-a-number)
%     eps  Floating point relative accuracy eps=2.2204e-016 
%
%   In the lower panel default values of the parameters and initial
%   values of the state variables are entered. 
%   Furthermore default commands can be entered here.
%
%   Example:
%     N=0.01;
%     r=0.5;
%     ax x N [0 10];
%   (the semicolon is not required, but suppresses unnecessary
%   output) 
%
%   Author: 
%
%   Egbert van Nes (<a href="mailto:egbert.vanNes@wur.nl">egbert.vanNes@wur.nl</a>)
%   Wageningen University
%   <a href="http://www.aew.wur.nl/">Aquatic Ecology and Water Quality Management group</a>
%   PO Box 47
%   6700 AA Wageningen 
%   The Netherlands
%    
%
%   GRIND for Matlab website:
%   <a href="http://www.sparcs-center.org/grind/">http://www.sparcs-center.org/grind/</a>
% 
%   See also installing, model
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('grind')">commands grind</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
%start
 
