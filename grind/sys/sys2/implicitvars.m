%IMPLICITVARS - define a fully implicit differential equation DAE
%   Differential-algebraic equations (DAEs) can be defined using this command. GRIND then uses
%   the MATLAB solver <a href="matlab:help ode15i">ode15i</a> to solve a fully 
%   implicit equation of the following form:
%   0=F(y,y',t)
%  
%   Example:
%   Weissing implicit differential equation:
%   First define the state variables:
%   implicitvars y
%   Then the equation should be given starting with 0=..
%   0 = t*y^2 * y'^3 - y^3 * y'^2 + t*(t^2 + 1)*y' - t^2 * y;
%  
%   in case you have more than one components the equation should have the 
%   same number of elements as the number of statet variables:
%   implicitvars x y
%   0=[f(x,x',y,y');...
%   g(y,y',x,x')];
%   If you would need to use strings like <a href="matlab:help val">val('A')</a> in your definition, use 
%   double quotes val("A") to avoid confusion with derivative.
%  
%  See also model, modelpanel
%
%   Reference page in Help browser:
%      <a href="matlab:commands('implicitvars')">commands implicitvars</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function implicitvars(varargin)
disp('Use this function only in a definition of a model');