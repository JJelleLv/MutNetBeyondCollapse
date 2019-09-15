%LAG   State variable with time lag in a delay differential equation.
%   Use a time lag for some state variable (e.g. May 1976). Use
%   this function in differential equations. The delay differential
%   equation (DDE) is solved with the dde23 or ddesd solver.
%   You can use the same function among others in time plots (see <a href="matlab:help out">out</a>) to plot 
%   lagged state variables, auxiliary variables or functions.
%
%   Usage:
%   LAG(statevar,TIMEDELAY)
%   statevar = name of state variable.
%   TIMEDELAY   = time delay (can be a equation including parameters, time (t) and state variables). If the time 
%   delay is variable <a href="matlab:help ddesd">ddesd</a> is used. Note that this equation should return the delays: GRIND adds
%   "t-" to these delays as needed in ddesd.
%
%See also model, <a href="matlab:help dde23">dde23</a>,<a href="matlab:help ddesd">ddesd</a> out
%
%   Reference page in Help browser:
%      <a href="matlab:commands('lag')">commands lag</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res=lag(var,timelag)
global g_t;
if ischar(var)
    if ~strncmp(var,'??',2)
       var=outfun(var);
    else
       help lag
       return;
    end
end
if size(var,1)~=size(g_t)
   error('GRIND:lag:ArgError','Error in lag function, var should be a string or a matrix of size g_t');
end
res=interp1(g_t,var,g_t-timelag);


