%EXTERNLAG   Get the value of an external variable with a time lag
%   This function can evaluate an external variable but with a time lag back in time. You can use this 
%   function in your model equation. First define external variables using <a href="matlab:help defextern">defextern</a> or in <a href="matlab:help vismod">vismod</a>.
%   The default value of the external variable is used when there is no value defined.
%   Note, if the first argument is not an external variable, the function will be ignored. If played back
%   in equations (for instance in <a href="matlab:help outfun">outfun</a>), externlag returns some NaN values at the start of the run
%   and will ignore variable time lags. Use <a href="matlab:help defpermanent">permanent variables</a> if you need the exact values
%   of the lagged variables. 
%
%   Usage:
%   EXTERNLAG(var,lag) - the function has two parameters: (1) the name of the external variable. 
%   2) the time lag in time steps (can be a parameter or an equation). A positive value of the time 
%   lag looks back in time.
%
%
%   See also defextern, externvar, setdata, loaddata, model, lag
%
%   Reference page in Help browser:
%      <a href="matlab:commands('externlag')">commands externlag</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res=externlag(var,timelag)
global g_t;
if ischar(var)
    if ~strncmp(var,'??',2)
       var=outfun(var);
    else
       %help externlag
       return;
    end
end
if size(var,1)~=size(g_t)
   error('GRIND:extenlag:ArgError','Error in externlag function, var should be a string or a matrix of size g_t');
end
res=interp1(g_t,var,g_t-timelag);

