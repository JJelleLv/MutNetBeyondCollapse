%DEFINEPARS   Define parameters. 
%   Function used to define parameters in a model. Only necessary if you 
%   wish to add extra parameters to the model that are not used elsewhere.
%
%   Usage:
%   DEFINEPARS PAR1,PAR2,PAR3 etc. - This line is used only to detect the parameter names, 
%   it is not a working function.
%   DEFINEPARS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'par' [identifier] - define parameter par.
%
%
%   See also model
%
%   Reference page in Help browser:
%      <a href="matlab:commands('definepars')">commands definepars</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function definepars(varargin)
fieldnams={'par','U1','define parameter par.',''}';
i_parseargs(fieldnams,'par(+)','',varargin,false,{@i_isid});
