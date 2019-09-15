%RANDINIT   Random initial values
%   Set random (uniformly distributed) initial values.
%  
%   Usage:
%   RANDINIT - chooses random initial points between 0 and 100.
%   RANDINIT MAX - chooses random initial points between 0 and MAX.
%   RANDINIT [MIN MAX] - chooses random initial points between MIN and MAX.
%   RANDINIT VAR - chooses random initial points of state variable VAR between 0 and 100.
%   RANDINIT VAR [MIN MAX] - chooses random initial points of VAR between MIN and MAX.
%   RANDINIT VAR MAX - chooses random initial points of VAR between 0 and MAX.
%   RANDINIT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'range' [number and length(number)<=2] - range for the initial points
%     'var' [state variable or empty] - state variable (empty=both)
%
%
%   See also ke, val
%
%   Reference page in Help browser:
%      <a href="matlab:commands('randinit')">commands randinit</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function randinit(varargin)
%(VarName,range)
global g_grind g_Y;
fieldnams={'var', 'v#E', 'state variable (empty=both)',[];...
   'range', 'n&length(n)<=2', 'range for the initial points',[0 100']}';
args=i_parseargs(fieldnams,'if(argtype(1,''v'')),deffields=''var,range'';else,deffields=''range'';end;','',varargin);
i_parcheck;
if ~isfield(args,'var')
    args.var=[];
end
if ~isfield(args,'range')||isempty(args.range)
    args.range=[0 100];
end
if numel(args.range)==1
    args.range=[0 args.range];
end

if isempty(args.var)
   N0=drawuniform(g_grind.statevars.dim,1,args.range(1),args.range(2));
else
   N0=i_initvar;
   [sfrom, sto] = i_statevarnos(args.var);
   if ~isempty(sfrom)
       N0(sfrom:sto)=drawuniform(sto-sfrom+1,1,args.range(1),args.range(2));
   end
end
i_keep(N0);
g_Y=[];
