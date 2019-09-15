%FUNCS   Edit functions
%   Edit/add "functions" (or auxiliary variables) for usage on the axes or for function plots
%   (<a href="matlab:help funplot">funplot</a>). The model is not changed. The functions may not have arguments.
%
%   Examples:
%   cons=g*Z*A/(A+H);
%   You may use IF FOR WHILE and SWITCH statements:
%   'Tent function:'
%   if (x(t)<=0.5)
%       fx=r*x(t)
%   else
%       fx=r-r*x(t)
%   end
%   x(t+1)=fx
%   Note that for modelling equations like this, the usage of intermediary variables
%   is required.  
%
%   Usage:
%   FUNCS - opens a screen to edit functions
%   FUNCS FUN1 FUN2 - adds FUN1 and FUN2 to the functions without opening the edit screen.
%   FUNCS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - the equation to add
%     
%
%   See also funplot, ax, out
%
%   Reference page in Help browser:
%      <a href="matlab:commands('funcs')">commands funcs</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function funcs(varargin)
global g_grind;
eoln = sprintf('\n');  %#ok<SPRINTFN>
if nargin>0
   fieldnams={'fun', 's', 'the equation to add',''}';
   args=i_parseargs(fieldnams,'fun(+)','',varargin);
   i_parcheck;
   if ~iscell(args.fun)
      args.fun={args.fun};
   end
   for i=1:length(args.fun)
      if isempty(strfind(args.fun{i},';'))
         g_grind.funcs=[g_grind.funcs eoln args.fun{i} ';'];
      else
         g_grind.funcs=[g_grind.funcs eoln args.fun{i}];
      end
   end
   updatefuncs;
   return;
else
    i_parcheck;
end
if isempty(g_grind.funcs)
   g_grind.funcs=eoln;
end
answer=inputdlg({'Edit functions'},'Functions for output (doesn''t change model)',20,{g_grind.funcs});
if ~isempty(answer)
   a = answer{1};
   g_grind.funcs = [];
   for i = 1:size(a, 1)
      a1 = strtrim(a(i, :));
      if ~isempty(a1)
         if ~strcontains(a1,';')
            g_grind.funcs = [g_grind.funcs a1 ';' eoln];
         else
            g_grind.funcs = [g_grind.funcs a1 eoln];
         end
      end
   end
   g_grind.funcs = g_grind.funcs(1:length(g_grind.funcs) - length(eoln));
   updatefuncs;
end
disp(g_grind.funcs);
function updatefuncs
global g_grind;
funcs=str2cell(g_grind.funcs);
if isfield(g_grind.funcnames,'dims')
   g_grind.funcnames=rmfield(g_grind.funcnames,'dims');
end
for i=1:length(funcs)
    f=strfind( funcs{i},'=');
   if ~isempty(f)
      funname = strtrim(funcs{i}(1:f(1) - 1));
      j=1;
      while (j<=length(g_grind.funcnames.names)) &&~strcmp(g_grind.funcnames.names{j},funname)
         j=j+1;
      end
      if j>length(g_grind.funcnames.names)
         g_grind.funcnames.names=[g_grind.funcnames.names {funname}];
      end
   end
end
clear g_funcs;
