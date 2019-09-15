%ANALUNITS   Analysis of the units in the model
%   This function finds if there are inconsistencies in the units of the model and 
%   can append undefined units to the model. Note that you have to save the model with
%   <a href="matlab:help savemodel">savemodel</a> or <a href="matlab:help savepar">savepars</a> to keep the changes.
%  
%   Limitations: it cannot assume a unit for constants and cannot evaluate x^p 
%   (only x^2 or so).  In these cases these parts are neglected and the remaining 
%   part of the equation can still be evaluated.
%   SI units can be compared. Units can entered in various ways, for example:
%   g/m3 = g m-3
%   mug m-3 = µg m-3 (ug is not allowed) 
%   g/m3/d = g m-3 d-1 (multiple / are analysed from left to right)
%   g.m-3.d-1 = g m-3 d-1 (dots are allowed)
%   g/(m^3 d) = g m-3 d-1 (brackets are allowed, nesting is not recommended)
%   Non-standard IS notation is not supported, but you can add some comments between square brackets. 
%   All text between square brackets is ignored. For instance: g P /m3 is not allowed 
%   but g[P]/m3 is allowed. Also inconsequent use of units kg + 1000*g will give an error.
%
%
%   Usage:
%   
%   ANALUNITS - the units can be added/ removed in a dialog box. Press the analyse 
%   units button to run the analysis.
%   [allvars, allunits]=ANALUNITS(equations, vars, units) - Check and add the units of eqns (' is 
%   derivative to time).
%   ANALUNITS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'checkonly' [logical] - true is check only no change of units
%     'debug' [logical] - true is debug mode
%     'equations' [cell] - equations to analyse
%     'silent' [logical] - true is no output
%     'units' [cell] - units
%     'vars' [cell] - variables
%   ANALUNITS('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - check the current model (without adding)
%     '-d' - debug mode, gives information on the iterations used
%     '-l' - list    
%     '-s' - silent
% 
%   Example:
%   analunits({'cons = g*Z*(A/(A+h))+tt;','A''=r*A*(1-A/K)-cons','Z''=e*cons-m*Z;'},
%   {'A','Z','K','r'},{'g/m3','g/m3','g/m3','d'})
%   returns:
%   A [g m-3],  K [g m-3], Z [g m-3], cons [g m-3 d], e [-], g [d], h [g m-3], m [d]
%   r [d], t [d-1], tt [g m-3 d] 
%
%   See also par, savemodel, savepar      
%    
%
%   Reference page in Help browser:
%      <a href="matlab:commands('analunits')">commands analunits</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function [un, vars] = analunits(varargin)
%(mod, varlist, unitlist)
global g_grind;
fieldnams={'equations','c','equations to analyse','';...
   'vars','c','variables','';...
   'units','c','units','';...
   'debug','l','true is debug mode',false;...
   'silent','l','true is no output',false;...
   'checkonly','l','true is check only no change of units',false}';
args=i_parseargs(fieldnams,'equations,vars,units,debug','-d,-c,-s,-l',varargin);
   if any(strcmp(args.opts, '-d'))
      args.debug = true;
   elseif ~isfield(args,'debug')
       args.debug=false;
   end
   if any(strcmp(args.opts, '-c'))
      args.checkonly = true;
   elseif ~isfield(args,'checkonly')
       args.checkonly=false;
   end
   if any(strcmp(args.opts, '-s'))
      args.silent = true; 
   elseif ~isfield(args,'silent')
       args.silent=false;
   end
dolist= any(strcmp(args.opts, '-l'));
if ~args.checkonly&&~isfield(args,'equations')&&~dolist
   if ~isempty(g_grind)
      if args.debug
         par('-e','-d')
      else
         par('-e')
      end
      return;
   else
      prompt={'Enter an equation(s) ; delimited','List of (some) variables (space delimited)','List of units of these variables  (space delimited)'};
      name = 'Analyse units';
      defaultanswer={'','',''};
      answer = inputdlg(prompt, name, 1, defaultanswer);
      s=strrep(answer{1},';','\n');
      args.equations = parsed_equation(str2cell(sprintf(s)));
      disp(args.equations)
      s=strrep(answer{2},' ','\n');
      s=strrep(s,',','\n');
      s=strrep(s,'\n\n','\n');
      args.vars = str2cell(sprintf(s));
      s=strrep(answer{3},' ','\n');
      s=strrep(s,',','\n');
      s=strrep(s,'\n\n','\n');
      args.units = str2cell(sprintf(s));
   end
end
if ~isfield(args,'equations')||isempty(args.equations)&&~dolist
    if g_grind.solver.isimplicit
        error('grind:analunits','Full analysis of units of implicit models (DAE) is not implemented, test the equations seperately')
    end
    args.equations=parsed_equation(g_grind.model);
end
if ~isfield(args,'units')||isempty(args.units)
   args.units = par('-u');
   if ischar(args.units)
       args.units={args.units};
   end
   n = length(args.units);
   if g_grind.statevars.vector
      for i = 1:length(g_grind.statevars.vectnames)
         args.units{n + i} = par('-u', g_grind.statevars.vectnames{i});
      end
   else
      for i = 1:length(g_grind.statevars.names)
         args.units{n + i} = par('-u', g_grind.statevars.names{i});
      end
   end
   if all(strcmp(args.units, '')) %shortcut if none of the units are defined
      if ~args.silent
         fprintf('Cannot check units, as none of units are defined\n');
      end
      if nargout > 0
         un = [];
         vars = {};
      end
      return;
   end
end



if ~isfield(args,'vars')||isempty(args.vars)
   args.vars = g_grind.pars;
   n = length(args.vars);
   if g_grind.statevars.vector
      for i = 1:length(g_grind.statevars.vectnames)
         args.vars(n + i) = g_grind.statevars.vectnames(i);
      end
   else
      for i = 1:length(g_grind.statevars.names)
         args.vars(n + i) = g_grind.statevars.names(i);
      end
   end
end
if dolist
   [~, ndx] = sort(upper(args.vars));
   args.vars = args.vars(ndx);
   args.units = args.units(ndx);
   for i = 1:length(args.vars)
      v = varunit(args.units{i});
      if strcmp(char(v), args.units{i})||(isempty(args.units{i})&&strcmp(char(v), '[]'))
         fprintf('%s [%s]\n', args.vars{i}, args.units{i});
      else
         fprintf('%s [%s], interpreted as [%s]\n', args.vars{i}, args.units{i}, char(v));
      end
   end
   for i = 1:length(args.vars)
      v = varunit(args.units{i});
      [~, a] = v.basic;
      if ~isempty(a)
         warning('analunits:unknown','Variable %s [%s] has unknown unit(s) (%s)',args.vars{i},args.units{i},strtrim(sprintf('%s ',a{:})));
      end
   end
   return;
end
if ~isa(args.equations, 'parsed_equation')
   args.equations = parsed_equation(args.equations);
end

%' (before = sign) is assumed to be derivative to t
for i = length(args.equations):-1:1
   if ~any(strcmp(args.equations(i).fs, '='))
      args.equations(i) = [];
   else
      fs1 = args.equations(i).fs;
      diffeq = strcmp(fs1, 't');
      if any(diffeq)
         %remove (t+1) (t) or (t-1)
         diffeq=diffeq|strcmp(fs1,'(')|strcmp(fs1,')')|strcmp(fs1,'-')|strcmp(fs1,'+')|strcmp(fs1,'1')|strcmp(fs1,'1');
         types = zeros(size(diffeq));
         for k = 1:length(diffeq)
            if diffeq(k)
               if n > 0||strcmp(fs1{k}, '(')
                  types(k) = 1;
                  n = n + 1;
                  if strcmp(fs1{k}, ')')&&n > 2
                     types(types == 1)=2;
                     n = 0;
                  end
               end
            else
               types(types == 1)=0;
               n = 0;
            end
         end
         fs1=fs1(types ~= 2);
      end
      f1=find(strcmp(fs1,''''),1);
      f2=find(strcmp(fs1, '='), 1);
      f3 = find(strcmp(fs1,';'), 1);
      if isempty(f3)
         f3 = length(fs1) + 1;
      end
      f4 = find(strncmp(fs1, '%', 1), 1);
      if ~isempty(f4)
         f3 = min(f3, f4);
      end
      if any(f1)
         fs1=[fs1(1) {'=','('} fs1(f2(1)+1:f3-1) {')','.*','t',';'}];
      else
         fs1 = [fs1(1) fs1(f2(1):end)];
      end
      args.equations(i).equation = sprintf('%s', fs1{:});
   end
end
 if args.debug
    analunits('-l',{},args.vars,args.units);
 end
[un1, vars1, iter] = analyseunits(args.equations, args.vars, args.units, args.debug);
if ~args.silent
   fprintf('%d iterations\n', iter);
   undef = sum(isundefined(un1));
   if undef == 0
      fprintf('No problems found\n');
   else
      fprintf('No problems found, but still %d units undefined\n',undef);
   end
   un2 = cell(size(un1));
   for i = 1:length(un1)
      un2{i} = char(un1(i));
   end
   analunits('-l',{},vars1, un2)
end
if nargout > 0
   un = un1;
   vars = vars1;
end
end

