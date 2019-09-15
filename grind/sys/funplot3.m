%FUNPLOT3   Plot a 3D surface function
%   Plot some equation in a surface plot (including parameters, state variables, time and
%   functions). For state variables the current initial value is used.
%
%   Usage:
%   FUNPLOT3 - without arguments, a dialog box appears.
%   FUNPLOT3 FUN VARX VARY XRANGE YRANGE - FUN is a function which is 
%   plotted on the y-axis ; VARX is the dependent variable for 
%   the X axis; VARY is the dependent variable for 
%   the Y axis; XRANGE is the range of which the variable VAR1 is 
%   varied. YRANGE sets the range of the Y axis.
%   FUNPLOT3('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - equation to be plotted
%     'hax' [handle] - the handle of the axes
%     'log' [logx | logy | logxy or empty] - makes the scale x-axis or the y-axis logarithmic:
%   values: 'logx' horizontal axis log scale, 'logy' vertical axis 
%   log scale, 'logxy' both log scales
%     'npoints' [integer>0 and length(integer)<=2] - the number of grid points
%     'varx' [string] - the dependent variable for the X axis;
%     'vary' [string] - the dependent variable for the y axis;
%     'xrange' [number and length(number)==2] - the range of the x-axis
%     'yrange' [number and length(number)==2] - the range of the y-axis
%
%   Examples:
%   funplot3 'x^2+y^2' x y [0 10] [0 10]
%  
%
%   See also funcs, implicitplot, funplot
%
%   Reference page in Help browser:
%      <a href="matlab:commands('funplot3')">commands funplot3</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function funplot3(varargin)
%(afun, varx, vary, axrange, ayrange, npoints)
global g_grind;
if ~exist('i_use', 'file')
   addpath([grindpath filesep 'sys2']);
end
fieldnams={'fun', 's', 'equation to be plotted';...
   'varx', 's', 'the dependent variable for the X axis;';...
   'vary', 's', 'the dependent variable for the y axis;';...
   'xrange', 'n&length(n)==2', 'the range of the x-axis';...
   'yrange', 'n&length(n)==2', 'the range of the y-axis';...
   'npoints', 'i>0&length(i)<=2', 'the number of grid points';...
   'log', 'e[logx|logy|logxy]#E', 'makes the scale x-axis or the y-axis logarithmic:';...
   'hax', 'h', 'the handle of the axes'}';
g_args=i_parseargs(fieldnams,...
    'fun,varx,vary,xrange,yrange,npoints', '',varargin);
if ~isfield(g_args,'npoints')
    g_args.npoints = 20;
end
g_defx={'x','y'};
if isfield(g_args,'fun')
    vars = symvar(g_args.fun);
    if (length(vars) == 2)
        if ~isfield(g_args,'varx')
           g_args.varx=vars{1};
        end
        if ~isfield(g_args,'vary')
           g_args.vary=vars{2};
        end
        %assume x=x and y=y
        if strcmp(g_args.varx,'y')&&strcmp(g_args.vary,'x')
            g_args.varx='x';
            g_args.vary='y'; 
        end
    elseif length(vars) > 1
        if ~any(strcmp(vars, g_defx{1}))
            g_defx = vars(1:2);
        end
    end
    clear vars;
else
    g_args.fun='';
end
if ~isfield(g_args,'varx')
    g_args.varx='';
end
if ~isfield(g_args,'vary')
    g_args.vary='';
end
if ~isfield(g_args,'xrange')
    g_args.xrange = [0 100];
end
if ~isfield(g_args,'yrange')
    g_args.yrange = [0 100];
end
if ~isfield(g_args,'log')
    g_args.log = '';
end

if isempty(g_args.fun)||isempty(g_args.varx)||isempty(g_args.vary)
   prompt={'3D function','variable on x axis','variable on y axis','range x','range y','number of raster points'};
   if isfield(g_grind, 'funplot3')
      answer = g_grind.funplot3;
      if ~isempty(g_args.fun)
         answer{1} = g_args.fun;
      end
   else
      if isempty(g_args.varx)
          g_args.varx=g_defx{1};
      end
      if isempty(g_args.vary)
          g_args.vary=g_defx{2};
      end
      answer = {g_args.fun,g_args.varx,g_args.vary,mat2str(g_args.xrange),mat2str(g_args.yrange),num2str(g_args.npoints)};
   end
   answer = inputdlg(prompt, '3D function plot', 1, answer);
   if isempty(answer) || isempty(answer{1})
      error('GRIND:funplot3:NoEquations','no equation entered');
   end
   g_grind.funplot3 = answer;
   g_args1=i_parseargs(fieldnams,'fun,varx,vary,xrange,yrange,npoints','',answer);
   g_args=mergestructs(g_args,g_args1);
end
i=strfind(g_args.fun, '=');
if ~isempty(i)
   g_args.fun=g_args.fun(i(1)+1:length(g_args.fun));
end
implicitplot(g_args.fun,g_args.varx,g_args.vary,g_args.xrange,g_args.yrange,g_args.npoints, 1);
if isfield(g_args,'log')
    if any(strcmp(g_args.log,{'logx','logxy'}))
        set(gca,'xscale','log');
    else
        set(gca,'xscale','linear');
    end
    if any(strcmp(g_args.log,{'logy','logxy'}))
        set(gca,'yscale','log');
    else
        set(gca,'yscale','linear');
    end
end
