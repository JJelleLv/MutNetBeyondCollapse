%IMPLICITPLOT   Implicit plot of an equation
%   Create a two-dimensional implicit plot of an equation, which may
%   include parameters and state variables (no 'functions'). The
%   function is drawn in the same plot that <a href="matlab:help funplot">funplot</a> uses.
%  
%   Usage:
%   IMPLICITPLOT FUN VARX VARY - create an implicit plot of 
%   VARX versus VARY of FUN which is a function of both VARX and
%   VARY.
%   IMPLICITPLOT FUN VARX VARY XRANGE YRANGE - a range for the
%   X-axis and/or the Y-axis may be supplied.
%   IMPLICITPLOT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - implicit function to plot
%     'npoints' [integer>0 and length(integer)<=2] - number of raster points
%     'surface' [logical] - make a surface instead of a zero line (used internally by funplot3)
%     'varx' [identifier] - variable for the x axis
%     'vary' [identifier] - variable for the y axis
%     'xrange' [number and length(number)==2] - range of the x axis
%     'yrange' [number and length(number)==2] - range of the y axis
%  
%   Examples:
%   IMPLICITPLOT X^2+Y^2=1 X Y [-1 1] - this function draws a circle.
%   IMPLICITPLOT A*r*(1-A/K)-g*Z*(A/(A+h)) A Z [0.001 10] [0.001 10] draws one of 
%   the nullclines of an algae-zooplankton model.
%  
%   See also funplot, funplot3
%
%   Reference page in Help browser:
%      <a href="matlab:commands('implicitplot')">commands implicitplot</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function implicitplot(varargin)
%(afun, varx, vary, axrange, ayrange, npoints, issurf)
global g_grind;
if ~exist('i_use', 'file')
   addpath([grindpath filesep 'sys2']);
end
fieldnams={'fun', 's', 'implicit function to plot','';...
   'varx', 'U1', 'variable for the x axis','x';...
   'vary', 'U1', 'variable for the y axis','y';...
   'xrange', 'n&length(n)==2', 'range of the x axis',[0 100];...
   'yrange', 'n&length(n)==2', 'range of the y axis',[0 100];...
   'npoints', 'i>0&length(i)<=2', 'number of raster points',100;...
   'surface', 'l', 'make a surface instead of a zero line (used internally by funplot3)',false}';
args=i_parseargs(fieldnams,'fun,varx,vary,xrange,yrange,npoints,surface','',varargin,false,{@i_isid});
if ~isfield(args,'npoints')
  args.npoints = 100;
end
if ~isfield(args,'surface')
   args.surface = false; %used internally by funplot3
end
if ~isfield(args,'fun')||~isfield(args,'varx')||~isfield(args,'vary')
   if ~isfield(args,'fun')
      args.fun = '';
   end
   prompt={'Implicit function','variable on x axis','variable on y axis','range x','range y','number of raster points'};
   if isfield(g_grind, 'implicitplot')
      answer = g_grind.implicitplot;
      if isempty(answer)
         error('GRIND:implicitplot:cancelled','Cancelled');
      elseif ~isempty(args.fun)
         answer{1} = args.fun;
      end
   else
      answer = {args.fun,'x','y','[0 100]','[0 100]','100'};
   end
   answer = inputdlg(prompt, 'Implicit function plot', 1, answer);
   if isempty(answer) || isempty(answer{1})
      error('GRIND:implicitplot:NoEquation','No equation entered');
   elseif isempty(answer{2})
      args.varx = symvar(answer{1});
      args.varx = args.varx{1};
   else
      args.varx = answer{2};
   end
   if isempty(answer{3})
      args.vary = symvar(answer{1});
      args.vary = args.vary{2};
   else
      args.vary = answer{3};
   end
   g_grind.implicitplot = answer;
   args.fun = answer{1};
   args.xrange = i_checkstr(answer{4});
   if isempty(args.xrange)
      args.xrange = [0 100];
   end
   args.yrange = i_checkstr(answer{5});
   args.npoints =  i_checkstr(answer{6});
else
   if ~isfield(args,'xrange')
      args.xrange = [0 100]; %used internally by funplot3
   end
   if ~isfield(args,'yrange')
      args.yrange = [0 100]; %used internally by funplot3
   end
end
fun = args.fun;
i=strfind(fun, '=');
if ~isempty(i)
   fun=[fun(1:i(1)-1) '- (' fun(i(1)+1:length(fun)) ')'];
end
if isempty(args.yrange)
   args.yrange = args.xrange;
end
eval(sprintf('global %s;', args.vary));
g_loc.comm=sprintf('%s=g_loc.y;', args.vary);
if numel(args.npoints)==1
    args.npoints=[args.npoints args.npoints];
end
g_l_ran = args.yrange(1):(args.yrange(2) - args.yrange(1))  / args.npoints(2):args.yrange(2);
eval(sprintf('%s=%g;', args.vary, g_l_ran(1)));
g_V = zeros(args.npoints(1) + 1,args.npoints(2)+1);
g_X = g_V;
g_Y = g_V;
for i = 1:args.npoints(2) + 1
   %  g_loc.y=g_l_ran(i);
   eval(sprintf('%s=%g;', args.vary, g_l_ran(i)));
   %   eval(sprintf('%s=%g;',vary,g_l_ran{i}));
   %   g_loc.comm=sprintf('global %s;\n%s=%g;',vary,g_l_ran(i));
   [xx, yy, g_s_par] = i_funcalc(fun, args.varx, args.xrange, args.npoints(1)); %, g_loc);
   g_V(:, i) = transpose(yy);
   g_X(:, i) = transpose(xx);
   g_Y(:, i) = g_l_ran(i);
end
[hfig, new] = i_makefig('funplot');
if new
   hold('on');
   nser = 0;
else
   nser = length(get(gca, 'children'));
end
settings = [];
for i = 1:length(g_s_par)
   settings.(g_s_par{i}) =  evalin('base', char(g_s_par{i}));
end
if isfield(g_grind, 'pen') && ~isempty(g_grind.pen)
   %  co = g_grind.pen.color2;
   %  pe = g_grind.pen.pen;
   pen = g_grind.pen;
else
   % co = [0 0 1];
   %  pe = '-';
   pen = i_nextpen([]);
end
pen.cycle = true;
for i = 1:nser
   pen = i_nextpen(pen);
end
plotedit('off');
set(hfig,'Name','Function plot');
oldhold = ishold;
hold('on')
if args.surface
   surf(g_X, g_Y, g_V);
   shading('flat');
   set(gca, 'View', [322.5 , 30]);
else
   [~, h] = contour(g_X, g_Y, g_V, [0 0], pen.pen, 'color', pen.color);
   i_grindlegend(1, h, args.fun, settings);
end
if oldhold
   hold('on');
else
   hold('off');
end
i_plotdefaults(hfig);
xlabel(args.varx);
ylabel(args.vary);
htitle = title(args.fun);
set(htitle,'fontweight','normal');
if ~args.surface
   i_grindlegend(11);
end
set(gca, 'YLim', args.yrange);
set(gca, 'XLim', args.xrange);
