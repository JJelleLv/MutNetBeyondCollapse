%RETURNTIME2D   Plot of the time necessary to reach a stable node
%   Create a 2D contour plot of the phase plane containing the time to reach any stable equilibrium.
%  
%   Usage:
%   RETURNTIME2d - estimates the return time based on the current <a href="matlab:help simtime">simtime</a> settings and a default
%   value for the maximum change in equilibrium (1E-8).
%   [V,X,Y]=RETURNTIME2D write the results also to V X and Y for further use for instance in SURF or
%   other MATLAB functions.
%   RETURNTIME2D ERR - use the maximum change of ERR.
%   RETURNTIME2D ERR MAXT - use ERR and simulate MAXT time units to find an equilibrium.
%   RETURNTIME2D('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'abstol' [number>0] - if the distance to the equlibrium < ABSTOL, the model is returned to equilibrium.
%     'method' [change | equil | both] - method to determine if equilibrium is reached 'both', 'change' or 'equil'
%     'ndays' [number>0] - number of days for run.
%     'npoints' [integer>0 and length(integer)<=2] - number of points for the contour plot
%  
%      
%   See also returntime, ru, simtime, findeq
%
%   Reference page in Help browser:
%      <a href="matlab:commands('returntime2d')">commands returntime2d</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [Vect, X, Y, Vect3d] = returntime2d(varargin)
%(err, npoints)
global g_grind;
fieldnams={'npoints', 'i>0&length(i)<=2', 'number of points for the contour plot',[50 50];...
   'abstol', 'n>0', 'if the distance to the equlibrium &lt; ABSTOL, the model is returned to equilibrium.',1E-8;...
   'ndays', 'n>0', 'number of days for run.',g_grind.ndays;...
   'method','e[change|equil|both]', 'method to determine if equilibrium is reached ''both'', ''change'' or ''equil''','both'}';
args=i_parseargs(fieldnams,'abstol,npoints,ndays,method','',varargin);
i_parcheck;
if ~isfield(args,'abstol')
   args.abstol = 1E-8;
end
if  ~isfield(args,'ndays')
   args.ndays = g_grind.ndays;
end
if  ~isfield(args,'npoints')
   args.npoints = [50, 50];
end
if numel(args.npoints)==1
    args.npoints=[args.npoints,args.npoints];
end
if ~isfield(args,'method')
    args.method='both';
end
if isfield(args,'ndays')
   args.ndays=g_grind.ndays;
end
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if (isempty(iX.no) || isempty(iY.no))
   ax -?;
   i_errordlg('Cannot returntime plot if there are no state variables/parameters on the axes.');
   error('GRIND:returntime2d:InvalidAxes','Cannot returntime plot if there are no state variables/parameters on the axes');
end
save3d=nargout>3;
oldX = evalin('base', g_grind.xaxis.var);
oldY = evalin('base', g_grind.yaxis.var);
try
   Vect = zeros(args.npoints(1), args.npoints(2));
   if save3d
       Vect3d=cell(args.npoints(1), args.npoints(2));
   end
   X = zeros(args.npoints(1), args.npoints(2));
   Y = zeros(args.npoints(1), args.npoints(2));
   %   t = 0;
   minX = g_grind.xaxis.lim(1);
   if minX < 0.0001, minX = 0.0001; end
   maxX = g_grind.xaxis.lim(2);
   minY = g_grind.yaxis.lim(1);
   if minY < 0.0001, minY = 0.0001; end
   maxY = g_grind.yaxis.lim(2);
   incrY = (maxY - minY) / (args.npoints(1)-1);
   incrX = (maxX - minX) / (args.npoints(2)-1);
   %   isdiffer = g_grind.solver.isdiffer;
   i_waitbar(0, args.npoints(1), 'returntime2d', 'Calculating',0.5)
  for y1 = 1:args.npoints(1)
      i_waitbar(1);
      py = (y1 - 1) * incrY + minY;
      evalin('base',sprintf('%s=%g;',g_grind.yaxis.var, py));
      %      assignin('base',g_grind.yaxis.var,py)
      for x1 = 1:args.npoints(2)
         px = (x1 - 1) * incrX + minX;
         evalin('base',sprintf('%s=%g;',g_grind.xaxis.var, px));
         %        assignin('base', g_grind.xaxis.var, px)
        
         [Vect(x1, y1),yy] = returntime('abstol',args.abstol,'ndays',args.ndays,'method',args.method);
         if save3d
            Vect3d{x1, y1}= yy; 
         end
         % Vect2=total length of the trajectory - is maybe closer to the potential
         X(x1, y1) = px;
         Y(x1, y1) = py;
      end
   end
   i_waitbar([]);
   evalin('base',sprintf('%s=%g;',g_grind.xaxis.var, oldX));
   evalin('base',sprintf('%s=%g;',g_grind.yaxis.var, oldY));
   %   assignin('base', g_grind.xaxis.var, oldX)
   %   assignin('base', g_grind.yaxis.var, oldY)
catch err
   %   err=lasterror;
   evalin('base',sprintf('%s=%g;',g_grind.xaxis.var, oldX));
   evalin('base',sprintf('%s=%g;',g_grind.yaxis.var, oldY));
   %assignin('base', g_grind.xaxis.var, oldX)
   %assignin('base', g_grind.yaxis.var, oldY)
   rethrow(err);
end
[H, new] = i_makefig('phase2');
if new
    set(H, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
    set(H, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
end
set(H, 'Name', 'Phase plane');
oldhold = ishold;
hold on;
surf(X, Y, Vect);
shading flat;
set(gca,'XLim', g_grind.xaxis.lim);
set(gca,'YLim', g_grind.yaxis.lim);
if ~oldhold
   hold off;
end



