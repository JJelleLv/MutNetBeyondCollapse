%VECTOR   direction vectors
%   Show the direction of change as vectors in the 2D phase plane.
%
%   Usage:
%   VECTOR - shows 20x20 arrows.
%   VECTOR N - shows NxN arrows.
%   VECTOR N -contour - makes a contour plot of the speed of change.
%   VECTOR N -contour -relative - makes a contour plot of the relative speed of change.
%   VECTOR('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'contour' [logical] - make contour plot
%     'method' [integer] - method=0 - normal; %method=1 - relative method=2 - accel (not working) 
%    method=4 - quasi steady state
%     'npoints' [integer>0 and length(integer)<=2] - size of the grid for the arrows
%   VECTOR('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - make contour plot instead of arrows
%     '-r' - plot relative change (change/state)
%
%   See also phas, null
%
%   Reference page in Help browser:
%      <a href="matlab:commands('vector')">commands vector</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [X1, Y1, Vect1]=vector(varargin)
%(igrid, opt, opt2)
fieldnams={'npoints', 'i>0&length(i)<=2', 'size of the grid for the arrows',20;...
   'contour', 'l', 'make contour plot',false;...
   'method', 'i', 'method=0 - normal; %method=1 - relative method=2 - accel (not working)',0}';
args=i_parseargs(fieldnams,'npoints','-c,-r',varargin);
i_parcheck;
if any(strcmp(args.opts,'-c'))
    args.contour=true;
end
if ~isfield(args,'contour')
    args.contour=false;
end
if ~isfield(args,'npoints')
    if args.contour
        args.npoints=100;
    else
        args.npoints=20;
    end
end
if any(strcmp(args.opts,'-r'))
    args.method=1;
end
if ~isfield(args,'method')
    args.method=0;
end

global g_grind;
iX = i_getno(g_grind.xaxis.var,true);
iY = i_getno(g_grind.yaxis.var,true);
if ~(iX.isvar||iX.ispar||iX.isext) || ~(iY.isvar||iY.ispar||iX.isext)
   ax('-?');
   i_errordlg('Cannot create vector field if there are no state variables on the axes.');
   error('GRIND:vector:NoStatevars','Cannot create vector field if there are no state variables on the axes');
end
%if ~(isempty(iX.no) || isempty(iY.no))
%method=2
[X, Y, Vect] = i_vector(args.npoints+1, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, [],args.method);
[hfig, new] = i_makefig('phase2');
%  for i=1:size(Vect,3);
%     Vect(:,:,i)=ln(Vect(:,:,i));
%  end
if new
  set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
  set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
  hold on;
end
%plotedit('off');
set(hfig, 'Name', 'Phase plane');
oldhold = ishold;
hold on;
if args.contour
   speed = Vect(:, :, 1).^2;
   for i = 2:size(Vect, 3)
      speed = speed + Vect(:, :, i).^2;
   end
   contourf(X,Y,real(log(sqrt(speed))),20)
   h=gca;
 %  shading flat;
   h1=colorbar;
   ylab = get(h1,'ylabel');
   set(ylab, 'FontAngle',  get(h, 'FontAngle'), ...
          'FontName',   get(h, 'FontName'), ...
          'FontSize',   get(h, 'FontSize'), ...
          'FontWeight', get(h, 'FontWeight'), ...
          'string',     'ln speed of change');
 %  axes(h);
else
   if iX.isvar && iY.isvar
      h = myquiver(X, Y, Vect(:, :, iX.no), Vect(:, :, iY.no));
   elseif iY.isvar
      h = myquiver(X, Y, zeros(size(Vect, 1), size(Vect, 2)), Vect(:, :, iY.no));
   elseif iX.isvar
      h = myquiver(X, Y, Vect(:, :, iX.no), zeros(size(Vect, 1), size(Vect, 2)));
      %else everything is zero
   end
   set(h, 'Color', [0.5 0.5 0.5]);
   h=gca;
end
i_plotdefaults(hfig);
set(h, 'XLim', g_grind.xaxis.lim);
set(h, 'YLim', g_grind.yaxis.lim);
if g_grind.statevars.dim > 2
   htitle=title(['Valid for ' i_othervars(i_initvar, iX.no, iY.no)]);
   set(htitle,'fontweight','normal');
end
if ~oldhold
   hold off;
end
if nargout>0
    X1=X;
    Y1=Y;
    Vect1=Vect;
end
%plotedit on;
%end
function  [h] = myquiver(X, Y, Vx, Vy)
maxlen = -1;
minlen = 999999;
for i = 1:size(Vx, 1)
   for j = 1:size(Vx, 2)
      len = norm([Vx(i, j), Vy(i, j)]);
      if len < minlen
         minlen = len;
      end
      if len > maxlen
         maxlen = len;
      end
   end
end
offset = (maxlen - minlen) / 10;

for i = 1:size(Vx, 1)
   for j = 1:size(Vx, 2)
      len = norm([Vx(i, j), Vy(i, j)]);
      Vx(i, j) = Vx(i, j) * (offset + len) / len;
      Vy(i, j) = Vy(i, j) * (offset + len) / len;
   end
end
[h] = quiver(X, Y, Vx, Vy);
%[Sx,Sy]=meshgrid(1:2:10,1:2:10);
%[h] = streamline(Y, X, Vy, Vx,Sx(:),Sy(:));

