%ATTRBASIN   Find basins of attraction by simulation
%   Vary the initial states of two state variables and run till an attractor
%   is reached. Plot the mean state of the attractor. Used to visualize the 
%   sensitive dependence of the initial state of systems with transient chaos. 
%   Use <a href="matlab:help null">null</a> and <a href="matlab:help perturb">perturb</a> or <a href="matlab:help manifolds">manifolds</a> to find a simple separatrix.
%
%   Usage:
%   ATTRBASIN - creates a 20x20 plot with the equilibrium state of the x-axis
%   ATTRBASIN NPOINTS VAR - creates a NPOINTSxNPOINTS plot with the equilibrium state of variable VAR.
%   RES=ATTRBASIN(NPOINTS,VAR) - save the simulation results in the 3 dimensional matrix RES.
%   ATTRBASIN NPOINTS VAR MATRIX - recreate the plot with the result matrix MATRIX.
%   ATTRBASIN('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'matrix' [number] - matrix with end values of runs to plot
%     'npoints' [integer>0 and length(integer)<=2] - number of points for plot (may be [Nx Ny] pair)
%     'var' [state variable] - name of the state variable used for plotting
%
%
%   See also ax, null, perturb, manifolds
%
%   Reference page in Help browser:
%      <a href="matlab:commands('attrbasin')">commands attrbasin</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function EndsVar = attrbasin(varargin)
%(n, var2, EndsX)
global g_Y t g_grind g_t;
fieldnams={'matrix','n','matrix with end values of runs to plot',[];...
   'npoints','i>0&length(i)<=2','number of points for plot (may be [Nx Ny] pair)',[20 20];...
   'var','v','name of the state variable used for plotting',g_grind.xaxis.var}';
args=i_parseargs(fieldnams,'npoints,var,matrix','',varargin);
i_parcheck;
if ~isfield(args,'npoints')
   args.npoints = [20 20];
end
if numel(args.npoints)==1
    args.npoints=[args.npoints args.npoints];
end
if ~isfield(args,'var')
   args.var = g_grind.xaxis.var;
   iRes=i_getno(args.var);
   if ~iRes.isvar
       args.var=i_statevars_names(1);
   end
end
if isfield(args,'matrix')
    siz=size(args.matrix);
    args.npoints=siz(1:2);
end
Xs = zeros(args.npoints);
Ys = zeros(args.npoints);
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if (isnan(iX.no)  ||  isnan(iY.no))
    ax('-?');
    error('GRIND:attrbasin:NoStatevars','Cannot create attrbasin if there are no state variables/parameters on the axes (see <a href="matlab:ax">ax</a>).');
end
oldX = i_getparvar(iX);
oldY = i_getparvar(iY);
try
   iRes = i_getno(args.var);
   if ~isfield(args,'matrix')
      EndsX = zeros(args.npoints(1), args.npoints(2), g_grind.statevars.dim);
      if ~g_grind.solver.isdiffer
         stabstep = 2;
      else
         stabstep = NaN;
      end
      oldstep = g_grind.tstep;
      hwait=i_waitbar(0, args.npoints(1),'attrbasin' ,'Calculating',0.5,true);
      N0 = i_initvar;
      for j = 1:args.npoints(1)
         i_waitbar(1)
         px = g_grind.xaxis.lim(1) + (g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) * j / args.npoints(1);
         N0 = i_setparvar(iX, N0, px);
         for k = 1:args.npoints(2)
            py = g_grind.yaxis.lim(1) + (g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) * k / args.npoints(2);
            N0 = i_setparvar(iY, N0, py);
            Xs( j,k) = px;
            Ys(j,k) = py;
            t1 = t;
            g_grind.tstep = stabstep;
            i_ru(t1, g_grind.ndays, N0, 0);
            N1 = transpose(g_Y(size(g_Y, 1), :));
            t1 = g_t(size(g_t, 1));
            g_grind.tstep = oldstep;
            i_ru(t1, g_grind.ndays / 10, N1, 0);
            mY = mean(g_Y,1);
            for l = 1:g_grind.statevars.dim
               EndsX(j, k, l) = mY(l);
            end
            if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
               disp('Cancelled by user');
               break;
            end
         end
         if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
            break;
         end
      end
      i_waitbar([]);
   else
      EndsX=args.matrix;      
      N0 = i_initvar;
      for j = 1:args.npoints(1)
         px = g_grind.xaxis.lim(1) + (g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) * j / args.npoints(1);
         N0 = i_setparvar(iX, N0, px);
         for k = 1:args.npoints(2)
            py = g_grind.yaxis.lim(1) + (g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) * k / args.npoints(2);
            N0 = i_setparvar(iY, N0, py);
            Xs(j,k) = px;
            Ys(j,k) = py;
         end
      end
   end
   colmap = 'i_grindcolormap';
   hfig = i_makefig('attrbasin');
   set(hfig,'DoubleBuffer','on');
   if ~isoctave&&verLessThan('matlab','8.6')
      set(hfig,'renderer','ZBuffer'); %in fact this is solving a window bug
   end
   set(hfig,'name','Basin of attraction');
   pcolor(Xs, Ys, EndsX(:, :, iRes.no));
   i_plotdefaults(hfig);
   ylabel(['initial ' g_grind.yaxis.var]);
   xlabel(['initial ' g_grind.xaxis.var]);
   shading flat;
   colormap(colmap);
   cb = colorbar;
   i_plotdefaults(hfig);
   ylab = get(cb, 'ylabel');
   set(ylab,'string',['Mean equilibrium ' args.var]);
   if nargout > 0
      EndsVar = EndsX;
   end
   i_plotdefaults(hfig);
   N0 = i_setparvar(iX, N0, oldX);
   i_setparvar(iY, N0, oldY);
catch err
%   err=lasterror;
   N0 = i_setparvar(iX, N0, oldX);
   i_setparvar(iY, N0, oldY);
   rethrow(err);
end
