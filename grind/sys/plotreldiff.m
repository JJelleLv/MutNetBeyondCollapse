%PLOTRELDIFF   Plot per capita growth of 1D differential equation
%   Create a plot of the per capita growth as defined by a differential equation of a one dimensional
%   system. The per capita growth/decay is plotted as function of the state variable.
%
%   Usage:
%   PLOTRELDIFF - plots the variable of the x axis with the range of the x-axis.
%   PLOTRELDIFF VAR1 LIM - plots the state variable VAR1 with a range of LIM.
%   PLOTRELDIFF 'parameter',value,...) - Valid parameter-value pairs [type]:
%     'hax' [handle] - handle to axis, used by replayall
%     'npoints' [integer>0] - number of points used for the plot.
%     'var' [state variable] - state variable to plot
%     'xlim' [number and length(number)==2] - the range of state variable.
%   PLOTRELDIFF('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-p' - makes it possible to replay a paranal session in combination with plotdiff.
%
%   See also plotdiff, null, phas, ax
%
%   Reference page in Help browser:
%      <a href="matlab:commands('plotreldiff')">commands plotreldiff</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [aX,aY,hplot] = plotreldiff(varargin)
%(avar, alim, npoints, hax)
global g_grind;
if isfield(g_grind,'xaxis')
    xaxis=g_grind.xaxis;
else
    xaxis=struct('var','','lim',[0 10]);
end
fieldnams={'var', 'v', 'state variable to plot',xaxis.var;...
   'xlim', 'n&length(n)==2', 'the range of state variable.',xaxis.lim;...
   'npoints', 'i>0', 'number of points used for the plot.',500;...
   'hax', 'h', 'handle to axis, used by replayall',[]}';
args=i_parseargs(fieldnams,'if ~isempty(str2num(args{1})),deffields=''xlim,npoints,hax'';else,deffields=''var,xlim,npoints,hax'';end;','-p',varargin);

defnpoints=500;

if ~isfield(args,'npoints')||isempty(args.npoints)
   args.npoints =defnpoints;
end
if ~isfield(args,'var')||isempty(args.var)
   args.var = i_statevars_names(1);
end
if ~isfield(args,'xlim')||isempty(args.xlim)
   if strcmp(g_grind.xaxis.var, args.var)
      args.xlim = g_grind.xaxis.lim;
   elseif strcmp(g_grind.yaxis.var, args.var)
      args.xlim = g_grind.yaxis.lim;
   elseif strcmp(g_grind.zaxis.var, args.var)
      args.xlim = g_grind.zaxis.lim;
   else
      args.xlim = [0 10];
   end
end

i_parcheck;
if ~isfield(args,'hax')
   args.hax = [];
end

if ischar(args.var)&&strncmpi(args.var, '-p', 2)
        [hfig, new] = i_makefig('phase1a');
        if new
            set(hfig,'tag','plotdiff')
            set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
            set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
            hold on;
        end
    plotdiff('-#plotreldiff','hax',args.hax);
    return;
end

[X, Y, args.var] = getplotreldiff(args.npoints, args.xlim, args.var);
if nargout  == 2
   aX = X;
   aY = Y;
else
   if isempty(args.hax)||~ishandle(args.hax)
      [hfig, new] = i_makefig('phase1a');
      args.hax=get(get(0,'currentfigure'),'currentaxes');
      if isempty(args.hax)
          args.hax=gca;
      end
    else
      hfig=get(args.hax,'parent');
      new=0;
    end
  % [H, new] = i_makefig('phase1');
   set(hfig, 'Name','Plot of 1D differential equation')
   ud=get(args.hax,'userdata');
   ud.meta=struct('func','plotreldiff','xname',args.var,'xlim',args.xlim,'yname',sprintf('%s''/%s',args.var,args.var),'zname','');  
   set(args.hax,'userdata',ud);
   if new
      set(hfig,'tag','plotdiff')
      set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
      set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
      hold on;
   end
   apen = g_grind.pen;
   apen.i = 1; %  do not cycle colors between runs
   apen.cycle = true;
   oldhold = ishold;
   if oldhold
      series = get(args.hax, 'children');
      for i = 1:length(series)
         apen = i_nextpen(apen);
      end
   end
   hplot = plot(args.hax,X, Y,apen.pen, 'Color', apen.color);
   i_plotdefaults(hfig)
   hold(args.hax,'on');
   axtext=sprintf('%s''/%s',i_disptext(args.var),i_disptext(args.var));
   i_grindlegend(2,hplot,axtext)
   xlabel(i_disptext(args.var));
   ylabel(axtext);
   if isempty(findobj(hfig,'tag','zeroline'))||~oldhold
      h = plot(args.hax,X, zeros(size(X)),'k');
      set(h,'tag','zeroline');
      i_grindlegend(-1, h);
      if ~oldhold
         hold(args.hax,'on');
      end
   end
   if g_grind.statevars.dim > 1
      N0=i_initvar;
      iX=i_getno(args.var);
      htitle=title(args.hax,i_disptext(['Valid for ' i_othervars(N0, iX.no)]));
      set(htitle,'fontweight','normal');
   end
   set(args.hax, 'XLim', args.xlim);
   if ~oldhold
      hold(args.hax,'off');
   end
   if nargout > 2
      aX = X;
      aY = Y;
   end
   i_grindlegend(11,args.hax);
end
if g_grind.solver.isdiffer
   i_warningdlg('GRIND:plotdiff:differenceeq','This function is designed for differential equations, please use <a href="itermap">itermap</a> instead');
end

function [X, Y, avar] = getplotreldiff(npoints, alim, avar)
global t;
iX = i_varno(avar);
if isempty(iX)
   warning('GRIND:plotreldiff:NoStatevar','Can only create plot if state variable is on the axis, plotted the first state variable');
   avar = i_statevars_names(1);
   iX = i_varno(avar);
end
N0 = i_initvar;
X = transpose((alim(1):(alim(2) - alim(1))  / npoints:alim(2)));
if length(N0)>1
    N0=repmat(transpose(N0),length(X),1);
    N0(:,iX)=X;
else
   N0=X;
end 
Y=i_runsinglestep(t,N0,true);
Y=Y(:,iX);
Y=Y./X;

