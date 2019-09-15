%PLOTDIFF   Plot 1D differential equation
%   Create a plot of the growth as defined by a differential equation of a one dimensional
%   system. The growth/decay is plotted as function of the state variable.
%
%   Usage:
%   PLOTDIFF - plots the variable of the x axis with the range of the x-axis.
%   PLOTDIFF VAR XLIM - plots the state variable VAR with a range of XLIM.
%   PLOTDIFF('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'hax' [handle] - handle to axis, used by replayall
%     'npoints' [integer>0] - number of points used for the plot.
%     'var' [state variable] - state variable to plot
%     'xlim' [number and length(number)==2] - the range of state variable.
%   PLOTDIFF('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-#plotreldiff' - used internally.
%     '-c' - not used, for compatibility with null
%     '-m=PAR' - find the Maxwell point by adapting parameter PAR.
%     '-p' - makes it possible to replay a paranal session in combination with plotdiff.
%     '-q' - not used, for compatibility with null
%     '-s' - calculate the areas before and after saddle point (for Maxwell point).
%
%  
%   See also null, phas, ax, plotreldiff, replayall, fokkerplanck
%
%   Reference page in Help browser:
%      <a href="matlab:commands('plotdiff')">commands plotdiff</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [aX,aY,hplot] = plotdiff(varargin)
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
args=i_parseargs(fieldnams,...
  'if ~isempty(str2num(args{1})),deffields=''xlim,npoints,hax'';else,deffields=''var,xlim,npoints,hax'';end;','-p,-m,-s,-c,-q,-#plotreldiff',varargin);
defnpoints = 500;
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
if isfield(g_grind.solver, 'switch_stochast')
   oldpars = cell(size(g_grind.solver.switch_stochast));
   for i = 1:length(g_grind.solver.switch_stochast)
      oldpars{i} = evalin('base', g_grind.solver.switch_stochast{i});
      aval = zeros(size(oldpars{i}));
      assignin('base', g_grind.solver.switch_stochast{i}, aval);
   end
end
if ~isfield(args,'hax')
   args.hax = [];
end
if any(strcmp(args.opts,'-p'))
   updateparanalreplay(args.hax, false, args.var, args.xlim);
   return;
elseif any(strcmp(args.opts,'-#plotreldiff'))  %used internally by plotreldiff
   updateparanalreplay(get(get(0,'currentfigure'),'currentaxes'),true,args.var,args.xlim);
   return;
elseif any(strncmp(args.opts,'-m',2))
   % -m=par
   %find maxwell point
   opt=args.opts{ strncmp(args.opts,'-m',2)};
   f=strfind(opt, '=');
   if isempty(f)
      error('grind:plotdiff','Unknown which parameter to change (use -m=par)')
   else
      parname = opt(f(1) + 1:end);
   end
   disp('Searching a Maxwell point...');
   p0 = evalin('base', parname);
   p = fminsearch(@(p0)findmaxwell(p0, parname), p0);
   fprintf('A Maxwell point is found at %s = %g\n', parname, p);
   plotdiff('-s');
elseif any(strcmp(args.opts,'-s'))
   ran=abs(args.xlim(2)-args.xlim(1));
   [X, Y] = getplotdiff(10000, args.xlim+[-1/10000*ran 1/10000*ran], args.var);
   posY=Y>=0;
   saddles=find(diff(posY)>0);
   nodes=find(diff(posY)<0);
   k=1;
   areas=struct('node1',{},'node2',[],'saddle',[],'area1',[],'area2',[],'diff',[]);
   for i = 1:length(saddles)
       node1=nodes(find(nodes<saddles(i),1,'last'));
       node2=nodes(find(nodes>saddles(i),1,'first'));
       if ~isempty(node1)&&~isempty(node2)
           area1.node1=X(node1);
           area1.node2=X(node2);
           area1.saddle=X(saddles(i));
           area1.area1 = trapz(Y(node1:saddles(i)));
           area1.area2 = trapz(Y(saddles(i):node2));
           area1.diff= abs(area1.area1+area1.area2);
           areas(k)=area1;
           k=k+1;
       end
   end
   aa=[areas.diff];
   if nargout > 0
      aX = areas;
      aY = aa;
   else
      for i = 1:length(areas)
           fprintf('Area between node %s=%g and saddle %s=%g is %g\n', g_grind.xaxis.var, areas(i).node1, g_grind.xaxis.var, areas(i).saddle, areas(i).area1);
           fprintf('Area between saddle %s=%g and node %s=%g is %g\n', g_grind.xaxis.var, areas(i).saddle, g_grind.xaxis.var, areas(i).node2,  areas(i).area2);
           fprintf('Difference = %g\n\n',areas(i).diff);
      end
   end
   return;
end

%if size(g_grind.statevars,2)>1
%   i_errordlg('Warning: This command is designed for 1D differential equations');
%end
[X, Y, arg.var] = getplotdiff(args.npoints, args.xlim, args.var);
if nargout  == 2
   aX = X;
   aY = Y;
else
   if isempty(args.hax)||~ishandle(args.hax)
      [hfig, new] = i_makefig('phase1');
      i_figure(hfig);
      args.hax=get(get(0,'currentfigure'),'currentaxes');
      if isempty(args.hax)
         args.hax = gca;
      end
   else
      hfig = get(args.hax, 'parent');
      new = 0;
   end
   % [H, new] = i_makefig('phase1');
   set(hfig, 'Name','Plot of 1D differential equation')
   if new
      set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
      set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
      set(hfig,'tag','plotdiff')
      hold on;
   end
   ud = get(args.hax, 'userdata');
   ud.var = arg.var;
   ud.meta=struct('func','plotdiff','xname',ud.var,'xlim',args.xlim,'yname',sprintf('%s''',ud.var),'zname','');  
   set(args.hax, 'userdata', ud);
   apen = g_grind.pen;
   apen.i = 1; %  do not cycle colors between runs
   apen.cycle = true;
   oldhold = ishold;
   if oldhold
      series = get(args.hax, 'children');
      for i = 1:length(series)
         apen = nextpen(apen);
      end
   end
   hplot = plot(args.hax,X, Y,apen.pen, 'Color', apen.color);
   i_plotdefaults(hfig)
   hold(args.hax, 'on');
   i_grindlegend(2,hplot,[i_disptext(args.var) ''''])
   xlabel(i_disptext(args.var));
   ylabel([i_disptext(args.var) '''']);
   if isempty(findobj(hfig,'tag','zeroline'))||~oldhold
      h = plot(args.hax,X, zeros(size(X)),'k');
      set(h,'tag','zeroline');
      i_grindlegend(-1, h);
      if ~oldhold
         hold(args.hax, 'on');
      end
   end
   if g_grind.statevars.dim > 1||g_grind.solver.nonautonomous
      N0 = i_initvar;
      iX = i_getno(args.var);
      htitle = title(args.hax,i_disptext(['Valid for ' i_othervars(N0, iX.no)]));
      set(htitle,'fontweight','normal');
   end
   set(args.hax, 'XLim', args.xlim);
   if ~oldhold
      hold(args.hax, 'off');
   end
   if nargout > 2
      aX = X;
      aY = Y;
   end
   i_grindlegend(11, args.hax);
end
if isfield(g_grind.solver, 'switch_stochast')
   for i = 1:length(g_grind.solver.switch_stochast)
      assignin('base', g_grind.solver.switch_stochast{i}, oldpars{i});
   end
end
if g_grind.solver.isdiffer
   i_warningdlg('GRIND:plotdiff:differenceeq','This function is designed for differential equations, please use <a href="itermap">itermap</a> instead');
end
function [X, Y, avar] = getplotdiff(npoints, alim, avar)
global t;
iX = i_varno(avar);
if isempty(iX)
   warning('GRIND:plotdiff:NoStatevar','Can only create plot if state variable is on the axis, plotted the first state variable');
   avar = i_statevars_names(1);
   iX = i_varno(avar);
end
N0 = i_initvar;
X = transpose(alim(1):(alim(2) - alim(1))  / npoints:alim(2));
if length(N0)>1
    N0=repmat(transpose(N0),length(X),1);
    N0(:,iX)=X;
else
   N0=X;
end 
Y=i_runsinglestep(t,N0,true);
Y=Y(:,iX);

function p = findmaxwell(apar, parname)
assignin('base', parname, apar);
[~, adiff] = plotdiff('-s');
if isempty(adiff)
   p = 9999;
else
   p = min(adiff);
end
function ud = updatenullclines(ud, plotted, avar, alim, hax)
global g_paranal ;
i_waitbar(0, size(g_paranal.run.parvalues,1), 'null', 'Preparing replay of paranal',0.5);
oldpar = evalin('base', g_paranal.run.pars{1});
try
   changed = false;
   for i = 1:size(g_paranal.run.parvalues, 1)
      i_waitbar(1);
      if isempty(g_paranal.nulldata.data{i})
         changed = true;
         multassignin('base', g_paranal.run.pars{1}, g_paranal.run.parvalues(i,1));
         [X, Y] = getplotdiff(500, alim, avar);
         g_paranal.nulldata.data{i} = [X, Y];
      end
      if ~plotted
         if ud.replay.isrelative
            [~, ~, ud.replay.hnull] = plotreldiff( avar,alim,[],hax);
         else
            [~, ~, ud.replay.hnull] = plotdiff(avar,alim,[],hax);
         end
         plotted = all(ishandle(ud.replay.hnull));
         if plotted
            set(ud.replay.hnull,'tag','nullp');
         end
      end
   end
   if changed||~isfield(g_paranal.nulldata, 'parvalues')
      g_paranal.nulldata.parvalues = g_paranal.run.parvalues;
   end
   multassignin('base', g_paranal.run.pars{1}, oldpar);
catch err
   multassignin('base', g_paranal.run.pars{1}, oldpar);
   rethrow(err);
end
hold(hax, 'on');
%ndx=g_paranal.run.p == g_paranal.run.p(1);
%ud.replay.HLine = plot(g_paranal.run.Y(ndx, ud.replay.iX.no), g_paranal.run.Y(ndx, ud.replay.iY.no), '-');
%set(ud.replay.HLine,'tag','nullp');
%i_grindlegend(-1, ud.replay.HLine);
i_waitbar([]);

function updateparanalreplay(hax, isrelative, avar, alim)
global g_paranal ;
if nargin < 1
   hax = [];
end
if isempty(hax)||~ishandle(hax)
   H = i_figno('phase1');
   if ishandle(H)
      hax=findobj(H,'type','axes');
      tags = get(hax, 'tag');
      if ~isempty(tags)
         hax = hax(~strcmp(tags, 'legend'));
      elseif isempty(hax)
         hax = gca;
      end
      ud = get(hax, 'userdata');
      hh=findobj(H,'Tag','nullp');
      delete(hh);
      i_grindlegend(12, hax);
   else
      H = i_makefig('phase1');
      hax=findobj(H,'type','axes');
      if isempty(hax)
         hax = gca;
      end
      tags = get(hax, 'tag');
      if ~isempty(tags)
         hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
      end
   end
else
   hh=findobj(hax,'Tag','nullp');
   delete(hh);
end
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
ud.replay.onturn = @i_replayparanalturn;
ud.replay.iX = i_getno(avar);
ud.replay.isrelative = isrelative;
ud.replay.hwhere = -1;
ud.replay.HLine = -1;
tstep = g_paranal.run.parvalues(2, 1)  - g_paranal.run.parvalues(1, 1);
ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
   g_paranal.nulldata.data = cell(size(g_paranal.run.parvalues));
end
ud = updatenullclines(ud, false, avar, alim, hax);
set(hax,'userdata', ud);
replayall('variable', g_paranal.run.pars{1});


function replaystart(hax)
global g_paranal ;
if ishandle(hax)&&~isempty(g_paranal)&&~isempty(g_paranal.run)
   i_figure(get(hax, 'parent'));
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      tstep = g_paranal.run.parvalues(2, 1)  - g_paranal.run.parvalues(1, 1);
      ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
      if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
         %pars = unique(g_paranal.run.p);
         %if ~strcmp(g_paranal.run.parname, ud.replay.settings.tvar)||(length(pars)~=length(ud.replay.pars))||(pars(1)~=ud.replay.pars(1))||pars(end)~=ud.replay.pars(end)
         if ishandle(ud.replay.HLine)
            delete(ud.replay.HLine);
         end
         if ud.replay.isrelative
            plotreldiff('-p')
         else
            plotdiff('-p')
         end
         return;
      end
      if ishandle(ud.replay.HLine(1))
         ax1 = get(ud.replay.HLine(1), 'parent');
         if ax1 ~= hax
            ud.replay.hwhere = -1;
            delete(ud.replay.HLine);
            ud.replay.HLine(1) = -1;
            hh=findobj(hax,'Tag','nullp');
            delete(hh);
            i_grindlegend(12, hax);
            ud = updatenullclines(ud, false, hax);
         end
      end
      if ~ishandle(ud.replay.hnull(1))||(get(ud.replay.hnull(1), 'parent')~=hax)
         hh=findobj(hax,'Tag','nullp');
         delete(hh);
         i_grindlegend(12, hax);
         ud = updatenullclines(ud, false, hax);
      end
      set(hax,'ylimmode','manual');
      ylim1 = get(hax, 'ylim');
      set(hax, 'ylim', [ - max(abs(ylim1)), max(abs(ylim1))])
      set(hax,'zlimmode','manual');
      oldnext = get(hax, 'NextPlot');
      set(hax,'NextPlot','add');
      % ud.replay.iX = i_getno(g_grind.xaxis.var);
      if ~ishandle(ud.replay.HLine)
         ud.replay.HLine = plot(hax,g_paranal.run.Y(1,ud.replay.iX.no,1), 0, '-');
         set(ud.replay.HLine,'tag','nullp');
      end
      if ~ishandle(ud.replay.hwhere)
         ud.replay.hwhere = plot(hax,g_paranal.run.Y(1,ud.replay.iX.no,1), 0, 'ro');
         set(ud.replay.hwhere, 'MarkerFaceColor', [1 0 0]);
         set(ud.replay.hwhere,'tag','nullp');
         i_grindlegend(-1, ud.replay.hwhere);
      end
      set(hax, 'NextPlot', oldnext);
   end
   set(hax, 'userdata', ud);
end

function replayend(hax, closedlg)
try
   if closedlg
      replaycallback(hax,'', 1);
   end
catch
end
if ishandle(hax)
   i_figure(get(hax, 'parent'));
   ud = get(hax, 'userdata');
   if ishandle(ud.replay.hwhere)
      delete(ud.replay.hwhere);
      ud.replay.hwhere = -1;
   end
   set(hax, 'userdata', ud);
end

function p = replaycallback(hax, avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&~(isempty(g_paranal)||isempty(g_paranal.run)||isempty(avar)||~strcmp(avar, g_paranal.run.pars{1}))
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
      [tndx, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
      p = g_paranal.run.parvalues(stepndx, 1);
      [~, pnow] = min(abs(g_paranal.nulldata.parvalues(:,1) - p));
      cdata = g_paranal.nulldata.data{pnow};
      xdata = g_paranal.run.Y(1:tndx, ud.replay.iX.no,stepndx);
      if ud.replay.isrelative
         set(ud.replay.hnull, 'XData',cdata(:,1) ,'YData', cdata(:,2)./cdata(:,1));
      else
         set(ud.replay.hnull, 'XData',cdata(:,1) ,'YData', cdata(:,2));
      end
      if ~isempty(xdata)
         if ishandle(ud.replay.hwhere)
            set(ud.replay.hwhere,'Xdata',xdata(end,1),'Ydata',0);
         end
         if ishandle(ud.replay.HLine)
            set(ud.replay.HLine,'Xdata',xdata(:,1),'Ydata',zeros(size(xdata)));
         end
      end
   end
end


