%ITERMAP   Iteration map for 1-D difference equation
%   Plot x(t) versus x(t+n) for a 1-D difference equation (default
%   n=1).
%
%   Usage: 
%   ITERMAP - plot x(t) versus x(t+1). This is simply a plot of
%   the right-hand-sides of the difference equation (with the
%   diagonal).
%   ITERMAP NITERS - plot x(t) versus x(t+N) of the state variable
%   on the x-axis (see <a href="matlab:help ax">ax</a>). NITERS can also be
%   negative for backwards simulations.
%   ITERMAP [NITERS1 NITERS2 .. NITERSn] - you can plot more than one itermaps at a time
%   by replacing the scalar NITERS by a vector with the N's.
%   ITERMAP NITERS VAR XLIM - plot VAR(t) versus VAR(t+NITERS) with a range
%   ITERMAP('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'hax' [handle] - handle for the axis (default = [])
%     'niters' [integer>0] - number of iterations of the map (default = 1)
%     'npoints' [integer>0] - number of points for plotting (default=200)
%     'var' [state variable] - variable to plot on the x-axis (default see ax)
%     'xlim' [number and length(number)==2] - limits of the x-axis (default: see ax)
%   ITERMAP('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - not used, added for compatibility with null
%     '-p' - makes it possible to replay a paranal in combination with itermap.
%     '-q' - not used, added for compatibility with null
%  
%   See also ru, null, ax, replayall
%
%   Reference page in Help browser:
%      <a href="matlab:commands('itermap')">commands itermap</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function hline = itermap(varargin)
%(n, avar, alim, hax)
global g_grind;
fieldnams={'niters', 'i>0', 'number of iterations of the map (default = 1)',g_grind.solver.iters;...
   'var', 'v', 'variable to plot on the x-axis (default see ax)',g_grind.xaxis.var;...
   'xlim', 'n&length(n)==2', 'limits of the x-axis (default: see ax)',g_grind.xaxis.lim;...
   'hax', 'h', 'handle for the axis (default = [])',[];...
   'npoints', 'i>0', 'number of points for plotting (default=200)',200}';
args=i_parseargs(fieldnams,...
  'if ~argtype(2,''v''),deffields=''niters,xlim,hax'';else,deffields=''niters,var,xlim,hax'';end;','-p,-c,-q',varargin);

i_parcheck;
if ~isfield(args,'niters')||isempty(args.niters)
   args.niters = g_grind.solver.iters;
end
if ~isfield(args,'npoints')||isempty(args.npoints)
    args.npoints=200;
end

if ~isfield(args,'hax')
    args.hax=[];
end
if ischar(args.hax)
    args.hax=i_checkstr(args.hax);
end
if any(strcmp(args.opts, '-p'))
   updateparanalreplay(args.niters, args.hax);
   return;
end

if ~isfield(args,'var')
   args.var = i_statevars_names(1);
end
if ~isfield(args,'xlim')
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
if ~(numel(args.niters) > 1)
    g_grind.solver.iters  = args.niters;
else
    iters=args.niters;
    for i=1:length(iters)
        args.niters = iters(i);
        itermap(args);
    end
    return;
end
if ~g_grind.solver.isdiffer
   err = i_warningdlg('GRIND:itermap:nodifference','This function is designed for difference equations, please use <a href="matlab:plotdiff">plotdiff</a> instead');
else
   err = -1;
end
if length(args.niters) > 1
   for i = 1:length(args.niters)
      itermap(args.niters(i), args.var, args.xlim, args.hax);
   end
   return;
end
[X, Y] = getitermap(args.niters, args.var, args.xlim, args.npoints);
if isempty(args.hax)||~ishandle(args.hax)
   [hfig, new] = i_makefig('phase1');
   args.hax=get(get(0,'currentfigure'),'currentaxes');
   if isempty(args.hax)
      args.hax = gca;
   end
else
   hfig = get(args.hax, 'parent');
   new = 0;
end
if new
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

if isempty(findobj(hfig,'tag','diagonal'))||~oldhold
   h = plot(args.hax,args.xlim, args.xlim, 'k');
   set(h,'tag','diagonal');
   i_grindlegend(-1, h); % Exclude line in legend
   i_plotdefaults(hfig);
   if ~oldhold
      hold on;
   end
end
if args.niters > 0
   h = plot(args.hax,X, Y, 'Color', apen.color2);
elseif args.niters < 0
   h = plot(args.hax,Y, X, 'Color', apen.color2);
else
   h = plot(args.hax,X,X, 'Color',apen.color2);
end
ylab = i_disptext(sprintf('%s_{t%+d}', char(args.var), args.niters));
i_grindlegend(2, h, ylab) %include in legend
if g_grind.statevars.dim > 1
   N0 = i_initvar;
   iX = i_getno(args.var);
   htitle = title(args.hax,i_disptext(['Valid for ' i_othervars(N0, iX.no)]));
   set(htitle,'fontweight','normal');
end
xlabel(args.hax,i_disptext( sprintf('%s_{t}', char(args.var))));
ylabel(args.hax, ylab);
set(args.hax, 'XLim', args.xlim);
i_grindlegend(11);
ud = get(args.hax, 'userData');
ud.var = args.var;
ud.meta=struct('func','itermap','xname',ud.var,'xlim',args.xlim,'yname',sprintf('%s(t+%d)',ud.var,args.niters),'zname','');  
set(args.hax, 'userData', ud);
if ~oldhold
   hold off
end
if ishandle(err)
   i_figure(err)
end
if nargout > 0
   hline = h;
end
function [X, Y] =  getitermap(~, avar, alim ,npoints)
global t;
N0 = i_initvar;
iX = i_getno(avar);
if ~iX.isvar
   i_errordlg('Need to have a state variable on the x-axis');
   error('GRIND:itermap:NoStatevars','Need to have a state variable on the x-axis');
end
X = transpose((alim(1):(alim(2) - alim(1)) / npoints:alim(2)));
if length(N0)>1
    N0=repmat(transpose(N0),length(X),1);
    N0(:,iX.no)=X;
else
   N0=X;
end  
Y=i_runsinglestep(t,N0,true);
Y=Y(:,iX.no);

function ud = updatenullclines(ud, plotted,hax)
global g_paranal g_grind;
g_grind.solver.iters = ud.replay.iters;
i_waitbar(0, size(g_paranal.run.parvalues,1), 'null', 'Preparing replay of paranal',0.5);
oldpar = evalin('base', g_paranal.run.pars{1});
try
   changed = false;
   hascomplex=false;
   for i = 1:size(g_paranal.run.parvalues, 1)
      i_waitbar(1);
      if isempty(g_paranal.nulldata)||length(g_paranal.nulldata.data)>=i&&(isempty(g_paranal.nulldata.data{i}))
         changed = true;
         multassignin('base', g_paranal.run.pars{1}, g_paranal.run.parvalues(i,1));
         [X, Y] = getitermap(g_grind.solver.iters, g_grind.xaxis.var, g_grind.xaxis.lim, 200);
         if any(~isreal(Y))||any(~isreal(X))
             hascomplex=true;
             X=real(X);
             Y=real(Y);
         end
         g_paranal.nulldata.data{i} = [X, Y];
      end
      if ~plotted
         ud.replay.hnull = itermap(g_grind.solver.iters, g_grind.xaxis.var, g_grind.xaxis.lim,hax);
         plotted = all(ishandle(ud.replay.hnull));
         if plotted
            set(ud.replay.hnull,'tag','nullp');
         end
      end
   end
   if changed
      g_paranal.nulldata.parvalues = g_paranal.run.parvalues;
   end
   if hascomplex
      warning('grind:itermap:complex','Imaginary parts of complex values ignored');
   end
   multassignin('base', g_paranal.run.pars{1}, oldpar);
catch err
   multassignin('base', g_paranal.run.pars{1}, oldpar);
   rethrow(err);
end
set(hax,'NextPlot','add')
%ndx=g_paranal.run.p == g_paranal.run.p(1);
%ud.replay.HLine = plot(g_paranal.run.Y(ndx, ud.replay.iX.no), g_paranal.run.Y(ndx, ud.replay.iY.no), '-');
%set(ud.replay.HLine,'tag','nullp');
%i_grindlegend(-1, ud.replay.HLine);
i_waitbar([]);

function updateparanalreplay(iters, hax)
global g_paranal g_grind;
if nargin == 0
   iters = g_grind.solver.iters;
end
if nargin < 2
   hax = [];
end
if isempty(g_paranal.run)||isempty(g_paranal.run.Y)
   error('grind:null:paranal','Paranal should be run before using this option');
end
if isempty(hax)||~ishandle(hax)
   H = i_figno('phase1');
   if ishandle(H)
      hax=findobj(H,'type','axes');
      tags = get(hax, 'tag');
      hax = hax(~strcmp(tags, 'legend'));
      ud = get(hax, 'userdata');
      hh=findobj(H,'Tag','nullp');
      delete(hh);
      i_grindlegend(12, hax);
   else
      H = i_makefig('phase1');
      hax=findobj(H,'type','axes');
      if isempty(hax)
         hax = gca;
         i_plotdefaults(H);
      end
      tags = get(hax, 'tag');
      hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
   end
else
   hh=findobj(hax,'Tag','nullp');
   delete(hh);
end
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
ud.replay.onturn = @i_replayparanalturn;
ud.replay.iX = i_getno(g_grind.xaxis.var);
ud.replay.hwhere = -1;
ud.replay.HLine = -1;
%ud.replay.pars = g_paranal.run.parvalues;
tstep = g_paranal.run.parvalues(2, 1) - g_paranal.run.parvalues(1, 1);
ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
ud.replay.iters = iters;
if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
   g_paranal.nulldata.data = cell(size(g_paranal.run.parvalues));
end
i_waitbar(0, size(g_paranal.run.parvalues,1), 'null', 'Preparing replay of paranal',0.5);
ud = updatenullclines(ud, false, hax);
set(hax,'userdata', ud);
i_waitbar([]);
replayall('variable', g_paranal.run.pars{1});


function replaystart(hax)
global g_paranal g_grind;
if ishandle(hax)
   i_figure(get(hax, 'parent'));
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      %   pars = unique(g_paranal.run.p);
      tstep = g_paranal.run.parvalues(2, 1) - g_paranal.run.parvalues(1, 1);
      ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',...
         tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
      if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
         if ishandle(ud.replay.HLine)
            delete(ud.replay.HLine);
         end
         itermap('-p')
         return;
      end
      %       if ~isfield(ud.replay.settings,'iters')
      %           ud.replay.settings.iters=g_grind.solver.iters;
      %       end
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
      set(hax,'zlimmode','manual');
      oldnext = get(hax, 'NextPlot');
      set(hax,'NextPlot','add');
      ud.replay.iX = i_getno(g_grind.xaxis.var);
      if ~ishandle(ud.replay.HLine)
         ud.replay.HLine = plot(hax,g_paranal.run.Y(1,ud.replay.iX.no,1), g_paranal.run.Y(2,ud.replay.iX.no,1), 'b-');
         set(ud.replay.HLine,'tag','nullp');
      end
      if ~ishandle(ud.replay.hwhere)
         ud.replay.hwhere = plot(hax,g_paranal.run.Y(1,ud.replay.iX.no,1), g_paranal.run.Y(2,ud.replay.iX.no,1), 'ro');
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
if ishandle(hax)&&isempty(avar)||strcmp(avar, g_paranal.run.pars{1})
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
      [tndx, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
      p = g_paranal.run.parvalues(stepndx, 1);
      [~, pnow] = min(abs(g_paranal.nulldata.parvalues(:,1) - p));
      %   ndx=(p == g_paranal.run.p) & (g_paranal.run.t<=t);
      cdata = g_paranal.nulldata.data{pnow};
      set(ud.replay.hnull, 'XData',cdata(:,1) ,'YData', cdata(:,2));
      xdata = g_paranal.run.Y(1:tndx, ud.replay.iX.no,stepndx);
      cobweb = i_getcobweb(xdata);
      if ~isempty(cobweb)
         if ishandle(ud.replay.hwhere)
            set(ud.replay.hwhere,'Xdata',cobweb(end,1),'Ydata',cobweb(end,2));
         end
         if ishandle(ud.replay.HLine)
            set(ud.replay.HLine,'Xdata',cobweb(2:end,1),'Ydata',cobweb(2:end,2));
         end
      end
   end
end
