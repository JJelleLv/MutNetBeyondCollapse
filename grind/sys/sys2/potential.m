%POTENTIAL   Create potentials
%   Another way to visualize the dynamics of a first-order system is 
%   based on the physical idea of potential energy (Strogatz, 1994).
%   The system moves downhill, equilibria are in a minimum. The landscape 
%   can be generated by solving the potential function of the model:
%   dx/dt=f(x)
%   transform this model to the potential equation:
%   -dV/dx=f(x)
%   The simulation of this model (initial condition=0) gives the potential plot.
%   The current state is often represented by a marble, but one should imagine
%   that this marble is heavily damped and sliding through a greasy substance 
%   (Strogatz, 1994).
%   The ball denotes the current initial state of the system. You can use <a href="matlab:help replayall">replayall</a>
%   to animate the last run.
%  
%
%   Usage:
%   POTENTIAL - creates potentials with full range and 50 points
%   POTENTIAL RANGE NPOINTS ONLY1D - sets range of Y-axis to RANGE and
%   uses NPOINTS for the calculations. If ONLY1D=1 then only one 
%   state variable is varied.
%   POTENTIAL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'hax' [handle] - handle of the axis (used by replayall)
%     'npoints' [integer>0] - number of points for the plot
%     'range' [number and length(number)<=2] - range of the state variable (x-axis)
%     'rangeV' [number] - range of potential (y-axis)
%     'value' [number] - value of the ball in the potential
%   POTENTIAL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-p' - make ready to replay the last <a href="matlab:help replayall">paranal</a> plot.
%     '-u' VALUE - update the position of the ball to VALUE
%
%   See also replayall, marbleplot, null, fokkerplanck, paranal, quasipot
%   
%
%   Reference page in Help browser:
%      <a href="matlab:commands('potential')">commands potential</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [xs,pot] = potential(varargin)
%(VRange, npoints, only2D)
global g_grind;
fieldnams={'rangeV', 'n', 'range of potential (y-axis)',[ -1E30, 1E30];...
   'npoints', 'i>0', 'number of points for the plot',50;...
   'range', 'n&length(n)<=2', 'range of the state variable (x-axis)',g_grind.xaxis.lim;...
   'hax', 'h', 'handle of the axis (used by replayall)',[];...
   'value', 'n', 'value of the ball in the potential',[]}';
args=i_parseargs(fieldnams,...
    'if(hasoption(''-p'')),deffields=''hax'';elseif(hasoption(''-u'')),deffields=''value'';else,deffields=''rangeV,npoints'';end;','-u,-p',varargin);
i_parcheck;
if ~isfield(args,'rangeV')
    args.rangeV = [ -1E30, 1E30];
end

if ~isfield(args,'npoints')
    args.npoints = 50;
end

iX = i_getno(g_grind.xaxis.var);
if iX.isvar
    var1 = g_grind.xaxis.var;
    lim1 = g_grind.xaxis.lim;
else
    var1 = g_grind.yaxis.var;
    iX = i_getno(var1);
    lim1 = g_grind.yaxis.lim;
end

if ~iX.isvar
    var1 = i_statevars_names(1);
    lim1 = [0 10];
end

if ~isfield(args,'range')
    args.range = lim1;
end

if any(strcmp(args.opts, '-p'))
    if ~isfield(args,'hax')
%            potential(args)
%            hfig=gcf;
         [hfig, new] = i_makefig('potent1');
         if new
             potential
         end

        args.hax = get(hfig, 'currentaxes');
    end

    addparanalreplay(args.hax)
    return;
end

if any(strcmp(args.opts, '-u'))
    if ~isfield(args,'value')
        args.value = i_initvar;
    end
    args.opts={};
    potential(args)
    hfig=gcf;
   % [hfig, new] = i_makefig('potent1');
   % if new
   %     potential
   % end
    
    hax = get(hfig, 'currentaxes');
    hdata=findobj(hax,'tag','potentialline');
    x = get(hdata, 'xdata');
    V = get(hdata, 'ydata');
    moveball(gca, args.value(1), x, V)
    return;
end

% if nargin < 3
%    only2D = 0;
% else
%    only2D = i_checkstr(only2D);
% end


if isfield(g_grind.solver, 'switch_stochast')
    oldpars = cell(size(g_grind.solver.switch_stochast));
    for i = 1:length(g_grind.solver.switch_stochast)
        oldpars{i} = evalin('base', g_grind.solver.switch_stochast{i});
        aval = zeros(size(oldpars{i}));
        assignin('base', g_grind.solver.switch_stochast{i}, aval);
    end

end


if (g_grind.statevars.dim == 1)&&~g_grind.solver.haslags&&~g_grind.solver.isimplicit
    % i_makepotfun(var1, potfun);
    hpotfun=i_getodehandle('potential','');
    [x, V] = ode45(hpotfun, linspace(lim1(1),lim1(2),args.npoints), 0);
    %  [x, V] = ode45(potfun, linspace(lim1(1),lim1(2),args.npoints), 0);
    if nargout == 0
        [hfig, new] = i_makefig('potent1');
        if new
            set(hfig,'name', 'Potential')
            set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
            set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
        end

        h1 = plot(x, V, g_grind.pen.pen, 'Color', g_grind.pen.color);
        set(h1,'tag','potentialline');
        i_plotdefaults(hfig);
        xlabel(i_disptext(var1));
        ylabel('potential');
        if args.rangeV(1) > -1E29
            set(gca, 'ylim', args.rangeV);
        end
        ud = get(gca, 'userdata');
        limy = get(gca, 'ylim');
        ud.meta=struct('func','potential','xname',var1,'xlim',lim1,'ylim',limy,'yname','potential','zname','');  
        set(gca, 'userdata', ud);
        set(gca, 'xlim', lim1);
        hold on
        h = fill([x(1); x; x(end)], [limy(1); V; limy(1)], [0.8 0.8 1]);
        set(h,'tag','potentialfill');
        xpos = i_initvar;
        drawball(gca, xpos, x, V);
        hold off
        addreplaydata(get(hfig, 'CurrentAxes'));
    else
        xs = x;
        pot = V;
    end
elseif g_grind.solver.isimplicit
   error('GRIND:potential:isimplicit','Error in potential, can not create potential for implicit differential equation (dae)');
elseif g_grind.solver.haslags
   error('GRIND:potential:haslags','Error in potential, can not create potential for delay differential equation');
elseif g_grind.statevars.dim > 1
   error('GRIND:potential:dims','Error in potential, can only calculate potential for 1 equation models (1D)');
else
    i_errordlg('Can only have state variables on the axes')
    error('GRIND:potential:NoStatevars','Error in potential, no state variables on the axes');
    % else
    %    i_makepotfun(var2, potfun);
    %    %initial conditions in the y direction
    %    assignin('base', var1, 0);
    %    [X, V] = ode45(potfun, lim2, 0);
    %    V0 = interp1(X, V, lim2(1):(lim2(2) - lim2(1)) / (npoints - 1):lim2(2), 'pchip');
    %    i_makepotfun(var1);
    %    Vs = zeros(npoints, npoints);
    %    Xs = zeros(npoints, npoints);
    %    Ys = zeros(npoints, npoints);
    %    i_waitbar(0, npoints, 'potential','Calculating',0.5)
    %    for i = 1:npoints
    %       waitbar(1);
    %       y = lim2(1) + (i - 1) * (lim2(2)-lim2(1)) / (npoints-1);
    %       assignin('base', var2, y);
    %       [X, V] = ode45(potfun, lim1, V0(i));
    %       Ys(i, :) = y;
    %       Xs(i, :) = lim1(1):(lim1(2) - lim1(1)) / (npoints - 1):lim1(2);
    %       Vs(i, :) = interp1(X, V, Xs(i, :), 'pchip');
    %       for j = 1:npoints
    %          if Vs(i, j) > VRange(2)
    %             Vs(i, j) = NaN;
    %          elseif Vs(i, j) < VRange(1)
    %             Vs(i, j) = NaN;
    %          end

    %       end

    %    end

    %    i_waitbar([]);
    %    hfig = figure(i_figno('potent3'));
    %    set(hfig, 'Name', 'Potentials surface');
    %    surf(Xs, Ys, Vs);
    %    xlabel(i_disptext(var1));
    %    ylabel(i_disptext(var2));
    %    zlabel('potential');
    %    shading flat;
    %    colorbar;
    %    if g_grind.statevars.dim > 2
    %       ti = ['Potentials valid for ' i_othervars(i_initvar, i_varno(var1), i_varno(var2))];
    %    else
    %       ti = 'Potentials';
    %    end

    %    htitle = title(ti);
    %    set(htitle,'fontweight','normal');
    %    figure (i_figno('potent2'));
    %    set(hfig, 'Name', 'Potentials contour');
    %    set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
    %    set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
    %    contourf(Xs, Ys, Vs, 50);
    %    htitle = title(ti);
    %    set(htitle,'fontweight','normal');
    %    xlabel(i_disptext(var1));
    %    ylabel(i_disptext(var2));
    %    colorbar;
end

% clear(potfun);
% delete([fullfile(grindpath, potfun) '.m']);
if isfield(g_grind.solver, 'switch_stochast')
    for i = 1:length(g_grind.solver.switch_stochast)
        assignin('base', g_grind.solver.switch_stochast{i}, oldpars{i});
    end

end


function drawball(hax, xpos, Xs, V)
oldhold = ishold(hax);
hold(hax, 'on')
pos = ballpos(hax, xpos, Xs, V);
h = zeros(2, 1);
pos1 = getpixelposition(gca);
h(1)=scatter(hax,pos(1,1),pos(1,2),200/325.5*pos1(4),'k','markerfacecolor','k','tag','bigball');
h(2)=scatter(hax,pos(2,1),pos(2,2),2.5/325.5*pos1(4),'w','markerfacecolor','w','tag','smallball');

if ~oldhold
    hold(hax, 'off');
end


function pos = ballpos(hax, xpos, Xs, V)
ylim = get(hax, 'ylim');
xlim = get(hax, 'xlim');
pos = zeros(2);
pos(1, 1) = xpos;
pos(1, 2) = interp1(Xs, V, xpos) + (ylim(2) - ylim(1)) * 0.031;
pos(2, 2) = interp1(Xs, V, xpos) + (ylim(2) - ylim(1)) * 0.04;
pos(2, 1) = xpos + (xlim(2) - xlim(1)) * 0.005;

function moveball(hax, xpos, Xs, V)
h=findobj(hax,'tag','bigball');
if isempty(h)||~any(ishandle(h))
    drawball(hax, xpos, Xs, V)
    h=findobj(hax,'tag','bigball');
end

h2=findobj(hax,'tag','smallball');
pos = ballpos(hax, xpos, Xs, V);
set(h, 'Xdata', pos(1, 1));
set(h, 'Ydata', pos(1, 2));
set(h2, 'Xdata', pos(2, 1));
set(h2, 'Ydata', pos(2, 2));

function addreplaydata(hax)
ud = get(hax, 'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
set(hax,'userdata', ud);

function replaystart(hax)
global g_t g_Y;
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
        hdata=findobj(hax,'tag','potentialline');
        if iscell(hdata)
            hdata = hdata{1};
        end

        ud.replay.pot.x = get(hdata, 'xdata');
        ud.replay.pot.V = get(hdata, 'ydata');
        %the balls must be refreshened else they are not working in a
        %combined figure
        h=findobj(hax,'tag','bigball');
        delete(h);
        h=findobj(hax,'tag','smallball');
        delete(h);
        ud.replay.ydata=g_Y(:,1);
        if any(~isreal(ud.replay.ydata))
             warning('grind:potential','Imaginary parts of complex values ignored');
             ud.replay.ydata=real(ud.replay.ydata);
        end
        moveball(hax, ud.replay.ydata(1), ud.replay.pot.x, ud.replay.pot.V)
        ud.replay.tdata = outfun('t');
    end

    set(hax, 'userdata', ud);
    if ~isoctave&&verLessThan('matlab','8.4.0')
        set(hax, 'drawmode','fast');
    else
        set(hax,'SortMethod','depth');
    end

end


function replayend(hax, closedlg)
if closedlg
    replaycallback(hax, 't', 1);
end


function t = replaycallback(hax,avar, relt)
t = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, 't')
    ud = get(hax, 'userdata');
    t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
    if isfield(ud, 'replay')
        if ~isfield(ud.replay, 'ydata')
            replaystart(hax);
            ud = get(hax, 'userdata');
        end

        ndx=find(ud.replay.tdata <= t);
        if ~isempty(ndx)
            xpos = ud.replay.ydata(ndx(end));
        end

        moveball(hax, xpos, ud.replay.pot.x, ud.replay.pot.V);
        %        if length(ndx)>3
        %               h=findobj(hax,'tag','bigball');
        %               prevballs=ud.replay.ydata(ndx(end-2:end));
        %               if prevballs(2)-prevballs(1)<prevballs(3)-prevballs(2)
        %                   set(h,'markerfacecolor','r');
        %               else
        %                   set(h,'markerfacecolor','k');
        %               end
        %        end

    end

end


function ud = updatepotentials(ud, ~, ~)
global g_paranal;
if  ~isfield(g_paranal, 'potdata')||isempty(g_paranal.potdata)
    g_paranal.potdata.data = cell(size(g_paranal.run.parvalues));
end

i_waitbar(0, size(g_paranal.run.parvalues,1), 'null', 'Preparing replay of potentials',0.5);
oldpar = evalin('base', g_paranal.run.pars{1});
try
    changed = false;
    for i = 1:size(g_paranal.run.parvalues, 1)
        i_waitbar(1);
        if isempty(g_paranal.potdata.data{i})
            changed = true;
            multassignin('base', g_paranal.run.pars{1}, g_paranal.run.parvalues(i,1));
            [X, V] = potential;
            g_paranal.potdata.data{i} = [X, V];
        end

    end

    if changed
        g_paranal.potdata.parvalues = g_paranal.run.parvalues;
    end
    multassignin('base', g_paranal.run.pars{1}, oldpar);
catch err
    multassignin('base', g_paranal.run.pars{1}, oldpar);
    rethrow(err);
end

i_waitbar([]);

function addparanalreplay(hax)
global g_paranal;
if nargin < 1
    hax = [];
end

if isempty(hax)||~ishandle(hax)
    hfig = i_figno('potent1');
    if ishandle(hfig)
        hax=findobj(hfig,'type','axes');
        tags = get(hax, 'tag');
        if ~isempty(tags)
            hax = hax(~strcmp(tags, 'legend'));
        elseif isempty(hax)
            hax = gca;
        end

        ud = get(hax, 'userdata');
        hh=findobj(hfig,'Tag','nullp');
        delete(hh);
        i_grindlegend(12, hax);
    else
        hfig = i_makefig('potent1');
        hax=findobj(hfig,'type','axes');
        if isempty(hax)
            hax = gca;
        end

        tags = get(hax, 'tag');
        if ~isempty(tags)
            hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
        end
    end

else
    hh=findobj(hax,'Tag','potentialline');
    delete(hh);
end

ud.replay.hball=[findobj(hax,'tag','bigball'),findobj(hax,'tag','smallball')];
ud.replay.callback = @replayparanalcallback;
ud.replay.onstart = @replayparanalstart;
ud.replay.onend = [];
ud.replay.onturn = @i_replayparanalturn;
tstep = g_paranal.run.parvalues(2, 1)  - g_paranal.run.parvalues(1, 1);
ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
    g_paranal.nulldata.data = cell(size(g_paranal.run.parvalues, 1));
end

ud = updatepotentials(ud, false, hax);
set(hax,'userdata', ud);
replayall('variable', g_paranal.run.pars{1});

function replayparanalstart(hax)
global g_paranal ;
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        tstep = g_paranal.run.parvalues(2, 1)  - g_paranal.run.parvalues(1, 1);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
        if isempty(g_paranal.potdata)||isempty(g_paranal.potdata.data)
            %pars = unique(g_paranal.run.p);
            %if ~strcmp(g_paranal.run.parname, ud.replay.settings.tvar)||(length(pars)~=length(ud.replay.pars))||(pars(1)~=ud.replay.pars(1))||pars(end)~=ud.replay.pars(end)
            updatepotentials(hax);
            return;
        end

        maxpot = max(cellfun(@(x)max(x(:, 2)), g_paranal.potdata.data));
        minpot = min(cellfun(@(x)min(x(:, 2)), g_paranal.potdata.data));
        ylimits = get(hax, 'ylim');
        maxpot = max(maxpot, ylimits(2));
        minpot = min(minpot, ylimits(1));
        set(hax,'ylimmode','manual');
        set(hax, 'ylim', [ minpot,maxpot])
        %       if any(~ishandle(ud.replay.hball))
        %           drawball(hax,g_paranal.run.Y(1,1,1),g_paranal.potdata.data{1}(:,1),g_paranal.potdata.data{1}(:,2));
        %           ud.replay.hball=[findobj(hax,'tag','bigball'),findobj(hax,'tag','smallball')];
        %       else
        h=findobj(hax,'tag','bigball');
        delete(h);
        h=findobj(hax,'tag','smallball');
        delete(h);
        moveball(hax, g_paranal.run.Y(1, 1, 1), g_paranal.potdata.data{1}(:, 1), g_paranal.potdata.data{1}(:, 2))
        %      end
        %       if ishandle(ud.replay.HLine(1))
        %           ax1 = get(ud.replay.HLine(1), 'parent');
        %           if ax1~=hax
        %               ud.replay.hwhere=-1;
        %               delete(ud.replay.HLine);
        %               ud.replay.HLine(1)=-1;
        %               hh=findobj(hax,'Tag','nullp');
        %               delete(hh);
        %               i_grindlegend(12,hax);
        %               ud=updatenullclines(ud, false, hax);
        %           end

        %       end

        %       if ~ishandle(ud.replay.hnull(1))||(get(ud.replay.hnull(1), 'parent')~=hax)
        %           hh=findobj(hax,'Tag','nullp');
        %           delete(hh);
        %           i_grindlegend(12,hax);
        %           ud=updatenullclines(ud,false,hax);
        %       end

        
    end

    set(hax, 'userdata', ud);
end



function p = replayparanalcallback(hax, avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, g_paranal.run.pars{1})
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
        [tndx, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
        p = g_paranal.run.parvalues(stepndx, 1);
        ballpos = g_paranal.run.Y(tndx, 1, stepndx);
        [~, pnow] = min(abs(g_paranal.potdata.parvalues(:,1) - p));
        Vdata = g_paranal.potdata.data{pnow};
        h=findobj(hax,'tag','potentialline');
        set(h, 'XData',Vdata(:,1) ,'YData', Vdata(:,2));
        ylim = get(hax, 'ylim');
        hfill=findobj(hax,'tag','potentialfill');
        set(hfill, 'XData',[Vdata(1,1); Vdata(:,1); Vdata(end,1)],'YData',[ylim(1); Vdata(:,2);ylim(1)]);
        
        moveball(hax, ballpos, Vdata(:, 1), Vdata(:, 2));
        %        if tndx>2
        %               h=findobj(hax,'tag','bigball');
        %               prevballs=g_paranal.run.Y(tndx-2:tndx,1,stepndx);
        %               if prevballs(2)-prevballs(1)<prevballs(3)-prevballs(2)
        %                   set(h,'markerfacecolor','r');
        %               else
        %                   set(h,'markerfacecolor','k');
        %               end
        %        end

        
    end

end


