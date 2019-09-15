%NULL   2D phase plane with nullclines
%   Creates nullclines in the 2D phase plane (it is also allowed to have a
%   parameter on the x or y-axes).
%   If there is no variable on the y-axis, a one-dimensional variant is used
%   (<a href="matlab:help itermap">itermap</a> or <a href="matlab:help plotdiff">plotdiff</a>).
%   If there are more than two state variables you can either make a cross-section
%   through the multi-dimensional nullcline or find a stable nullcline for the
%   other variables ("quasi-steady state mode").
%   
%   Usage: 
%   NULL - creates nullclines with default accuracy
%   NULL NPOINTS - creates nullclines using a grid of NPOINTS x NPOINTS points (default=50)
%   NULL [N1 N2]  - creates nullclines using a grid of N1 x N2 points 
%   NULL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'hax' [handle] - Handle of the plot (used in replayall)
%     'npoints' [all(integer>=10) and length(integer)<=2] - Size of the grid used to calculate the nullclines.
%     'vect' [number] - Optionally you can add the matrices for calculating
%   NULL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - cross-section mode for high dimensional models: keep all other variables at their current values.
%     '-i' - inflection points (beta)
%     '-p' - makes it possible to replay a paranal session. (see also <a href="matlab:help replayall">replayaall</a>)
%     '-q' - force quasi-steady state mode for high dimensional models (uses <a href="matlab:help findeq">findeq</a>).
%      
%
%  
%   Remark:
%   NULL is also a MATLAB function to determine the null space of a matrix, for help on this function see:
%   <a href="matlab:help matfun/null">help matfun/null</a>
%   You can still use this function (if the first argument is a matrix).
%  
%  
%   See also itermap, plotdiff, null3, ax, replayall, growths, findeq, <a href="matlab:help matfun/null">matfun/null</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('null')">commands null</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function Z = null(varargin)
%(npoints, opt, opt2)
global g_grind ;

if nargout > 0 && nargin>0 && ~ischar(varargin{1})
    %there is a problem that NULL was already defined in MATLAB,
    %if the first argument is a matrix or if there is output
    %it is assumed that the matlab NULL is meant.
    cur = pwd;
    cd(fullfile(matlabroot,'toolbox','matlab','matfun'));
    Z = null(varargin{:});
    cd(cur);
    return;
end

iX = i_getno(g_grind.xaxis.var,true);
iY = i_getno(g_grind.yaxis.var,true);
if g_grind.statevars.dim==1&& (~(iX.isvar || iX.ispar || iX.isext)  ||  ~(iY.isvar || iY.ispar || iY.isext))
    if isempty(g_grind)  ||  g_grind.solver.isdiffer
        itermap(varargin{:});
    else
        plotdiff(varargin{:});
    end
    return;
end
if strcmp(g_grind.solver.opt.Vectorized,'on')
    npoint=250;
else
    npoint=25;
end
fieldnams={'npoints', 'all(i>=10)&length(i)<=2', 'Size of the grid used to calculate the nullclines.',npoint;...
   'hax', 'h', 'Handle of the plot (used in replayall)',[];...
   'vect', 'n', 'Optionally you can add the matrices for calculating',[]}';
[args,defaults]=i_parseargs(fieldnams,'npoints', '-p,-q,-c,-i',varargin);
i_parcheck;
hasdata = 0; %if data are entered i_vector is not evaluated
if (isempty(iX.no)  ||  isempty(iY.no))
    ax('-?');
    i_errordlg('Cannot create null-isoclines if there are no state variables on the axes, use "phas 2" instead.');
    error('GRIND:null:NoStatevars','null: Cannot create null-isoclines if there are no state variables on the axes, use <a href="matlab:phas 2">phas 2</a> instead or change the axes with  <a href="matlab:ax">ax</a>.');
end
quasiss = [];
inflectionpoints  = false;
args=mergestructs(defaults,args);
% if ~isfield(args,'npoints')||isempty(args.npoints)
%     args.npoints = 250;
% end
% 
% if ~isfield(args,'hax')
%     args.hax=[];
% elseif ischar(args.hax)
%     args.hax=i_checkstr(hax);
% end

if any(strcmp(args.opts, '-q'))
    quasiss = true;
elseif  any(strcmp(args.opts, '-c'))
    quasiss = false;
end

if  any(strcmp(args.opts, '-p'))
    updateparanalreplay(args.npoints, args.hax, quasiss);
    return;
elseif  any(strcmp(args.opts, '-i'))
    inflectionpoints = true;
end

if  ~isempty(args.vect)
    %plug in data from a user-defined function (used for R* evaluations)
    Vect = args.vect;
    hasdata = 1;
    args.npoints = size(Vect);
end

if ~hasdata&&g_grind.statevars.dim > 2&&isempty(quasiss)&&~inflectionpoints
    quasiss= strcmp(questdlg('There are too many state variables to plot, do you want to assume quasi-steady states for the other?','null','Yes','No','Yes'),'Yes');
else
    quasiss = false;
end

if hasdata
    [X, Y] = meshgrid(linspace(g_grind.xaxis.lim(1), g_grind.xaxis.lim(2), args.npoints(1)), linspace(g_grind.yaxis.lim(1), g_grind.yaxis.lim(2), args.npoints(2)));
else
    if quasiss
        i_waitbar(0, prod(args.npoints)^2, 'null', 'Finding quasi-steady states',0.5,true)
        [X, Y, Vect, qssvalues] = i_vector(args.npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, [],4);
        if isempty(Vect)
            disp('Cancel pressed');
            return;
        end
        i_waitbar([]);
    else
        if inflectionpoints
            [X, Y, Vect] = i_vector(args.npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, [],2);
        else
            [X, Y, Vect] = i_vector(args.npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, [],0);
        end

    end

end

if isempty(args.hax)||~ishandle(args.hax)
    [hfig, new] = i_makefig('phase2');
    args.hax = gca;
    if ~new
        set(args.hax, 'View', [0 90]);
    end

else
    hfig = get(args.hax, 'parent');
    new = 0;
end
if new
    set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
    set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
    nser = 0;
    hold(args.hax, 'on');
    oldhold = 1;
else
    oldhold = ishold;
    if ~oldhold
        %milder replace children than default
        delete(get(args.hax, 'children'))
        set(args.hax, 'userdata', [])
        legend off
        hold on;
    end
    if ~i_grindlegend(4, args.hax)&&~inflectionpoints %not settings changed?
        ser = i_grindlegend(3, args.hax);
        uds = get(ser, 'userdata');
        if iscell(uds)
            uds1 = zeros(size(uds));
            for i = 1:length(uds)
                uds1(i) = uds{i}.nr;
            end

            uds = uds1;
        else
            uds = uds.nr;
        end

        nser=sum(uds ~= max(uds));
        delete(ser(uds == max(uds)));
    else
        nser = length(i_grindlegend(3, args.hax));
    end

end
ud = get(args.hax, 'userdata');
ud.meta=struct('func','null','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname',g_grind.yaxis.var,'ylim',g_grind.yaxis.lim,'zname','');
set(args.hax, 'userdata', ud);
set(hfig, 'Name', 'Phase plane');
% args.hax=findobj(H,'type','axes');
% tags = get(args.hax, 'tag');
% args.hax=args.hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
% hfig = H;



%hold(args.hax,'on');
%Hs = zeros(g_grind.statevars.dim, 1);
colors = [0, 0.5, 0; 0, 0, 1; 1, 0, 0; 0, 0.75, 0.75; 0.75, 0, 0.75; 0.75, 0.75, 0; 0.25, 0.25, 0.25];
%colors = {'m', 'b', 'r', 'g','k'};
if isoctave||~verLessThan('matlab','8.4.0')
    lines = {'--','-.',':','-'};
else
    lines = {':','-.','--','-'};
end
maxi = size(colors, 1) * length(lines);
for i = 1:size(Vect, 3)
    if ~quasiss||i==iX.no||i==iY.no
        j = mod(i + nser, maxi);
        % pen = [char(lines{floor(j / size(colors, 2) + 1)}) char(colors{mod(j, size(colors, 2)) + 1})];
        pen  = char(lines{floor(j / size(colors, 1) + 1)});
        col = colors(mod(j, size(colors, 1)) + 1,:);
        set(hfig,'doublebuffer','on');
        if ~isoctave&&verLessThan('matlab','8.4.0')
            set(args.hax, 'drawmode','fast');
        else
            set(args.hax,'SortMethod','depth');
        end
        if g_grind.solver.iscomplex
           h= plotcomplexnullclines(args.hax,X,Y,Vect,i,pen,col);
        else
        [c] = getnullclines(X, Y,  Vect(:, :, i));
        if isempty(c)
            c = [g_grind.xaxis.lim(1), g_grind.yaxis.lim(1); g_grind.xaxis.lim(1), g_grind.yaxis.lim(1)];
        end
        if g_grind.solver.isdiffer
            leg = i_disptext([i_statevars_names(i) '_{t} = ' i_statevars_names(i) '_{t+1}']);
        else
            leg = i_disptext([i_statevars_names(i) '''' '=0']);
        end

        h = plot(args.hax(1),c(1,:), c(2,:), pen, 'color', col);
        set(h, 'linewidth', g_grind.pen.linewidth);
        i_grindlegend(2, h, leg);
        end
        if nargout > 0
          Z(i) = h;
        end
    end

end

if new
    i_plotdefaults(hfig);
    if ~isoctave&&verLessThan('matlab','8.4.0')
        set(args.hax, 'drawmode','fast');
    else
        set(args.hax,'SortMethod','depth');
    end

end

if new||isempty(get(get(gca,'xlabel'),'string'))
    h = xlabel(i_disptext(g_grind.xaxis.var));
    set(h, 'fontsize', g_grind.pen.fontsize);
    h = ylabel(i_disptext(g_grind.yaxis.var));
    set(h, 'fontsize', g_grind.pen.fontsize);
end

i_grindlegend(11);

if ~hasdata  &&  (g_grind.statevars.dim > 2) &&iX.isvar&&iY.isvar
    if quasiss
        ud = get(gca, 'userdata');
        ud.quasiss = true;
        ud.qssvalues = permute(qssvalues, [2 1 3]); %Vect is also transposed
        ud.mask = true(g_grind.statevars.dim, 1);
        ud.mask(iX.no) = false;
        ud.mask(iY.no) = false;
        set(gca, 'userdata', ud);
        htitle = title(args.hax, 'Other variables in steady state');
    else
        htitle = title(args.hax,i_disptext(['Valid for ' i_othervars(i_initvar, iX.no, iY.no)]));
    end

    set(htitle,'fontweight','normal');
end

set(args.hax, 'XLim', g_grind.xaxis.lim);
set(args.hax, 'YLim', g_grind.yaxis.lim);
if ~oldhold
    hold(args.hax, 'off');
end
function h=plotcomplexnullclines(hax,X,Y,Vect,ii,pen,col)
global g_grind
[c] = getnullclines(X, Y, real( Vect(:, :, ii)));
if isempty(c)
    c = [g_grind.xaxis.lim(1), g_grind.yaxis.lim(1); g_grind.xaxis.lim(1), g_grind.yaxis.lim(1)];
end
if g_grind.solver.isdiffer
    leg = i_disptext(['real ' i_statevars_names(ii) '_{t} = ' i_statevars_names(ii) '_{t+1}']);
else
    leg = i_disptext(['real ' i_statevars_names(ii) '''' '=0']);
end

h(1) = plot(hax,c(1,:), c(2,:), pen, 'color', col);
set(h(1), 'linewidth', g_grind.pen.linewidth);
i_grindlegend(2, h(1), leg);

[c] = getnullclines(X, Y, imag( Vect(:, :, ii)));
if isempty(c)
    c = [g_grind.xaxis.lim(1), g_grind.yaxis.lim(1); g_grind.xaxis.lim(1), g_grind.yaxis.lim(1)];
end
if g_grind.solver.isdiffer
    leg = i_disptext(['imag ' i_statevars_names(ii) '_{t} = ' i_statevars_names(ii) '_{t+1}']);
else
    leg = i_disptext(['imag ' i_statevars_names(ii) '''' '=0']);
end

h(2) = plot(hax,c(1,:), c(2,:), pen, 'color', 'r');
set(h(2), 'linewidth', g_grind.pen.linewidth);
i_grindlegend(2, h(2), leg);

function c = getnullclines(X, Y, Vecti)
if size(Vecti, 3) > 1
    c = cell(size(Vecti, 3), 1);
    for i = 1:size(Vecti, 3)
        c{i} = getnullclines(X, Y, Vecti(:, :, i));
    end

else
    if ~all(isinf(Vecti(:))|isnan(Vecti(:)))
      if ~any(isreal(Vecti(:)))
          Vecti=real(Vecti);
          warning('grind:null:complex','Imaginary parts of complex values ignored');
      end

      [c] = contourc(X(:,1), transpose(Y(1,:)),  transpose(Vecti), [0, 0]);
    else
        c = [];
    end

    %I do not use contour as it gives lines that are not clippable in Adobe
    %Illustrator.
    if ~isempty(c)
        limit = size(c, 2);
        k = 1;
        while (k < limit)
            npoints = c(2, k);
            c(:, k) = nan; %nans in between parts of the nullclines (if they consist of some parts)
            k = k + npoints + 1;
        end
    end

end

function ud = updatenullclines(ud, plotted, hax, quasiss)
global g_paranal g_grind;
npoints = ud.replay.npoints;
i_waitbar(0, size(g_paranal.run.parvalues,1), 'null', 'Preparing replay of paranal',0.5);
oldpar = evalin('base', g_paranal.run.pars{1});
try
    changed = false;
    for i = 1:size(g_paranal.run.parvalues, 1)
        i_waitbar(1);
        if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data{i})
            multassignin('base', g_paranal.run.pars{1}, g_paranal.run.parvalues(i,1));
            [X, Y, Vect] = i_vector(npoints, ud.replay.iX, g_grind.xaxis.lim, ud.replay.iY, g_grind.yaxis.lim, [],0);
            changed = true;
            g_paranal.nulldata.data{i} = getnullclines(X, Y, Vect);
        end

        if ~plotted
            if isempty(quasiss)
                ud.replay.hnull = null(num2str(npoints), 'hax', hax);
            elseif quasiss
                ud.replay.hnull = null(num2str(npoints),'-q', 'hax', hax);
            else
                ud.replay.hnull = null(num2str(npoints),'-c', 'hax', hax);
            end

            plotted = all(ishandle(ud.replay.hnull));
            if plotted
                set(ud.replay.hnull,'tag','nullp');
            else
                delete(ud.replay.hnull(ishandle(ud.replay.hnull)));
                i_grindlegend(12, hax);
            end

        end

    end

    if changed
        g_paranal.nulldata.parvalues = g_paranal.run.parvalues;
    end
    multassignin('base', g_paranal.run.pars{1}, oldpar);
catch err
    multassignin('base', g_paranal.run.pars{1}, oldpar);
    rethrow(err);
end

set(hax,'nextplot','add');
if ud.replay.iX.isvar&&ud.replay.iY.isvar
    ud.replay.HLine = plot(hax,g_paranal.run.Y(:, ud.replay.iX.no,1), g_paranal.run.Y(:, ud.replay.iY.no,1), 'b-');
    set(ud.replay.HLine,'tag','nullp');
    i_grindlegend(-1, ud.replay.HLine);
end

i_waitbar([]);

function updateparanalreplay(npoints, hax, quasiss)
global g_paranal g_grind;
if isempty(g_paranal.run)||isempty(g_paranal.run.Y)
    error('grind:null:paranal','Paranal should be run before using this option');
end

if nargin == 0
    npoints = 50;
end

if nargin<3
    quasiss=[];
end

if nargin < 2
    hax = [];
end

if isempty(hax)||~ishandle(hax)
    H = i_figno('phase2');
    if ishandle(H)
        hax=findobj(H,'type','axes');
        tags = get(hax, 'tag');
        hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
        ud = get(hax, 'userdata');
        hh=findobj(hax,'Tag','nullp');
        delete(hh);
        i_grindlegend(12, hax);
    else
        H = i_makefig('phase2');
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
ud.replay.iX = i_getno(g_grind.xaxis.var,true);
ud.replay.iY = i_getno(g_grind.yaxis.var,true);
ud.replay.hwhere = -1;
ud.replay.HLine = -1;
tstep = g_paranal.run.parvalues(2, 1) - g_paranal.run.parvalues(1, 1);
ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
ud.replay.npoints = npoints;
if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
    g_paranal.nulldata.data = cell(size(g_paranal.run.parvalues));
end

%i_waitbar(0, size(g_paranal.run.parvalues,1), 'null', 'Preparing replay of paranal',0.5);
plotted = false;
ud = updatenullclines(ud, plotted, hax, quasiss);
% hax=findobj(H,'type','axes');
% tags = get(hax, 'tag');
% hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
set(hax,'userdata', ud);
%i_waitbar([]);
replayall('variable', g_paranal.run.pars{1});

function replaystart(hax)
global g_paranal g_grind;
if ishandle(hax)
    hfig = get(hax, 'parent');
    i_figure(hfig);
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        tstep = g_paranal.run.parvalues(2, 1) - g_paranal.run.parvalues(1, 1);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',size(g_paranal.run.parvalues,1)*10,'ndata',size(g_paranal.run.parvalues,1));
        %pars=unique(g_paranal.run.p);
        if isempty(g_paranal.nulldata)||isempty(g_paranal.nulldata.data)
            if ishandle(ud.replay.HLine)
                delete(ud.replay.HLine);
            end

            null('-p', 50, hax);
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
                ud = updatenullclines(ud, false, hax, []);
            end

        end

        if ~ishandle(ud.replay.hnull(1))||(get(ud.replay.hnull(1), 'parent')~=hax)
            hh=findobj(hax,'Tag','nullp');
            delete(hh);
            i_grindlegend(12, hax);
            ud = updatenullclines(ud, false, hax, false);
        end

        set(hax,'ylimmode','manual');
        set(hax,'zlimmode','manual');
        oldnext = get(hax, 'NextPlot');
        set(hax,'NextPlot','add');
        ud.replay.iX = i_getno(g_grind.xaxis.var,true);
        ud.replay.iY = i_getno(g_grind.yaxis.var,true);
        if ~ishandle(ud.replay.HLine)
            ud.replay.HLine = plot(hax,g_paranal.run.Y(1,ud.replay.iX.no), g_paranal.run.Y(1,ud.replay.iY.no), '-');
            set(ud.replay.HLine,'tag','nullp');
        end

        if ~ishandle(ud.replay.hwhere)
            ud.replay.hwhere = plot(hax,g_paranal.run.Y(1,ud.replay.iX.no), g_paranal.run.Y(1,ud.replay.iY.no), 'ro');
            set(ud.replay.hwhere, 'MarkerFaceColor', [1 0 0]);
            set(ud.replay.hwhere,'tag','nullp');
            i_grindlegend(-1, ud.replay.hwhere);
        end

        set(hax, 'NextPlot', oldnext);
    end

    set(hax, 'userdata', ud);
end

function replayend(hax, closedlg)
if closedlg
    replaycallback(hax,'', 1);
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


function p = replaycallback(hax,avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, g_paranal.run.pars{1})
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
        %  t = g_paranal.run.t(ndx1);
        [tndx, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
        p = g_paranal.run.parvalues(stepndx);
        [~, pnow] = min(abs(g_paranal.nulldata.parvalues(:,1) - p));
        cc = g_paranal.nulldata.data{pnow};
        for i = 1:length(cc)
            cdata = cc{i};
            if ~isempty(cdata)
                set(ud.replay.hnull(i), 'XData',cdata(1,:) ,'YData', cdata(2,:));
            end

        end

        % ndx=(p == g_paranal.run.parvlues) & (g_paranal.run.t<=t);
        xdata = g_paranal.run.Y(1:tndx, ud.replay.iX.no,stepndx);
        ydata = g_paranal.run.Y(1:tndx, ud.replay.iY.no,stepndx);
        if any(~isreal(xdata(:)))||any(~isreal(ydata))
            warning('grind:null:complex','Imaginary parts of complex values ignored');
            xdata=real(xdata);
            ydata=real(ydata);
        end
        if ~isempty(xdata)
            if ishandle(ud.replay.hwhere)
                set(ud.replay.hwhere,'Xdata',xdata(end),'Ydata',ydata(end));
            end

            if ishandle(ud.replay.HLine)
                set(ud.replay.HLine, 'XData',xdata ,'YData',ydata);
            end

        end

    end

end

