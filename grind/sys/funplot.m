%FUNPLOT   Plot any 2D equation (including parameters)
%   Plot some equation (including parameters, state variables, time and
%   auxiliary variables). For state variables the current initial value is used.
%
%   Usage:
%   FUNPLOT FUN VAR XRANGE YRANGE LOG - FUN is an equation which is 
%   plotted on the y-axis ; VAR is the dependent variable for 
%   the X axis; XRANGE is the range of which the variable VAR1 is 
%   varied. YRANGE (optional) sets the range of the Y axis. LOG 
%   (optional) makes the scale x-axis or the y-axis logarithmic: 
%   values: 'logx' horizontal axis log scale, 'logy' vertical axis 
%   log scale, 'logxy' both log scales.
%   FUNPLOT FUN VAR MIRROR XRANGE YRANGE LOG - if MIRROR is 
%   Y then the function is plotted on the x axis and the dependent
%   variable on the y-axis.
%   FUNPLOT VARY FUN YRANGE XRANGE - If the variable and function
%   are exchanged, the function is plotted on the X axis and VARY is
%   the dependent variable for the Y axis. The first range is the 
%   range for the -axis now.
%   FUNPLOT FUN VARX - if AXRange is omitted, a range of [0.1
%   100] is assumed.
%   FUNPLOT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - equation to be plotted
%     'hax' [handle] - the handle of the axes (used in replayall)
%     'log' [logx | logy | logxy or empty] - makes the scale x-axis or the y-axis logarithmic:
%    values: 'logx' horizontal axis log scale, 'logy' vertical axis log scale, 'logxy' both log scales
%     'mirror' [logical] - x and y axes are exchanged if mirror is true
%     'npoints' [integer>0] - number of points in the plot
%     'var' [string] - the name of the dependent variable
%     'xrange' [number and length(number)==2] - the range of the x-axis
%     'yrange' [number and length(number)==2 or empty] - the range of the y-axis
%   FUNPLOT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-u' - update all lines using the current parameters (this is done
%   automatically if a <a href="matlab:help paranal">paranal</a> plot is replayed (<a href="matlab:help replayall">replayall</a>).
%
%   Examples:
%   funplot 'A*r*(1-A/K)' A [0 10] ('normal' plot)
%   funplot 'A*r*(1-A/K)' A Y [0 10] ('mirrored' plot)
%
%   See also funcs, implicitplot, funplot3, replayall
%
%   Reference page in Help browser:
%      <a href="matlab:commands('funplot')">commands funplot</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function l_HH = funplot(varargin)
%function funplot(afun,g_l_avar,[exchange],[axrange, ayrange],logxy, npoints);
%
global g_grind;
if ~exist('i_use', 'file')
    addpath([grindpath filesep 'sys2']);
end
ndx=strncmp(varargin,'logx',4)|strncmp(varargin,'logy',4);
if any(ndx)
    logs=varargin(ndx);
    varargin=[varargin(~ndx) {'log'} logs];
end
fieldnams={'fun', 's', 'equation to be plotted','';...
   'var', 's', 'the name of the dependent variable','';...
   'mirror', 'l', 'x and y axes are exchanged if mirror is true', false;...
   'xrange', 'n&length(n)==2', 'the range of the x-axis',[0.01 100];...
   'yrange', 'n&length(n)==2#E', 'the range of the y-axis',[];...
   'log', 'e[logx|logy|logxy]#E', 'makes the scale x-axis or the y-axis logarithmic:','';...
   'npoints', 'i>0', 'number of points in the plot',1000;...
   'hax', 'h', 'the handle of the axes (used in replayall)',[]}';
g_args=i_parseargs(fieldnams,'if(hasoption(''-u'')),deffields=''hax'';elseif(argtype(3,''l'')),deffields=''fun,var,mirror,xrange,yrange'';else,deffields=''fun,var,xrange,yrange'';end', '-u',varargin);
%use weird names to avoid that a parameter is changed
if ~isfield(g_args,'npoints')
    g_args.npoints = 1000;
end
g_defx = 'x';
if isfield(g_args,'fun')&&isfield(g_args,'var')
    %exchange axes if fun and var are exchanged
    if length(regexp(g_args.var,'[a-zA-Z_][a-zA-Z_0-9]*', 'match'))<length(g_args.var)&&length(regexp(g_args.fun,'[a-zA-Z_][a-zA-Z_0-9]*', 'match'))==length(g_args.fun)
        h=g_args.var;
        g_args.var=g_args.fun;
        g_args.fun=h;
        g_args.mirror = true;
    end
end
if isfield(g_args,'fun')
    g_args.fun = checkeq(g_args.fun);
    vars = symvar(g_args.fun);
    if (length(vars) == 1)&&~isfield(g_args,'var')
        g_args.var=vars{1};
    elseif length(vars) > 1
        if ~any(strcmp(vars, g_defx))
            g_defx = vars{1};
        end
    end
    clear vars;
else
    g_args.fun='';
end
if ~isfield(g_args,'var')
    g_args.var='';
end
if ~isfield(g_args,'mirror')
    g_args.mirror = false;
end
if ~isfield(g_args,'xrange')
    g_args.xrange = [0.01 100];
end
if ~isfield(g_args,'yrange')
    g_args.yrange = [];
end
if ~isfield(g_args,'log')
    g_args.log = '';
end
if any(strcmp(g_args.opts, '-u'))
    if ~isfield(g_args,'hax')
        hfig = i_figno('funplot');
        if ~ishandle(hfig)
            error('funplot:update','There is no plot to update, run funplot first');
        end
        figure(hfig);
        if ~isempty(g_grind)&&~isempty(g_grind.figopts)
          set(hfig,g_grind.figopts{:});
       end
        g_args.hax = get(hfig, 'children');
        tags = get(g_args.hax, 'tag');
        g_args.hax=g_args.hax(~(strcmp(tags,'legend')|strcmp(tags,'GRIND menu')));
    end
    ser = get(g_args.hax(1), 'children');
    for i = 1:length(ser)
        ud = get(ser(i), 'userdata');
        if isfield(ud,'fun')
            [g_l_ran, g_l_res] = i_funcalc(ud.fun, ud.var, ud.axrange, g_args.npoints);
            if any(~isreal(g_l_ran))||any(~isreal(g_l_res))
                g_l_ran=real(g_l_ran);
                g_l_res=real(g_l_res);
                warning('grind:funplot:complex','Imaginary parts of complex values ignored') 
            end
            if ud.doexchange
                set(ser(i),'XData',g_l_res,'YData',g_l_ran);
            else
                set(ser(i),'XData',g_l_ran,'YData',g_l_res);
            end
        end
    end
    return;
end

if isempty(g_args.var)||isempty(g_args.fun)
%     prompt = {'Function', ...
%         'Independent variable','Exchange x and y axis? (y/n)','Range for independent variable','Range for function',...
%         'Logarithmic scales (logx/logy/logxy)?'};
    if isfield(g_grind, 'funplot')
        answer = g_grind.funplot;
        if ~isempty(g_args.fun)
            answer{1} = g_args.fun;
        end
    else
        answer={g_args.fun,g_defx,iif(g_args.mirror,'Y','N'),mat2str(g_args.xrange),mat2str(g_args.yrange),g_args.log,mat2str(g_args.npoints)};
    end
    answer=i_funplotdlg(answer);
  %  answer = inputdlg(prompt, 'Function plot', 1, answer);
    if isempty(answer)
        error('GRIND:funplot:cancelled','Cancelled by user');
    elseif isempty(answer{1})
        error('GRIND:funplot:NoEquation','No equation entered');
    else
        if isempty(answer{4})
            answer{4}='[]';
        end
        if isempty(answer{5})
            answer{5}='[]';
        end
        g_args1=i_parseargs(fieldnams,'fun,var,mirror,xrange,yrange,log,npoints','',answer);
        g_args=mergestructs(g_args,g_args1);
        %        g_nargin = 6;
    end
    g_grind.funplot = answer;
    clear answer prompt g_defx;
end
g_args.fun = checkeq(g_args.fun);
g_args.log = i_checklog(g_args.log);

%doexchange = 0;
if isempty(g_args.mirror)
    isfunction=false; %should check if there is a function
    [~, ~, g_args.mirror] = checkfun(g_args.var, g_args.fun, isfunction);
end
[g_l_ran, g_l_res, g_l_pars] = i_funcalc(g_args.fun, g_args.var, g_args.xrange, g_args.npoints);
if g_args.mirror
    h1 = i_funplot(g_l_res,g_l_ran, g_args.fun,g_args.var,g_args.yrange,g_args.xrange,g_args.log,g_l_pars);
else
    h1 = i_funplot(g_l_ran,g_l_res, g_args.var,g_args.fun,g_args.xrange,g_args.yrange,g_args.log,g_l_pars);
end
ud = get(h1, 'userdata');
ud.fun = g_args.fun;
ud.var = g_args.var;
ud.axrange = g_args.xrange;
ud.doexchange = g_args.mirror;
set(h1, 'userdata', ud);
if nargout == 1
    l_HH = h1;
end
%++++++++++FUNCTIONS++++++++++
function [g_l_avar, g_l_afun, isfunction] = checkfun(g_l_avar, g_l_afun, isfunction)
global g_grind;

%check if there is an operand in g_l_avar
for i = 1:length(g_l_avar)
    if ~isempty(strfind('=*-/+^([.\|&', g_l_avar(i))) %#ok<STREMP> % No convert
        isfunction = 1;
        [g_l_avar, g_l_afun] = exchange(g_l_avar, g_l_afun);
        break;
    end
end
%check if g_l_avar is numeric
isnum = 1;
for i = 1:length(g_l_avar)
    if ~strcontains('1234567890-.', g_l_avar(i)) 
        isnum = 0;
    end
end
if isnum
    isfunction  = 1;
    [g_l_avar, g_l_afun] = exchange(g_l_avar, g_l_afun);
end
%check if g_l_avar is a defined function
if isfield(g_grind, 'funcs')
    afuncs = str2cell(g_grind.funcs);
    for i = 1:length(afuncs)
        f = strfind(char(afuncs{i}), '=');
        if ~isempty(f) && strcmp(strtrim(afuncs{i}(1:f(1) - 1)), g_l_avar)
            isfunction = 1;
            [g_l_avar, g_l_afun] = exchange(g_l_avar, g_l_afun);
            break;
        end
    end
end
g_l_afun = checkeq(g_l_afun);
return;

function [g_l_afun] = checkeq(g_l_afun)
f=strfind(g_l_afun, '=');
if ~isempty(f)
    g_l_afun = g_l_afun(f(1) + 1:length(g_l_afun));
end
return;

function h = i_funplot(x, y, xaxis, yaxis, axrange, ayrange, logxy, pars)
global g_grind;
if isempty(logxy)
    logxy = [0 0];
end
[hfig, new] = i_makefig('funplot');
if new
    %  hax=findobj(hfig,'type','axes');
    hax = gca;
    tags = get(hax, 'tag');
    ndx=~(strcmp(tags,'Colorbar')|strcmp(tags,'legend'));
    if ~isempty(ndx)
        hax = hax(ndx);
    end
    ud = get(hax, 'userdata');
    ud.replay.callback = @replaycallback;
    ud.replay.onstart = @replaystart;
    ud.replay.onend = @replayend;
    ud.replay.onturn = @i_replayparanalturn;
    set(hax, 'userdata', ud);
end
settings = [];
for i = 1:length(pars)
    settings.(pars{i}) =  evalin('base', char(pars{i}));
end
if new
    hold('on');
    nser = 0;
else
    nser = length(get(gca, 'children'));
end
%plotedit('off');
set(hfig,'Name','Function plot');
oldhold = ishold;
hold('on');
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
if length(y) == 1
    y1 = y(1);
    y = zeros(size(x)) + y1;
end
if length(x) == 1
    x1 = x(1);
    x = zeros(size(y)) + x1;
end
if any(~isreal(y))||any(~isreal(x))
    warning('grind:funplot:complex','Imaginary parts of complex values ignored');
    x=real(x);
    y=real(y);
end
h = plot(x, y, pen.pen, 'Color', pen.color);
i_plotdefaults(hfig);
xlabel(i_disptext(xaxis));
if nser > 1
    ylabel('');
else
    ylabel(i_disptext(yaxis));
end
i_grindlegend(1, h, yaxis, settings);
if ~isempty(axrange)
    set(gca, 'XLim', axrange);
end
% if isfield(g_grind, 'pen') && ~isempty(g_grind.pen)
%    nextpen;
% end
if ~isempty(ayrange)
    set(gca, 'YLim', i_checkstr(ayrange));
end
%plotedit('on');
if logxy(1)
    set(gca,'XScale','Log');
else
    set(gca,'XScale','Linear');
end
if logxy(2)
    set(gca,'YScale','Log');
else
    set(gca,'YScale','Linear');
end
i_grindlegend(11);
if ~oldhold
    hold('off');
end

function replaystart(hax)
global g_paranal;
if ishandle(hax)
    ud = get(hax, 'userdata');
    if isfield(g_paranal.run, 'parvalues')&&~isempty(g_paranal.run.parvalues)
        ud.replay.pars = g_paranal.run.parvalues(:, 1);
        ud.replay.oldpar = evalin('base', g_paranal.run.pars{1});
        tstep = ud.replay.pars(2) - ud.replay.pars(1);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(1,end)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',numel(g_paranal.run.t));
    else
        ud.replay.pars = {};
        ud.replay.oldpar = [];
        ud.replay.settings=struct('tvar',[],'tlim',[0 10],'tstep',1,'numt',1000);
    end
    set(hax, 'userdata', ud);
    i_figure(get(hax, 'parent'));
end
function replayend(hax, closedlg)
global g_paranal;
if ishandle(hax)&&closedlg
    i_figure(get(hax, 'parent'));
    ud = get(hax, 'userdata');
    if ~isempty(ud.replay.oldpar)
        assignin('base', g_paranal.run.pars{1}, ud.replay.oldpar);
        funplot('-u', hax);
    end
    set(hax, 'userdata', ud);
end
function p = replaycallback(hax,avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&isempty(avar)||~isempty(g_paranal.run)&&(isfield(g_paranal.run, 'parvalues')&&~isempty(g_paranal.run.parvalues)&&strcmp(avar, g_paranal.run.pars{1}))
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')&&~isempty(ud.replay.oldpar)
        ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
        [~, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
        p = g_paranal.run.parvalues(stepndx);
        [~, pnow] = min(abs(ud.replay.pars - p));
        if ~isempty(pnow)
            oldpar = evalin('base', g_paranal.run.pars{1});
            try
                assignin('base', g_paranal.run.pars{1}, ud.replay.pars(pnow));
                funplot('-u', hax);
                assignin('base', g_paranal.run.pars{1}, oldpar);
            catch err
                assignin('base', g_paranal.run.pars{1}, oldpar);
                rethrow(err);
            end
        end
    end
end


