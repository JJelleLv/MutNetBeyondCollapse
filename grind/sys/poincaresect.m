%POINCARESECT   Construct a Poincare section
%   All trajectories that cross a certain surface in the state variable space
%   are plotted.
%   POINCARESECT analyses the results of the last run. (if there is
%   no last run, it calls <a href="matlab:help time">time</a>). Use <a href="matlab:help ru">ru</a> if you want to
%   update the last run.
%
%   Usage:
%   POINCARESECT - opens a dialog to enter the variables.
%   POINCARESECT APLANE- analyses the plane APLANE. APLANE is a function that 
%   defines the plane (for instance X=9). Only increasing trajectories 
%   are analysed.
%   POINCARESECT APLANE 0 - Analyses only decreasing trajectories.
%   POINCARESECT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - variable or function for poincaremap.
%     'hax' [handle] - handle to axis, used by replayall
%     'increasing' [logical] - increasing or decreasing trajectories
%     'plane' [string] - function that defines the plane (for instance X=9)
%   POINCARESECT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-addreplaydata' - add replaydata (used internally)
%
%  
%   See also poincaremap, lorenzmap, takens, lyapunov
%
%   Reference page in Help browser:
%      <a href="matlab:commands('poincaresect')">commands poincaresect</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function poincar2 = poincaresect(varargin)
%(aplane, increasing)
global g_grind g_t;
fieldnams={'fun', 's', 'variable or function for poincaremap.',i_statevars_names(1);...
   'plane', 's', 'function that defines the plane (for instance X=9)',sprintf('%s==0',i_statevars_names(1));...
   'increasing', 'l', 'increasing or decreasing trajectories',true;...
   'hax', 'h', 'handle to axis, used by replayall',get(get(0,'CurrentFigure'),'CurrentAxes')}';
args=i_parseargs(fieldnams,'if(hasoption(''-addreplaydata'')),deffields=''hax'';else,deffields=''plane,increasing'';end;','-addreplaydata',varargin);
i_parcheck;
if ~isfield(args,'plane')&&~any(strcmp(args.opts,'-addreplaydata'))
   prompt={'Variable for mapping','Formula of the plane (e.g. x==4)','Increasing trajectories?(1/0)'};
   if isfield(g_grind,'poincare')
       answer=g_grind.poincare;
   else
       answer={i_statevars_names(1),sprintf('%s==0',i_statevars_names(1)),'1'};
   end
   answer=inputdlg(prompt,'Poincare Map',1,answer);
   args.fun=answer{1};
   args.plane=answer{2};
   args.increasing=i_checkstr(answer{3});
%   errordlg('Not enough arguments, Usage: poincaremap avar1 avar avalue');
%   error('GRIND:poincaremap:ArgError','Not enough arguments, Usage: poincaremap avar1 avar avalue');
end
if ~isfield(args,'increasing')
    args.increasing=true;
end
if ~isfield(args,'fun')
    args.fun=i_statevars_names(1);
end
if any(strcmp(args.opts,'-addreplaydata'))
    addreplaydata(args.hax)
    return;
end
g_grind.poincare = {args.fun, args.plane, int2str(args.increasing)}; % also used by poincaremap
t_poincare = i_poincare(args.plane, args.increasing);
if isempty(t_poincare)
    error('GRIND:poincaresect:intersection','Not any points intersect the plane');
end

% i = 1;
% if i == ivar
%    i = i + 1;
% end
% j = i + 1;
% if (j == ivar) && (j < g_grind.statevars.dim)
%    j = j + 1;
% end
% if j <= g_grind.statevars.dim
f=strfind(args.plane, '=');
if length(f) <= 2
    avar = args.plane(1:f - 1);
    avalue = str2double(args.plane(f + 2:end));
    if isnan(avalue)
        avar = '';
    end
else
    avar = '';
end

us.plane = args.plane;
us.increasing = args.increasing;
us.t_poincare = t_poincare;
XX = outfun(g_grind.xaxis.var);
XX = interp1(g_t, XX, t_poincare, 'pchip');
if strcmp(avar, g_grind.xaxis.var)
    XX = zeros(size(XX)) + avalue;
end
if ~isempty(g_grind.yaxis.var)
    YY = outfun(g_grind.yaxis.var);
    YY = interp1(g_t, YY, t_poincare, 'pchip');
    if strcmp(avar, g_grind.yaxis.var)
        YY = zeros(size(YY)) + avalue;
    end
else
    YY=[];
end
if ~isempty(g_grind.zaxis.var)
    ZZ = outfun(g_grind.zaxis.var);
    ZZ = interp1(g_t, ZZ, t_poincare, 'pchip');
    if strcmp(avar, g_grind.zaxis.var)
        ZZ = zeros(size(ZZ)) + avalue;
    end
else
    ZZ = [];
end
if ~isempty(YY)
    [H, isnew] =  i_makefig('poinsec');
    if isnew
        hold on;
    end
    if isempty(ZZ)
        h1 = plot(XX, YY, 'k.');
    else
        h1 = plot3(XX, YY, ZZ, 'k.');
    end
    if isnew
        i_plotdefaults;
    end
    if ~args.increasing
        set(h1, 'Color', [1 0 0])
    end
    if args.increasing
        set(h1,'displayname','increasing')
    else
        set(h1,'displayname','decreasing')
    end
    set(h1, 'MarkerSize', 5)
    ud=get(h1,'userdata');
    ud.tdata = t_poincare;
    ud.plane = args.plane;
    ud.xdata = XX;
    ud.ydata = YY;
    ud.zdata = ZZ;
    set(h1, 'userdata', ud);
    set(H,'name',['Poincare section through plane ' args.plane])
    set(H, 'userdata', us);
    htitle = title(['Poincare section S: ' args.plane ]);
    set(htitle,'fontweight','normal');
    xlabel(i_disptext(g_grind.xaxis.var));
    ylabel(i_disptext(g_grind.yaxis.var));
    if ~isempty(ZZ)
        zlabel(i_disptext(g_grind.zaxis.var));
    end
    ch=findobj(gca,'type','line');
    if length(ch) > 1
        legend(gca, 'show');
    end
    poincaresect('-addreplaydata',gca);
end
% else
%    error('GRIND:poincaresect:TooFewDims','Not enough dimensions to draw a 2D poincare section');
% end
if nargout > 0
    poincar2 = t_poincare;
end

function addreplaydata(hax)
global g_t;
ud = get(hax, 'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
ud.replay.hwhere  = -1;
if isempty(g_t)
    ud.replay.settings=struct('tvar','t','tlim',[0 g_grind.ndays],'numt',0);
else
    ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
end
set(hax,'userdata', ud);


function replaystart(hax)
global g_grind g_t;
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
        if g_grind.solver.isdiffer&&g_grind.statevars.dim==1 %for cobwebbing!
            ud.replay.settings.numt = ud.replay.settings.numt * 2 - 1;
        end
    end
    if ~ishandle(ud.replay.hwhere)
        hold(hax,'on');
        ud.replay.meantstep = mean(diff(g_t));
        ud.replay.hwhere = plot(hax,g_t(1), 0, 'ro');
        set(ud.replay.hwhere, 'MarkerFaceColor', [1 0 0]);
        set(ud.replay.hwhere,'tag','hwhere');
        set(ud.replay.hwhere,'visible','off');
        i_grindlegend(-1, ud.replay.hwhere);
    end
    set(hax, 'userdata', ud);
    if ~isoctave&&verLessThan('matlab','8.4.0')
        set(hax, 'drawmode','fast');
    else
        set(hax,'SortMethod','depth');
    end
    set(hax,'xlimmode','manual');
    set(hax,'ylimmode','manual');
end

function replayend(hax, closedlg)
if closedlg
    replaycallback(hax, 't', 1);
end
if ishandle(hax)
    ud=get(hax,'userdata');
    if ishandle(ud.replay.hwhere)
        delete(ud.replay.hwhere);
        ud.replay.hwhere=-1;
    end
    set(hax,'userdata',ud);
    i_figure(get(hax, 'parent'));
    set(hax,'xlimmode','auto');
    set(hax,'ylimmode','auto');
end

function t = replaycallback(hax,avar, relt)
t = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, 't')
    ud = get(hax, 'userdata');
    t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
    if isfield(ud, 'replay')
        ser = get(hax, 'children');
        for i = 1:length(ser)
            udser = get(ser(i), 'userdata');
            if isfield(udser, 'tdata')
                ndx=udser.tdata <= t;
                if ishandle(ud.replay.hwhere)
                    lastndx = find(ndx, 1, 'last');
                    if ~isempty(lastndx)%&&abs(t - ud.tdata(lastndx)) < ud.replay.meantstep
                        set(ud.replay.hwhere,'visible','on');
                        if isempty(udser.zdata)
                            set(ud.replay.hwhere, 'XData',  udser.xdata(lastndx), 'YData', udser.ydata(lastndx));
                        else
                            set(ud.replay.hwhere, 'XData',  udser.xdata(lastndx), 'YData', udser.ydata(lastndx), 'ZData', udser.zdata(lastndx));
                        end
                    end
                end
                if isempty(udser.zdata)
                    set(ser(i), 'XData',  udser.xdata(ndx), 'YData', udser.ydata(ndx));
                else
                    set(ser(i), 'XData',  udser.xdata(ndx), 'YData', udser.ydata(ndx), 'ZData', udser.zdata(ndx));
                end
            end
        end
    end
end




