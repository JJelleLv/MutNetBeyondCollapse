%REPLAYALL  Replay the model
%    Opens a dialog box that makes it possible to "replay" the model, i.e. 
%    the figures are redrawn as if the model runs. This is especially important 
%    for commands like <a href="matlab:help viewcells">viewcells</a> and <a href="matlab:help vectplot">vectplot</a> that are updated in time. 
%    All figures are updated simultaneously (including time graphs and
%    phase plane plots). It is also possible to save the animation as movie file (*.avi).
%    It is also possible to replay <a href="matlab:help paranal">paranal</a> sessions in combination with nullcline plots
%    (<a href="matlab:help null">null -p</a> or time plots <a href="matlab:help time">time -p</a>
%       
%    Usage:
%    REPLAYALL - opens the dialog screen.
%    REPLAYALL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'avifig' [handle] - handle of the figure for making a video.
%     'close' [logical] - close dialog after the run
%     'createavi' [logical] - set number or title of the figure to create an avi file (use in combination with 'createavi') 
%       (default = current figure).
%     'delay' [number>=0] - set the delay of the replay (s) (default = 0)
%     'forward' [logical] - set to 1 for forward replay and to 0 for backward replay (applies only
%       to replaying paranal in two directions) (default is current direction).
%     'nodialog' [logical] - suppress the replayall dialog (default = 0)
%     'numsteps' [integer>=0] - set the number of steps (instead of the step size)
%     'play' [logical] - set to 1 to replay now (default=0)
%     'relvalue' [number>=0 and number<=1] - set the relative start value of the variable as a value between
%       0=min and 1=max. You can also specify the range for replaying ([start end]). (default [current 1]).
%     'step' [number>0] - set the step size of replay (units of variable) (default depends)
%     'value' [number or start | end] - set the start value of the variable, either a value of the variable.
%       You can also specify the range for replaying ([start end]) or 'start' for the start 
%       value or 'end' for the last value (default='last');
%     'variable' [string] - set the name of the variable to change (default = 't')
%     'videowriter' [general] - set to a VideoWriter object to write the frames to a file. see
%       <a href="matlab:help videowriter">VideoWriter</a> for options).
%   REPLAYALL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-pl' HAX - add replaydata to a plot.
%     '-removeinfo' - remove replaydata from all figures.
%         
%      Examples:
%      replayall('value','start','variable','t','play',1) - replay the time plots from the start
%      replayall('value',1000,'variable','t','close',1) - set the time to 1000 and close the replayall dialog
%       
%       
%    See also vectplot, viewcells, time, phas, null, funplot, <a href="matlab:help videowriter">VideoWriter</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('replayall')">commands replayall</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function res = replayall(varargin)
global g_grind;
% if isempty(which('i_replayall'))
%     addpath([grindpath filesep 'sys2']);
% end
if nargin < 1            
    if ~isempty(g_grind)&&~isempty(g_grind.figopts)
        opts=g_grind.figopts;
    else
        opts={};
    end
    i_replayalldlg('init','',opts);
    args.opts={};
else
    fieldnams={'delay', 'n>=0', 'set the delay of the replay (s) (default = 0)',0;...
        'step', 'n>0', 'set the step size of replay (units of variable) (default depends)',1;...
        'nodialog', 'l', 'suppress the replayall dialog (default = 0)',false;...
        'numsteps', 'i>=0', 'set the number of steps (instead of the step size)',100;...
        'variable', 's', 'set the name of the variable to change (default = ''t'')','t';...
        'value', 'n#e[start|end]', 'set the start value of the variable, either a value of the variable.','start';...
        'relvalue', 'n>=0&n<=1', 'set the relative start value of the variable as a value between',1;...
        'forward', 'l', 'set to 1 for forward replay and to 0 for backward replay (applies only',true;...
        'play', 'l', 'set to 1 to replay now (default=0)',false;...
        'close', 'l', 'close dialog after the run',false;...
        'createavi', 'l', 'set number or title of the figure to create an avi file (use in combination with ''createavi'')',1;...
        'avifig', 'h', 'handle of the figure for making a video.',get(0,'CurrentFigure');...
        'videowriter', '', 'set to a VideoWriter object to write the frames to a file. see',[]}';
    args=i_parseargs(fieldnams,'variable','-removeinfo,-pl',varargin);
end
if any(strcmp(args.opts,'-removeinfo'))
    hfigs=findall(0,'type','figure');
    tagndx=strcmp(get(hfigs,'tag'),'i_replayalldlg');
    delete(hfigs(tagndx));
    hfigs(tagndx) = [];
    haxes=findobj(hfigs,'type','axes');
    tags = get(haxes, 'tag');
    if ~isempty(tags)
        haxes=haxes( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
        for i  = length(haxes):-1:1
            ud = get(haxes(i), 'userdata');
            if isfield(ud, 'replay')
                ud = rmfield(ud, 'replay');
                set(haxes(i), 'userdata', ud);
            end
        end
    end
    return;
end
if any(strcmp(args.opts,'-pl'))
    if nargin < 2
        par1 = 't';
    else
        par1 = varargin{2};
    end
    if nargin < 3
        par2 = gca;
    else
        par2 = varargin{3};
    end
    replayplot(par2, par1);
end
%'PropertyName', PropertyValue pairs
%delay = the delay in replay (s)
%step = the stepsize of replay (units of variable)
%variable = the name of the variable to change
%value = set the value of the variable
%play = set to 1 to replay now
%close = set to 1 to close the replay now (cannot be combined with
%run);
%   command = struct(varargin{:});

fields = fieldnames(args);
%    for i = 1:length(fields)
%       if ~any(strcmpi(fields{i},{'delay','step','nodialog','numsteps','variable','value','relvalue','forward','play','close','createavi','avifig','videowriter'}))
%          error('grind:replayall','Unknown property of replayall "%s"',fields{i});
%       end
%    end
figs=findall(0,'type','figure');
tagndx=strcmp(get(figs,'tag'),'i_replayalldlg');
if ~any(tagndx)
    if nargout > 0
        res = false;
    end
    if length(fields)==1&&strcmp(fields, 'close')
        return;
    end
    if ~isempty(g_grind)&&~isempty(g_grind.figopts)
        opts=g_grind.figopts;
    else
        opts={};
    end
    hfig = i_replayalldlg('init','', opts);
    drawnow;
else
    if nargout > 0
        res = true;
    end
    hfig = figs(tagndx);
    i_figure(hfig);
end
if any(strcmp(fields, 'delay'))
    h=findobj(hfig, 'Tag','DelayEdit');
    set(h,'string',num2str(args.(fields{strcmpi(fields,'delay')})));
end
if any(strcmp(fields, 'step'))
    h=findobj(hfig, 'Tag','StepEdit');
    v = args.(fields{strcmpi(fields, 'step')});
    if ischar(v)
        set(h, 'string', v);
    else
        set(h, 'string', num2str(v));
    end
end
if any(strcmp(fields, 'numsteps'))
    h=findobj(hfig, 'Tag','StepEdit');
    ud = get(hfig, 'userdata');
    v = args.(fields{strcmpi(fields, 'numsteps')});
    if ischar(v)
        v = str2double(v);
    end
    if isnan(v)
        v = 1000;
    end
    v = (ud.tlim(2) - ud.tlim(1)) / (v + 1);
    set(h, 'string', num2str(v));
end
if any(strcmp(fields, 'variable'))
    h=findobj(hfig,'Tag','VarList');
    s = get(h, 'string');
    v = args.(fields{strcmpi(fields, 'variable')});
    f = strcmp(v, s);
    if ~any(f)
        warning('replayall:variable','variable "%s" cannot be replayed',v);
        s = [s {v}];
        set(h, 'string', s)
        f = strcmp(v, s);
    end
    set(h, 'value', find(f));
end
if any(strcmp(fields, 'forward'))
    v = args.(fields{strcmpi(fields, 'forward')});
    if ischar(v)
        v=any(strcmpi(v,{'1','on','true'}));
    end
    h1 = i_replayalldlg('funhandles');
    h1.doturn(hfig, [], v);
end
valueset = false;
if any(strcmp(fields, 'relvalue'))
    v = args.(fields{strcmpi(fields, 'relvalue')});
    valueset = true;
    if ischar(v)
        if strcmp(v, 'end')
            v = 1;
        elseif strcmp(v, 'start')
            v = 0;
        else
            v = evalin('base', v);
        end
    end
    if isempty(v)
        v = 0;
    end
    if v(1) < 0
        v(1) = 0;
        warning('replayall:relvalue','relvalue too low, set the value to 0');
    end
    if v(end) > 1
        v(end) = 1;
        warning('replayall:relvalue','relvalue too high, set the value to 1');
    end
    h=findobj(hfig, 'Tag','TheSlider');
    set(h, 'value', v(1));
    if length(v) == 2
        ud = get(hfig, 'userdata');
        ud.endvalue = v(end);
        set(hfig, 'userdata', ud);
    end
    callb = get(h, 'callback');
    if ~isempty(callb)
        callb(h);
    end
end

if any(strcmpi(fields, 'value'))
    v = args.(fields{strcmpi(fields, 'value')});
    valueset = true;
    if ischar(v)&&any(strcmpi(v,{'end','start'}))
        h=findobj(hfig, 'Tag','TheSlider');
        if strcmpi(v, 'end')
            set(h, 'value', 1)
        else
            set(h, 'value', 0);
        end
        callb = get(h, 'callback');
        if ~isempty(callb)
            callb(h);
        end
    else
        h=findobj(hfig, 'Tag','TimeEdit');
        if ischar(v)
            evalin('base', v);
        end
        set(h, 'string', num2str(v(1)));
        if length(v) == 2
            h1 = i_replayalldlg('funhandles');
            ud = get(hfig, 'userdata');
            [~, ~, fullud] = h1.makeconnections(hfig,0);
            if isempty(fullud)
                ud.endvalue = v(2);
            else
                ud.endvalue = (v(2) - fullud.settings.tlim(1)) / (fullud.settings.tlim(2) - fullud.settings.tlim(1));
            end
            set(hfig, 'userdata', ud);
        end
        callb = get(h, 'callback');
        if ~isempty(callb)
            callb(h);
        end
    end
end

if any(strcmpi(fields, 'play'))
    v = args.(fields{strcmpi(fields, 'play')});
    if ischar(v)
        v=any(strcmpi(v,{'1','on','true'}));
    end
    %    ud=get(hfig,'userdata');
    if v
        if any(strcmpi(fields, 'nodialog'))&& args.(fields{strcmpi(fields, 'nodialog')})==0
            movegui(hfig, 'onscreen');
            set(hfig,'visible','on');
            i_figure(hfig);
            set(hfig,'visible','off');
        end
        pushbutton('playButton');
    end
end
if any(strcmpi(fields, 'createavi'))
    v = args.(fields{strcmpi(fields, 'createavi')});
    if ischar(v)
        v=any(strcmpi(v,{'1','on','true'}));
    end
    if v
        havifig=findobj(0,'type','figure');
        f=strcmp(get(havifig,'Tag'),'i_replayalldlg');
        if any(f)
            havifig = havifig(~f);
        end
        if isempty(havifig)
            error('replayall:nofigure','No figure to replay for creating movie');
        end
        if (length(havifig) > 1)
            if any(strcmpi(fields, 'avifig'))
                avifig = args.(fields{strcmpi(fields, 'avifig')});
                if ~ischar(avifig)
                    avifig = num2str(avifig);
                end
                titles = get(havifig, 'name');
                f = strcmp(avifig, titles);
                if any(f)
                    havifig = havifig(f);
                    havifig = havifig(1);
                else
                    titles = str2double(titles);
                    f=titles == havifig;
                    if any(f)
                        havifig = havifig(f);
                        havifig = havifig(1);
                    else
                        error('replayall:avifig','Figure "%s" not found',avifig)
                    end
                end
            end
            if length(havifig) > 1
                havifig =  get(0, 'currentfig');
            end
        end
    end
    if v
        %             movegui(hfig,'onscreen');
        %             set(hfig,'visible','on');
        %             figure(hfig);
        if valueset
            h=findobj(hfig, 'Tag','TheSlider');
            ud = get(hfig, 'userdata');
            rangev = [get(h, 'value') ud.endvalue];
        else
            rangev = [0 1];
        end
        h1 = i_replayalldlg('funhandles');
        mov = h1.makeavi(hfig, havifig, rangev);
    end
    if any(strcmpi(fields, 'videowriter'))
        vidObj = args.(fields{strcmpi(fields, 'videowriter')});
    else
        vidObj = VideoWriter('video.avi');
        warning('replayall:videowriter','No "videowriter" object specified: writing to "video.avi" using default settings')
    end
    h1.writevideoobj(vidObj, mov);
else
    
end

if any(strcmpi(fields, 'close'))
    v = args.(fields{strcmpi(fields, 'close')});
    if ischar(v)
        v=any(strcmpi(v,{'1','on','true'}));
    end
    if v
        h=findobj(hfig, 'Tag','TheSlider');
        h1 = i_replayalldlg('funhandles');
        h1.closeconnections(h, 2);
        %           ud =  get(hfig, 'userdata');
        %           if isfield(ud, 'hListener')&&ishandle(ud.hListener)
        %              delete(ud.hListener);
        %           end
        delete(hfig)
    end
end

function pushbutton(btn)
figs=findall(0,'type','figure');
tagndx=strcmp(get(figs,'tag'),'i_replayalldlg');
if any(tagndx)
    h1 = figs(find(tagndx, 1));
else
    h1 = i_replayalldlg('init');
end
%set(h1,'visible','off');
h = findobj(h1, 'tag', btn);
callb = get(h, 'callback');
if ~isempty(callb)
    callb(h);
end

function replayplot(H1, avar)
ud = get(H1, 'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
data = figdata('-all', H1);
ud.replay.xdata = data(:, 1);
xlim([min(ud.replay.xdata) max(ud.replay.xdata)]);
ud.replay.ydata = data(:, 2:end);
ud.replay.settings=struct('tvar',avar,'tlim',[min(ud.replay.xdata) max(ud.replay.xdata)],'numt',length(ud.replay.xdata));
set(H1,'userdata', ud);

function replaystart(hax)
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    if ~isoctave&&verLessThan('matlab','8.4.0')
        set(hax, 'drawmode','fast');
    else
        set(hax,'SortMethod','depth');
    end
end

function replayend(hax, closedlg)
ud = get(hax, 'userdata');
if closedlg
    replaycallback(hax, ud.replay.settings.tvar, 1);
end
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    set(hax,'xlimmode','auto');
    set(hax, 'userdata', ud);
end

function t = replaycallback(hax,avar, relt)
ud = get(hax, 'userdata');
t = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar,  ud.replay.settings.tvar)
    t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
    if isfield(ud, 'replay')
        ser = get(hax, 'children');
        tt = ud.replay.xdata;
        xx=tt(tt <= t);
        for i = 1:length(ser)
            if ~isempty(tt)
                set(ser(i), 'XData', xx, 'YData', ud.replay.ydata(tt <= t,i));
            end
        end
    end
end

