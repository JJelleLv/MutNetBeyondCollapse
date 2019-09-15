function res = i_replayalldlg(flag, par,opts)
if nargin < 1
    flag = 'init';
end

switch flag
    case 'funhandles'
        res.closeconnections = @closeconnections;
        res.makeconnections = @makeconnections;
        res.doturn = @doturn;
        res.makeavi = @makeavi;
        res.writevideoobj =  @writevideoobj;
    case 'init'
        if nargin < 2
            par = '';
        end
        
        if nargin < 3
            opts = {};
        end
        
        [~, ~, fullud] = makeconnections(-1,1);
        if isempty(fullud)||fullud.settings.numt==1
            SliderStep = [0.001 0.1];
            StepEdit = 1;
            tlim = [0 1];
        else
            SliderStep = [1 / (fullud.settings.numt - 1) 0.1];
            StepEdit = (fullud.settings.tlim(2) - fullud.settings.tlim(1)) / (fullud.settings.numt - 1);
            tlim = fullud.settings.tlim;
            if isnan(StepEdit)
                StepEdit = 1;
            end
            
        end
        
        figs=findall(0,'type','figure');
        tagndx=strcmp(get(figs,'tag'),'i_replayalldlg');
        if any(tagndx)
            hfig = figs(find(tagndx, 1));
            i_figure(hfig,opts{:});
            ud = get(hfig, 'userdata');
            ud.tlim = tlim;
            if nargin > 1
                h=findobj(hfig,'tag','VarList');
                varlist = get(h, 'string');
                f = find(strcmp(par, varlist));
                if ~isempty(f)
                    set(h, 'value', f(1));
                else
                    varlist = [varlist {par}];
                    set(h, 'string', varlist);
                    set(h, 'value', length(varlist));
                end
                
            end
            
            h=findobj(hfig,  'Tag','StepEdit');
            set(h, 'String', num2str(StepEdit));
            h=findobj(hfig,'Tag','TheSlider');
            set(h, 'SliderStep', SliderStep);
            movegui(hfig, 'onscreen');
            if ~any(strcmpi(opts,'visible'))
                set(hfig,'visible','on');
            end
            
        else
            ud = struct('active', 0,'endvalue',1,'playbutton','Play >','tlim',tlim);
            hfig = figure( ...
                'WindowStyle','Normal',...
                'PaperUnits',get(0,'defaultfigurePaperUnits'),...
                'Color', [0.941176470588235 0.941176470588235 0.941176470588235], ...
                'Colormap',get(0,'defaultfigureColormap'),...
                'IntegerHandle','off',...
                'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
                'MenuBar','none',...
                'Name','Replay Settings',...
                'NumberTitle','off',...
                'PaperPosition',get(0,'defaultfigurePaperPosition'),...
                'PaperSize',get(0,'defaultfigurePaperSize'),...
                'PaperType',get(0,'defaultfigurePaperType'),...
                'Position', [544 518 469 107], ...
                'HandleVisibility','on',...
                'ResizeFcn',@resizereplay,...
                'Visible','off',...
                'CloseRequestFcn', @closebtn_Callback, ...
                'Tag','i_replayalldlg');
            try
                p1 = get(0, 'screensize');
                p2 = get(0, 'monitorpositions');
                if (numel(p1)==numel(p2))&&all(p1==p2)
                    movegui(hfig, 'south');
                    movegui(hfig, 'onscreen');
                end
                
                if ~isempty(opts)
                    set(hfig,opts{:});
                else
                    set(hfig,'visible','on')
                end
                
            catch
                set(hfig,'visible','on')
            end
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [15 83 61 18], ...
                'String','Delay (s):',...
                'Style','text',...
                'Tag','DelayText');
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [190 83 61 18], ...
                'String','Replay step:',...
                'Style','text',...
                'Tag','ReplayText');
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'BackgroundColor', [1 1 1], ...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [70 80 60 27], ...
                'String','0',...
                'Style','edit',...
                'Tag','DelayEdit');
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'BackgroundColor', [1 1 1], ...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [260 80 60 27], ...
                'String', num2str(StepEdit), ...
                'Style','edit',...
                'Tag','StepEdit');
            uicontrol( 'Parent', hfig, ...
                'FontUnits','pixels',...
                'Callback', @stepbtn_Callback, ...
                'FontSize', 10.6666666666667, ...
                'Position', [321 80 20 27], ...
                'String','...',...
                'TooltipString','Edit as number of steps',...
                'Tag','stepButton' );
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'Callback', @closebtn_Callback, ...
                'FontSize', 10.6666666666667, ...
                'Position', [400 83 70 27], ...
                'String','Close',...
                'Tag','closeButton' );
            uicontrol('Parent', hfig, ...
                'Units','pixels',...
                'FontUnits','pixels',...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [15 40 123 22 ], ...
                'String','Inactive',...
                'Style','text',...
                'Tag','CurrentText');
            uicontrol('Parent', hfig, ...
                'Units','pixels',...
                'Callback', @varlist_Callback, ...
                'FontUnits','pixels',...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Style','popupmenu',...
                'String', {par}, ...
                'Value', 1, ...
                'Position', [190 46 50 22], ...
                'Tag','VarList');
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'Callback', @timeedit_Callback, ...
                'BackgroundColor', [1 1 1], ...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [260 46 60 27], ...
                'String','0',...
                'Style','edit',...
                'Tag','TimeEdit');
            
            uicontrol( 'Parent', hfig, ...
                'FontUnits','pixels',...
                'Callback', @avibtn_Callback, ...
                'FontSize', 10.6666666666667, ...
                'Position', [400 50 70 27], ...
                'String','Create AVI',...
                'Tag','aviButton' );
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'BackgroundColor', [0.9 0.9 0.9], ...
                'Callback', @slider_Callback, ...
                'FontSize', 10.6666666666667, ...
                'Position', [15 7 305 27], ...
                'Max', 1, ...
                'Min', 0, ...
                'Value', 0, ...
                'SliderStep',SliderStep, ...
                'String',{  'Slider' },...
                'Style','slider',...
                'Tag','TheSlider' );
            
            %          ud.hListener = addlistener(hSlider,'Value','PostSet',@slider_Listener);
            
            uicontrol('Parent', hfig, ...
                'FontUnits','pixels',...
                'Callback', @playbtn_Callback, ...
                'FontSize', 10.6666666666667, ...
                'Position', [331 7 70 28], ...
                'String', ud.playbutton, ...
                'TooltipString','Play',...
                'Tag','playButton' );
            
            uicontrol( 'Parent', hfig, ...
                'Units','pixels',...
                'FontUnits','pixels',...
                'FontSize', 10.6666666666667, ...
                'HorizontalAlignment','left',...
                'Position', [140 40 50 22], ...
                'String','Variable:',...
                'Tag','textvar',...
                'Style','text');
            
            %       i_figure(hfig);
            makeconnections(hfig, 0);
            drawnow;
        end
        
        set(hfig, 'userdata', ud);
        if nargout>0
            res=hfig;
        end
        
end



function avibtn_Callback(hobject, ~)
hdlg = getparentfig(hobject);
hfigs=findobj(0,'type','figure');
hfigs=hfigs(~strcmp(get(hfigs,'Tag'),'i_replayalldlg'));
hfigs = sort(hfigs);
if length(hfigs) > 1
    titles = get(hfigs, 'name');
    numtitles = strcmp(get(hfigs,'NumberTitle'),'on');
    for i = 1:length(hfigs)
        if isa(hfigs, 'matlab.ui.Figure')
            no = get(hfigs(i), 'Number');
        else
            no = hfigs(i);
        end
        if numtitles(i)
            titles{i} = sprintf('Figure %d: %s',no , titles{i});
        end
        
    end
    
    s = listdlg('PromptString','Select a figure:',...
        'SelectionMode','single',...
        'ListString', titles);
    hfigs = hfigs(s);
end

if ~isempty(hfigs)
    mov = makeavi(hdlg, hfigs);
    [vidobj, mov] = i_videodlg(mov);
    if ~isempty(vidobj)
        writevideoobj(vidobj, mov);
    end
    vidobj=1;
    while ~isempty(vidobj)&&strcmp(questdlg('Do you want to make another video file with these frames?', ...
            'Video Question', ...
            'Yes', 'No', 'No'),'Yes')
        [vidobj, mov] = i_videodlg(mov);
        if ~isempty(vidobj)
            writevideoobj(vidobj, mov);
        end
    end
    
else
    disp('No figures to convert');
end

function resizereplay(hobj,~)
hfig=getparentfig(hobj);
data={...
    'DelayText',[1 0 0 1],[10.5 61.5 45.75 13.5];...
    'ReplayText',[1 0 0 1],[141.75 61.5 45.75 13.5];...
    'DelayEdit',[1 0 0 1],[51.75 59.25 45 20.25];...
    'StepEdit',[1 0 0 1],[194.25 59.25 45 20.25];...
    'stepButton',[1 0 0 1],[240 59.25 15 20.25];...
    'closeButton',[0 1 0 1],[299.25 61.5 52.5 20.25];...
    'CurrentText',[1 0 0 1],[10.5 29.25 92.25 16.5];...
    'VarList',[1 0 0 1],[141.75 33.75 37.5 16.5];...
    'TimeEdit',[1 0 0 1],[194.25 33.75 45 20.25];...
    'aviButton',[0 1 0 1],[299.25 36.75 52.5 20.25];...
    'TheSlider',[1 1 0 1],[10.5 4.5 228.75 20.25];...
    'playButton',[0 1 0 1],[247.5 4.5 52.5 21];...
    'textvar',[1 0 0 1],[104.25 29.25 37.5 16.5];...
    };
i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[336 27.75 351.75 80.25]);


function mov = makeavi(hdlg, hfig, relfromto)
%global g_t;

% create AVI object

%create movie
if nargin < 3
    relfromto = [0 1];
end

i_figure(hfig)
haxes=findobj(hfig,'type','axes');
if ~isempty(haxes)
    tags = get(haxes, 'tag');
    haxes=haxes( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
    uds = get(haxes, 'userdata');
else
    uds = {};
end

if ~iscell(uds)
    uds = {uds};
end

k = 1;
for i = 1:length(uds)
    ud = uds{i};
    if isfield(ud, 'replay')
        replayopt{k} = ud.replay;
        hs(k) = haxes(i);
        k = k + 1;
        h=findobj(hdlg,'Tag','StepEdit');
        step = str2double(get(h, 'string'));
        if isnan(step)
            step = 1;
        end
        
        numsteps = round(abs((ud.replay.settings.tlim(2) - ud.replay.settings.tlim(1)) / step)) + 1;
        avar = ud.replay.settings.tvar;
        tt = linspace(relfromto(1), relfromto(2), round(numsteps .* (relfromto(2) - relfromto(1))));
    end
    
end

if k == 1
    error('grind:replayall','Cannot make a movie of this figure');
end

mov(1:length(tt)) = struct('cdata', [], ...
    'colormap', []);
oldcol = get(hfig, 'color');
set(hfig, 'Color','white')
try
    set(hdlg,'Visible','off');
    for i = 1:length(tt)
        updateall(hdlg, hs, replayopt, avar, tt(i));
        mov(i) = getframe(hfig);
    end
    
    set(hfig, 'color', oldcol);
    set(hdlg,'Visible','on');
catch err
    set(hdlg,'Visible','on');
    set(hfig, 'color', oldcol);
    rethrow(err);
end

drawnow;

function writevideoobj(vidobj, frames)
i_waitbar(1,100,'writeing video file','writing (cannot show progress)');
i_waitbar(1);
open(vidobj);
i_waitbar(1);
writeVideo(vidobj, frames);
close(vidobj);
i_waitbar([]);
function setstatustext(hdlg, astatus, at, fullud)
if ~isempty(at)
    tt = at;
    if ~isempty(fullud)
        tt = fullud.settings.tlim(1) + at * (fullud.settings.tlim(2) - fullud.settings.tlim(1));
        if isfield(fullud.settings, 'tstep')
            tt = fullud.settings.tlim(1) + floor((tt - fullud.settings.tlim(1)) / fullud.settings.tstep) * fullud.settings.tstep;
        end
        
        %         atext=sprintf('%s, %s=%g',astatus,fullud.settings.tvar,tt);
        %     else
        %         atext=sprintf('Not connected, slider value =%g',tt);
    end
    
    h=findobj(hdlg,'tag','TimeEdit');
    set(h, 'string', num2str(tt));
end

if ~isempty(astatus)
    h=findobj(hdlg,'tag','CurrentText');
    set(h, 'string', astatus);
end


function stepbtn_Callback(hObject, ~)
hdlg = get(hObject, 'parent');
[~, ~, fullud] = makeconnections(hdlg,0);
if ~isempty(fullud)
    h=findobj(hdlg,'tag','StepEdit');
    thestep = str2double(get(h, 'string'));
    numsteps = round(abs((fullud.settings.tlim(2) - fullud.settings.tlim(1)) / thestep)) + 1;
    if ~isfield(fullud.settings, 'ndata')
        fullud.settings.ndata = fullud.settings.numt;
    end
    
    prompt={sprintf('Number of steps for replay (n data=%g):', fullud.settings.ndata)};
    name = 'Enter the number of steps';
    answer = {num2str(numsteps)};
    answer = inputdlg(prompt, name, 1, answer);
    if ~isempty(answer)
        numsteps = str2double(answer{1});
        sgn = sign(thestep);
        if isnan(sgn)
            sgn = 1;
        end
        
        thestep = sgn * abs(fullud.settings.tlim(2) - fullud.settings.tlim(1)) / numsteps;
        set(h, 'string', num2str(thestep));
    end
    
end


function closebtn_Callback(hObject, ~)
try
    closeconnections(hObject, 1);
    disp('<a href="matlab:replayall">replayall</a> opens the replay dialog again');
    if ~strcmpi(get(hObject,'type'),'figure')
        hObject = get(hObject, 'parent');
    end
    
    if ishandle(hObject)
        %  ud =  get(hObject, 'userdata');
        %         if isfield(ud, 'hListener')&&ishandle(ud.hListener)
        %             delete(ud.hListener);
        %         end
        
        delete(hObject);
    end
    
catch err %in case of an error close the dialog anyway
    delete(hObject);
    rethrow(err);
end


function varlist_Callback(hObject, ~, ~)
hdlg = get(hObject, 'parent');
ud = get(hdlg, 'userdata');
wasactive = 0;
if ud.active
    wasactive = 1;
    ud.active = 0;
    set(hdlg,'userdata',ud);
end

[~, ~, fullud] = makeconnections(hdlg,0);
if ~isempty(fullud)
    StepEdit = (fullud.settings.tlim(2) - fullud.settings.tlim(1)) / (fullud.settings.numt - 1);
    h=findobj(hdlg,'tag','StepEdit');
    set(h, 'string', num2str(StepEdit));
end

h=findobj(hdlg,'tag','TheSlider');
at = get(h, 'value');
setstatustext(hdlg, '', at, fullud);
if wasactive
    h=findobj(hdlg,'tag','playButton');
    callb = get(h, 'callback');
    callb(h);
end

function timeedit_Callback(hObject, ~, ~)
hdlg = get(hObject, 'parent');
[~, ~, ud] = makeconnections(hdlg,0);
tt = str2double(get(hObject, 'string'));
if ~isnan(tt)
    if isempty(ud)
        at = tt;
    else
        at = (tt - ud.settings.tlim(1)) / (ud.settings.tlim(2) - ud.settings.tlim(1));
    end
    
    if at < 0
        if ~isempty(ud)
            warning('replayall:value','Value too %s, set the value to %g\n',iif(ud.settings.tlim(1)>ud.settings.tlim(2),'high','low'),ud.settings.tlim(1))
        end
        
        at = 0;
    end
    
    if at > 1
        if ~isempty(ud)
            warning('replayall:value','Value too %s, set the value to %g\n',iif(ud.settings.tlim(1)>ud.settings.tlim(2),'low','high'),ud.settings.tlim(2))
        end
        
        at = 1;
    end
    
    hslider=findobj(hdlg,'tag','TheSlider');
    set(hslider, 'value', at);
    slider_Callback(hslider);
end


function slider_Callback(hObject, ~, ~)
hdlg = get(hObject, 'parent');
ud = get(hdlg, 'userdata');
if ~ud.active
    hslider=findobj(hdlg,'tag','TheSlider');
    at = get(hslider, 'value');
    [handles, replayopt, fullud] = makeconnections(hdlg,0);
    for i = 1:length(handles)
        hfig = get(handles(i), 'parent');
        i_figure(hfig);
    end
    
    setstatustext(hdlg, 'Inactive', at, fullud);
    aVar = '';
    if ~isempty(fullud)
        aVar = fullud.settings.tvar;
    end
    
    updateall(hdlg, handles, replayopt,aVar, at)
    % else
    % set(hslider,'value',at);
end

function slider_Listener(~, ~, ~)
figs = allchild(0);
tagndx=find(strcmp(get(figs,'Tag'),'i_replayalldlg'),1);
hdlg = figs(tagndx);
ud = get(hdlg, 'userdata');
if ~ud.active
    hslider=findobj(hdlg,'tag','TheSlider');
    at = get(hslider, 'value');
    [handles, replayopt, fullud] = makeconnections(hdlg,0);
    setstatustext(hdlg, 'Inactive', at, fullud);
    
    for j = 1:length(handles)
        if ~isempty(replayopt{j}.callback)
            replayopt{j}.callback(handles(j), at);
        end
        
    end
    
    %    updateall(hdlg, handles, replayopt, at, atext)
    % else
    % set(hslider,'value',at);
end




function playbtn_Callback(hObject, ~, ~)
hdlg = get(hObject, 'parent');
ud = get(hdlg, 'userdata');
if ud.active
    closeconnections(hObject, 0);
    return;
end

ud.active = 1;
set(hdlg, 'userdata', ud);
[handles, replayopt, fullud] = makeconnections(hdlg,1);
hslider=findobj(hdlg,'tag','TheSlider');
tstart = get(hslider, 'value');
if ~isempty(fullud)&&fullud.settings.numt~=1
    set(hslider, 'SliderStep', [1 / (fullud.settings.numt - 1) 0.1]);
    h=findobj(hdlg,'Tag','playButton');
    ud = get(hdlg, 'userdata');
    ud.playbutton = get(h, 'string');
    set(hdlg, 'userdata', ud);
    set(h,'string','[Stop]');
    setstatustext(hdlg, 'Playing', [], fullud);
    at = tstart;
    while ud.active&&(at>=0) &&(at < ud.endvalue)
        if ~ishandle(hdlg)
            return;
        end
        
        h=findobj(hdlg,'Tag','DelayEdit');
        delay = str2double(get(h, 'string'));
        h=findobj(hdlg,'Tag','StepEdit');
        step = str2double(get(h, 'string'));
        if isnan(step)
            step = 1;
        end
        
        step = step / (fullud.settings.tlim(2) - fullud.settings.tlim(1));
        at = get(hslider, 'Value') + step;
        if at > ud.endvalue
            at = ud.endvalue;
        end
        
        updateall(hdlg, handles, replayopt, fullud.settings.tvar, at);
        if ~isnan(delay)&&(delay > 0)
            pause(delay);
        else
            drawnow;
        end
        
        if ishandle(hdlg)
            ud = get(hdlg, 'userdata');
        else
            ud.active = 0;
        end
        
    end
    
    if ud.active
        doturn(hdlg, fullud);
    end
    
else
    set(hslider, 'Value', 1);
    h=findobj(hdlg,'tag','TimeEdit');
    set(h,'String','1');
end

if ishandle(hdlg)
    ud.endvalue = 1;
    set(hdlg, 'userdata', ud);
end

closeconnections(hObject, 0);

function doturn(hdlg, fullud, goforward)
if nargin==1||isempty(fullud)
    [~, ~, fullud] = makeconnections(-1,1);
end

if nargin > 2&&isfield(fullud, 'settings')
    if goforward&&fullud.settings.tlim(1) < fullud.settings.tlim(2)
        return
    elseif ~goforward&&fullud.settings.tlim(1) > fullud.settings.tlim(2)
        return
    end
    
end

if isfield(fullud,'onturn')&&isa(fullud.onturn,'function_handle') %it should be a function_handle
    [ok, btntext] = fullud.onturn();
    if ok
        makeconnections(hdlg, 1);
        h=findobj(hdlg,'Tag','StepEdit');
        aval = str2double(get(h, 'string'));
        set(h, 'string', num2str(-aval));
        hslider=findobj(hdlg,'tag','TheSlider');
        set(hslider, 'value',0);
        h=findobj(hdlg,'Tag','playButton');
        set(h, 'string', btntext);
    end
    
end


function [handles, replayopt, fullud] = makeconnections(hdlg,runonstart)
handles = getallaxes;
replayopt = cell(size(handles));
fullud = [];
onturn = [];
varlist = replayopt;
for i  = length(handles):-1:1
    replayopt{i} = [];
    if ishandle(handles(i))
        ud = get(handles(i), 'userdata');
    else
        ud = [];
    end
    
    if isfield(ud, 'replay')
        if runonstart&&~isempty(ud.replay.onstart)
            %     try
            ud.replay.onstart(handles(i));
            %           catch
            %           end
            ud = get(handles(i), 'userdata');
        end
        
        replayopt{i} = ud.replay;
        if isfield(ud.replay, 'settings')
            fullud = ud.replay;
            varlist{i} = ud.replay.settings.tvar;
        end
        
        if isfield(ud.replay, 'onturn')
            onturn = ud.replay.onturn;
        end
        
    else
        handles(i) = [];
        replayopt(i) = [];
        varlist(i) = [];
    end
    
end

if ~isempty(fullud)
    fullud.onturn = onturn;
end

if ishandle(hdlg)
    if isempty(varlist)
        varlist = {'slider'};
    else
        ndx = zeros(size(varlist));
        for i = 1:length(varlist)
            ndx(i) = ~isempty(varlist{i});
        end
        
        varlist = varlist(logical(ndx));
        if isempty(varlist)
            varlist = {'slider'};
        end
        
    end
    
    uvarlist = unique(varlist);
    h=findobj(hdlg(1),'tag','VarList');
    aval = get(h, 'value');
    lst = get(h, 'String');
    if (aval > 0)&&aval<=length(lst)
        oldvar = lst{aval};
    else
        oldvar = '';
    end
    
    j = find(strcmp(oldvar, uvarlist));
    if isempty(j)
        maxn = -1;
        j = 1;
        for i = 1:length(uvarlist)
            n = sum(strcmp(uvarlist{i}, varlist));
            if n > maxn
                maxn = n;
                j = i;
            end
            
        end
        
    end
    
    set(h, 'Value', j);
    set(h, 'String', uvarlist);
    j =   find(strcmp(oldvar, varlist),1);
    if ~isempty(j)&&j<=length(replayopt)
        fullud = replayopt{j(1)};
    end
    
end


function updateall(hdlg, handles, replayopt,avar, at)
hslider=findobj(hdlg,'tag','TheSlider');
if at<0
    at=0;
end

if at>1
    at=1;
end

set(hslider, 'value', real(at));
tt = [];
for j = 1:length(handles)
    if ishandle(handles(j))&&~isempty(replayopt{j}.callback)
        t1 = replayopt{j}.callback(handles(j),avar, at);
        if ~isempty(t1)
            tt = t1;
        end
        
    end
    
end

h=findobj(hdlg,'tag','TimeEdit');
set(h, 'string', num2str(tt));
%drawnow('update');
drawnow;

function closeconnections(hObject, runonend)
if ishandle(hObject)
    hdlg = get(hObject, 'parent');
    ud = get(hdlg, 'userdata');
    ud.active = 0;
    set(hdlg, 'userdata', ud);
    [handles, replayopt, fullud] = makeconnections(hdlg,0);
    if runonend > 0
        for j = 1:length(handles)
            if ~isempty(replayopt{j}.onend)
                replayopt{j}.onend(handles(j), runonend == 1);
            end
            
        end
        
    end
    
    hslider=findobj(hdlg,'tag','TheSlider');
    at = get(hslider, 'value');
    h=findobj(hdlg,'tag','CurrentText');
    if ~isempty(fullud)
        setstatustext(h, 'Inactive', at, fullud);
    end
    
    h=findobj(hdlg,'Tag','playButton');
    s = get(h, 'string');
    if strcmp(s, '[Stop]')
        if isfield(ud, 'playbutton')
            set(h, 'string', ud.playbutton);
        else
            set(h,'string','Play >');
        end
        
    end
    
end


function haxes = getallaxes
haxes=findobj(get(0,'children'),'type','axes');
if ~isempty(haxes)
    tags = get(haxes, 'tag');
    haxes=haxes( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
end

