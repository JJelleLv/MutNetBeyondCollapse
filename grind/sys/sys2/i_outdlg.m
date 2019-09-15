function fig = i_outdlg(flag2)
global g_grind;
if nargin <1
    if ~isempty(g_grind)&&isfield(g_grind,'outdlg')
        flag2=g_grind.outdlg;
    else
        flag2=1;
    end

end

if ~isempty(g_grind)
    N=length(g_grind.timevars);
    while (N>0)&&isempty(g_grind.timevars{N})
        N=N-1;
    end

end

if N<1
    N=1;
end

if flag2>N
    N=flag2;
end

plotlst=cell(1,N);
for i = 1:N
    plotlst{i} = sprintf('Time Plot %d', i);
end

hfig=findobj(0,'tag','outdlg');
if ~isempty(hfig)
    delete(hfig);
end

hfig = figure('Units','points', ...
    'Color',[0.914 0.914 0.914], ...
    'MenuBar','none', ...
    'Name','Edit Time Plots', ...
    'NumberTitle','off', ...
    'PaperPosition',[18 180 576 432], ...
    'PaperType','A4', ...
    'PaperUnits','points', ...
    'Position',[296.25 270.75 366 234], ...
    'Tag','outdlg', ...
    'ToolBar','none',...
    'ResizeFcn',@outresize,...
    'CreateFcn',@(h,evnt)movegui(h, 'center'));
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17.25 160 303.75 15], ...
    'String','X-axis: one variable/function', ...
    'Style','text', ...
    'Tag','StaticText1');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17.25 110 303.75 15], ...
    'String','Y-axis: list of variables/functions, delimited with spaces', ...
    'Style','text', ...
    'Tag','StaticText2');

uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17.25 197.25 72.75 15], ...
    'String','Time Plot No.', ...
    'Style','text', ...
    'Tag','StaticText3');
hplotno = uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'Callback',@c_changeplot, ...
    'ListboxTop',0, ...
    'Position',[81 195 113.25 22.5], ...
    'String',plotlst, ...
    'Style','popupmenu', ...
    'Tag','PlotNo', ...
    'TooltipString','Select Time Plot here',...
    'Value',flag2);
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_newplot, ...
    'ListboxTop',0, ...
    'Position',[243 196.5 93.75 22.5], ...
    'String','Add Time Plot', ...
    'TooltipString','Add new plot to the list of time plots',...
    'Tag','AddBtn');

uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17.25 145 318 16.5], ...
    'String','t', ...
    'Style','edit', ...
    'TooltipString','Edit the variable for the time axis here',...
    'Tag','EditTime');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Max',1, ...
    'Position',[17.25 95 318.75 16.5], ...
    'Style','edit', ...
    'TooltipString','Enter here a list of variables or functions for the Y-axis',...
    'Tag','EditList');


uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_ok, ...
    'ListboxTop',0, ...
    'Position',[31.5 20 93.75 22.5], ...
    'String','OK', ...
    'TooltipString','Close and save',...
    'Tag','OKBtn');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_cancel, ...
    'ListboxTop',0, ...
    'Position',[141 20 93.75 22.5], ...
    'String','Cancel', ...
    'TooltipString','Close and undo changes',...
    'Tag','CancelBtn');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback','commands outdlg', ...
    'ListboxTop',0, ...
    'Position',[250.5 20 93.75 22.5], ...
    'String','Help', ...
    'TooltipString','Open help window',...
    'Tag','HelpBtn');


if N==1
    set(hplotno,'Enable','off');
end

%    global g_grind;
if ~isempty(g_grind)
    ud.oudlist=g_grind.timevars;
    ud.oudoutt=g_grind.outt;
else
    ud.oudlist={};
    ud.oudoutt={};
end

ud.curr=flag2;
set(hfig,'userdata',ud);
setplot(ud.curr,hfig);
if nargout > 0, fig = hfig; end
function c_ok(hobject,~)
global g_grind;
hfig=getparentfig(hobject);
ud=get(hfig,'userdata');
try
    savecurrent(hfig,ud.curr);
catch err
    msgbox(err.message)
    rethrow(err)
end
g_grind.outdlg=ud.curr;
delete(hfig);
disp(' ');
out('-?');
function c_cancel(hobject,~)
hfig=getparentfig(hobject);
ud=get(hfig,'userdata');
global g_grind;
g_grind.timevars=ud.oudlist;
g_grind.outt=ud.oudoutt;
g_grind.outdlg=ud.curr;
delete(hfig);
disp('Cancelled');
function c_newplot(hobject,eventdata)
hfig=getparentfig(hobject);
h=findobj(hfig,'Tag','PlotNo');
plotlst = get(h,'string');
No=length(plotlst)+1;
plotlst{No} = sprintf('Time Plot %d', No);
if No>1
    set(h,'Enable','on');
end

set(h,'string',plotlst);
set(h,'value',No);
c_changeplot(hobject,eventdata);
function c_changeplot(hobject,~)
hfig=getparentfig(hobject);
ud=get(hfig,'userdata');
savecurrent(hfig,ud.curr);
h=findobj(hfig,'tag','PlotNo');
v=get(h,'Value');
setplot(v,hfig);
ud.curr=v;
set(hfig,'userdata',ud);

function outresize(hobj,~)
hfig=getparentfig(hobj);
data={'StaticText1',[1 0 1 0],[17.25 160 303.75 15];...
'StaticText2',[1 0 1 0],[17.25 110 303.75 15];...
'StaticText3',[1 0 1 0],[17.25 197.25 72.75 15];...
'PlotNo',[1 1 1 0],[81 195 113.25 22.5];...
'AddBtn',[0 1 1 0],[243 196.5 93.75 22.5];...
'EditTime',[1 1 1 0],[17.25 145 318 16.5];...
'EditList',[1 1 1 0],[17.25 95 318.75 16.5];...
'OKBtn',[1 0  0 1],[31.5 20 93.75 22.5];...
'CancelBtn',[1 0  0 1],[141 20 93.75 22.5];...
'HelpBtn',[1 0  0 1],[250.5 20 93.75 22.5]};
i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[296.25 162 366 234]);
i_resizedlg(hfig,'horizcenter',{'OKBtn','CancelBtn','HelpBtn'});

function savecurrent(hfig,No)
h=findobj(hfig,'Tag','EditList');
lst=get(h,'String');
if ~isempty(lst)
    eval(sprintf('out -silent -%d %s', No, lst));
else
    eval(sprintf('out -silent -%d -clear',No));
end

h=findobj(hfig,'tag','EditTime');
tim=strtrim(get(h,'string'));
if isempty(tim)
    tim='t';
end

global g_grind
g_grind.outt{No}=tim;

function setplot(No,f)
global g_grind;
if isempty(g_grind)||(No>length(g_grind.timevars))
    s='';
    tim='t';
else
    vars = g_grind.timevars{No};
    for i = 1:length(vars)
        if ~isempty(strfind(vars{i}, ' '))
            vars{i}=['''' vars{i} ''''];
        end

    end

    s = sprintf('%s ', vars{:});
    tim=g_grind.outt{No};
end

h=findobj(f,'tag','EditList');
set(h,'string',s);
h=findobj(f,'tag','EditTime');
set(h,'string',tim);
