function fig = i_viewcellsdlg
global g_grind;
if ~isempty(g_grind) && ~isfield(g_grind, 'viewcells')
    viewcells('-defaults');
end

if ~isempty(g_grind)
    N = length(g_grind.viewcells.vars);
    while (N > 0) && isempty(g_grind.viewcells.vars{N})
        N = N - 1;
    end

end

if N < 1
    N = 1;
end

plotlst = cell(N);
for i = 1:N
    plotlst{i} = sprintf('Plot No. %d', i);
end

%delete double dialogs
ch=get(0,'children');
ifig=find(strcmp(get(ch,'tag'),'viewcellsdlg'));
if ~isempty(ifig)
    close(ch(ifig));
end

hfig = figure('Units','points', ...
    'Color',[0.914 0.914 0.914], ...
    'MenuBar','none', ...
    'Name','Plots for Viewcells', ...
    'NumberTitle','off', ...
    'PaperPosition',[18 180 576 432], ...
    'PaperType','A4', ...
    'PaperUnits','points', ...
    'Position',[199 232 366 205], ...
    'Tag','viewcellsdlg', ...
    'ToolBar','none',...
    'CreateFcn',@(h,evnt)movegui(h, 'center'));
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17 180 - 40 303.75 15], ...
    'String','Auxliiary or state variable', ...
    'Style','text', ...
    'Tag','StaticText2');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17.25 138.75 - 40 303.75 15], ...
    'String','Size of variable', ...
    'Style','text', ...
    'Tag','StaticText2');

uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[17 206 - 40 72.75 15], ...
    'String','Viewcells plot No.', ...
    'Style','text', ...
    'Tag','StaticText1');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'Callback',@c_changeplot, ...
    'ListboxTop',0, ...
    'Position',[105 205.5 - 40 113.25 22.5], ...
    'String',plotlst, ...
    'Style','popupmenu', ...
    'Tag','PlotNo', ...
    'TooltipString','Select Plot number here', ...
    'Value', 1);
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_newplot, ...
    'ListboxTop',0, ...
    'Position',[243 205.5 - 40 93.75 22.5], ...
    'String','Add Viewcells Plot', ...
    'Tag','AddBtn', ...
    'TooltipString','Add new plot to the list of viewcells');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'HorizontalAlignment','left', ...
    'Callback',@c_xchanged, ...
    'ListboxTop',0, ...
    'Position',[17 168 - 40 318 16.5], ...
    'Style','edit', ...
    'Tag','EditX', ...
    'TooltipString','Edit the variable for the plot here');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Max',1, ...
    'Position',[17 126 - 40 317 17], ...
    'Style','edit', ...
    'Tag','EditSize', ...
    'TooltipString','Enter here the size of the variable');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_ok, ...
    'ListboxTop',0, ...
    'Position',[31.5 8.25 93.75 22.5], ...
    'String','OK', ...
    'Tag','OKBtn', ...
    'TooltipString','Close and save');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_cancel, ...
    'ListboxTop',0, ...
    'Position',[141 8.25 93.75 22.5], ...
    'String','Cancel', ...
    'Tag','CancelBtn', ...
    'TooltipString','Close and undo changes');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback','commands viewcells', ...
    'ListboxTop',0, ...
    'Position',[250.5 8.25 93.75 22.5], ...
    'String','Help', ...
    'Tag','HelpBtn', ...
    'TooltipString','Open help window');


%if N==1
%   set(hplotno,'Enable','off');
% end

%    global g_grind;
if ~isempty(g_grind)
    ud.oudlist = g_grind.viewcells;
else
    ud.oudlist = {};
end

ud.curr = g_grind.viewcells.currno;
set(hfig, 'userdata', ud);
setplot(hfig,ud.curr);
if nargout > 0, fig = hfig; end
function c_ok(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
savecurrent(hfig,ud.curr);
delete(hfig);
disp(' ');
viewcells('-list');
viewcells();
function c_cancel(hobject,~)
global g_grind;
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
g_grind.viewcells = ud.oudlist;
delete(hfig);
disp(' ');
viewcells('-list');
function c_newplot(hobject,~)
hfig=getparentfig(hobject);
h=findobj(hfig,'Tag','PlotNo');
plotlst = get(h, 'string');
No = length(plotlst) + 1;
plotlst{No} = sprintf('Plot No. %d', No);
if No > 1
    set(h,'Enable','on');
end

set(h, 'string', plotlst);
set(h, 'value', No);
c_changeplot(hfig);
function c_xchanged(hobject,~)
global g_grind;
hfig=getparentfig(hobject);
h=findobj(hfig,'tag','EditX');
nam = get(h, 'string');
iX = i_getno(nam);
d = [];
if iX.isvar
    d = g_grind.statevars.dims{iX.vecno};
elseif iX.isfun
    d = g_grind.funcnames.dims{iX.vecno};
end

if ~isempty(d)
    h=findobj(hfig,'tag','EditSize');
    set(h,'string',sprintf('[%d %d]',[d.dim1,d.dim2]));
end

function c_changeplot(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
savecurrent(hfig,ud.curr);
h=findobj(hfig,'tag','PlotNo');
v = get(h, 'Value');
setplot(hfig,v);
ud.curr = v;
set(hfig, 'userdata', ud);


function savecurrent(hfig,No)
global g_grind;
h=findobj(hfig,'tag','EditX');
g_grind.viewcells.vars{No}.name = get(h, 'string');
h=findobj(hfig,'tag','EditSize');
g_grind.viewcells.vars{No}.dims = str2num(get(h, 'string')); 
g_grind.viewcells.currno = No;

function setplot(hfig,No)
global g_grind;
if isempty(g_grind) || (No > length(g_grind.viewcells.vars))
    x = '';
    siz = [1 1];
else
    x = g_grind.viewcells.vars{No}.name;
    siz = g_grind.viewcells.vars{No}.dims;
end

h=findobj(hfig,'tag','EditX');
set(h, 'string', x);
h=findobj(hfig,'tag','EditSize');
set(h,'string',sprintf('[%d %d]',siz));
h=findobj(hfig,'tag', 'PlotNo');
set(h, 'value', No);
