function i_seteventdlg()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
%
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.


h1 = figure( ...
   'Units','characters',...
   'PaperUnits',get(0,'defaultfigurePaperUnits'),...
   'Color', [0.941176470588235 0.941176470588235 0.941176470588235], ...
   'Colormap', [0 0 0.5625; 0 0 0.625; 0 0 0.6875; 0 0 0.75; 0 0 0.8125; 0 0 0.875; 0 0 0.9375; 0 0 1; 0 0.0625 1; 0 0.125 1; 0 0.1875 1; 0 0.25 1; 0 0.3125 1; 0 0.375 1; 0 0.4375 1; 0 0.5 1; 0 0.5625 1; 0 0.625 1; 0 0.6875 1; 0 0.75 1; 0 0.8125 1; 0 0.875 1; 0 0.9375 1; 0 1 1; 0.0625 1 1; 0.125 1 0.9375; 0.1875 1 0.875; 0.25 1 0.8125; 0.3125 1 0.75; 0.375 1 0.6875; 0.4375 1 0.625; 0.5 1 0.5625; 0.5625 1 0.5; 0.625 1 0.4375; 0.6875 1 0.375; 0.75 1 0.3125; 0.8125 1 0.25; 0.875 1 0.1875; 0.9375 1 0.125; 1 1 0.0625; 1 1 0; 1 0.9375 0; 1 0.875 0; 1 0.8125 0; 1 0.75 0; 1 0.6875 0; 1 0.625 0; 1 0.5625 0; 1 0.5 0; 1 0.4375 0; 1 0.375 0; 1 0.3125 0; 1 0.25 0; 1 0.1875 0; 1 0.125 0; 1 0.0625 0; 1 0 0; 0.9375 0 0; 0.875 0 0; 0.8125 0 0; 0.75 0 0; 0.6875 0 0; 0.625 0 0; 0.5625 0 0], ...
   'IntegerHandle','off',...
   'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
   'MenuBar','none',...
   'Name','setevent',...
   'NumberTitle','off',...
   'PaperPosition',get(0,'defaultfigurePaperPosition'),...
   'PaperSize',get(0,'defaultfigurePaperSize'),...
   'PaperType',get(0,'defaultfigurePaperType'),...
   'Position', [103.8 32.77 109.4 28.69], ...
   'Resize','off',...
   'HandleVisibility','callback',...
   'UserData', [], ...
   'Tag','eventdlg',...
   'Visible','on',...
   'CreateFcn',@(h,evnt)movegui(h, 'center'));

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @(hObject, eventdata)addline(hObject, eventdata), ...
   'FontSize', 10.6666666666667, ...
   'Position', [10 2 15 2.5], ...
   'String','Add',...
   'Tag','pushbutton2' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @(hObject, eventdata)deleteline(hObject, eventdata), ...
   'FontSize', 10.6666666666667, ...
   'Position', [30 2 15 2.5], ...
   'String','Delete',...
   'Tag','DeleteBtn',...
   'CreateFcn', @(hObject, eventdata)createbtn(hObject, eventdata));

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @(hObject, eventdata)clearall(hObject, eventdata), ...
   'FontSize', 10.6666666666667, ...
   'Position', [50 2 15 2.5], ...
   'String','Clear all',...
   'Tag','ClearallBtn' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @(hObject, eventdata)closedlg(hObject, eventdata), ...
   'FontSize', 10.6666666666667, ...
   'Position', [70 2 15 2.5], ...
   'String','OK',...
   'Tag','pushbutton6' );

uitable( ...
   'Parent', h1, ...
   'Units','characters',...
   'BackgroundColor', [1 1 1; 0.96078431372549 0.96078431372549 0.96078431372549], ...
   'CellEditCallback', @(hObject, eventdata)updategrind(hObject, eventdata), ...
   'CellSelectionCallback', @(hObject, eventdata)selectcell(hObject, eventdata), ...
   'ColumnFormat',{  [] 'logical' 'numeric' [] [] [] },...
   'ColumnEditable', true, ...
   'ColumnName',{  'Event name'; 'Enabled'; 'First time'; ' '; ' '; ' ' },...
   'ColumnWidth', {  100 60 60 150 150 150 }, ...
   'Position', [ - 0.2 5.92307692307692 110.2 23.1538461538462], ...
   'UserData', [], ...
   'Tag','eventtable',...
   'CreateFcn', @(hObject, eventdata)updatedlg(hObject, eventdata));


function closedlg(~, ~)
global g_grind;
%remove empty events
for i = length(g_grind.event.events):-1:1
   if isempty(g_grind.event.events(i).event)||(strcmp(g_grind.event.events(i).event, 'simpleevent')...
           &&(isempty(g_grind.event.events(i).arg)||~ischar(g_grind.event.events(i).arg{1})||isempty(strtrim(g_grind.event.events(i).arg{1}))))
      g_grind.event.events(i) = [];
   end
end
g_grind.event.t = 0;
g_grind.event.queue = g_grind.event.events;
i_setevent('sortqueue');
g_grind.checks.lastsettings = [];
g_grind.solver.hasevents = ~isempty(g_grind.event.queue);
close;

function updategrind(hObject, ~)
global g_grind;
hfig=getparentfig(hObject);
htable=findobj(hfig,'tag','eventtable');
dat = get(htable, 'data');
g_grind.event.events = [];
for i = 1:size(dat, 1)
   g_grind.event.events(i).event = dat{i, 1};
   g_grind.event.events(i).enabled = dat{i, 2};
   if isnumeric(dat{i, 3})
      g_grind.event.events(i).t = dat{i, 3};
      g_grind.event.events(i).strt = num2str(dat{i, 3});
   else
      g_grind.event.events(i).t = str2num(dat{i, 3}); 
      g_grind.event.events(i).strt = dat{i, 3};
   end

   args = dat(i, 4:end);
   if strcmp(dat{i, 1}, 'simpleevent')&&isempty(dat(i, 4))
      args{1} = ' ';
   end

   if strcmp(dat{i, 1}, 'simpleevent')&&isempty(dat(i, 5))
      args{2} = NaN;
   end

   ndx = cellfun(@isempty, args);
   g_grind.event.events(i).arg = args(~ndx);
end


function clearall(hObject, eventdata)
hfig=getparentfig(hObject);
htable=findobj(hfig,'tag','eventtable');
button = questdlg('Are you sure that you want to delete all events?', ...
   'Events', 'Yes', 'No', 'No');
if strcmp(button, 'Yes')
   set(htable, 'data', []);
   updategrind(hObject, eventdata)
   h=findobj(hfig,'Tag','ClearallBtn');
   set(h,'enable','off');
   h=findobj(hfig,'tag','DeleteBtn');
   set(h,'enable','off');
end


    
function deleteline(hObject, eventdata)
hfig=getparentfig(hObject);
htable=findobj(hfig,'tag','eventtable');
ud = get(htable, 'userdata');
if isfield(ud, 'Indices')&&~isempty(ud.Indices)
   button = questdlg(sprintf('Are you sure that you want to delete line %d?',ud.Indices(1)), ...
      'Events', 'Yes', 'No', 'No');
   if strcmp(button, 'Yes')
      dat = get(htable, 'data');
      dat(ud.Indices(1), :) = [];
      set(htable, 'data', dat);
      set(htable, 'userdata', []);
      drawnow;
      h=findobj(hfig,'tag','DeleteBtn');
      set(h,'enable','off');
      if isempty(dat)
         h=findobj(hfig,'Tag','ClearallBtn');
         set(h,'enable','off');
      end

      updategrind(hObject, eventdata)
   end

end


function selectcell(hObject, eventdata)
hfig=getparentfig(hObject);
htable=findobj(hfig,'tag','eventtable');
ud = get(htable, 'userdata');
ud.Indices = eventdata.Indices;
set(htable, 'userdata', ud);
h=findobj(hfig,'tag','DeleteBtn');
set(h,'enable','on');


function addline(hObject, eventdata)
hfig=getparentfig(hObject);
htable=findobj(hfig,'tag','eventtable');
dat = get(htable, 'data');
dat{end + 1, 1} = 'simpleevent';
dat{end, 2} = true;
dat{end, 3} = 0;
dat{end, 4} = ' ';
dat{end, 5} = NaN;
set(htable, 'data', dat);
updategrind(hObject, eventdata)
h=findobj(hfig,'Tag','ClearallBtn');
set(h,'enable','on');

function createbtn(hObject, ~)
set(hObject,'enable','off');

function updatedlg(hObject, ~)
global g_grind;
hfig=getparentfig(hObject);
htable=findobj(hfig,'tag','eventtable');
ncol = 6;
onlysimple = true;
dat={};
if isfield(g_grind, 'event')
   dat = cell(length(g_grind.event.events), ncol);
   for i = 1:length(g_grind.event.events)
      dat{i, 1} = g_grind.event.events(i).event;
      onlysimple = onlysimple&&strcmp(dat{i, 1}, 'simpleevent');
      dat{i, 2} = logical(g_grind.event.events(i).enabled);
      dat{i, 3} = g_grind.event.events(i).t;
      for j = 1:length(g_grind.event.events(i).arg)
         dat{i, 3 + j} =  g_grind.event.events(i).arg{j};
      end

   end
   set(htable, 'data', dat);
end

if isempty(dat)
    h=findobj(hfig,'Tag','ClearallBtn');
    set(h,'enable','off');
end

tab = get(htable, 'columnname');
if onlysimple
   tab{4} = 'Commands';
   tab{5} = 'Time interval';
else
   tab{4} = ' ';
   tab{5} = ' ';
end

set(htable, 'columnname', tab)


