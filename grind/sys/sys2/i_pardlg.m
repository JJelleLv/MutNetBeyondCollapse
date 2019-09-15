function i_pardlg(debugunits)
global g_grind;
if nargin==0
    debugunits=false;
end

defgroup = '-';
if ~isfield(g_grind, 'pargroups')
   g_grind.pargroups = cell(size(g_grind.pars));
   for i = 1:length(g_grind.pargroups)
      g_grind.pargroups{i} = defgroup;
   end

end


if ~g_grind.statevars.vector
   pardata = cell(length(g_grind.pars) + length(g_grind.statevars.names), 5);
   for i = 1:length(g_grind.statevars.names)
      avar = g_grind.statevars.names{i};
      pardata{i, 1} = avar;
      pardata{i, 2} = 'initial condition';
%      pardata{i, 3} = par('-d', g_grind.statevars.names{i});
      valu = evalin('base', g_grind.statevars.names{i});
      if isnumeric(valu) || islogical(valu)
         if numel(valu) > 1000
            pardata{i, 4} = sprintf('[%dx%d %s]', size(valu), class(valu));
         else
            pardata{i, 4} =  mat2str(valu);
         end
      end

      [pardata{i, 3},pardata{i, 5}] = par('-d', g_grind.statevars.names{i});
      pardata{i, 6} = valu;
      if numel(valu) == 1
         pardata{i, 7}=sprintf('%s = %g;',avar,valu);
      else
         if max(max(valu)) - min(min(valu)) < eps
            pardata{i, 7}=sprintf('%s = zeros(%d,%d)+%g;',avar,size(valu),mean(mean(valu)));
         else
            pardata{i, 7}=sprintf('%s = %s;',avar,pardata{i,4});
         end

      end

   end

   n = length(g_grind.statevars.names);
else
   pardata  = cell(length(g_grind.pars), 5);
   n = 0;
end

for i = 1:length(g_grind.pars)
   pardata{i + n, 1} = g_grind.pars{i};
   pardata{i + n, 2} = g_grind.pargroups{i};
  % pardata{i + n, 3} = par('-d', g_grind.pars{i});
   valu = evalin('base', g_grind.pars{i});
   if isnumeric(valu) || islogical(valu)
      if numel(valu) > 1000
         pardata{i + n, 4} = sprintf('[%dx%d %s]', size(valu), class(valu));
      else
         pardata{i + n, 4} =  mat2str(valu);
      end
   end

   [ pardata{i + n, 3},pardata{i + n, 5}] = par('-d', g_grind.pars{i});
   pardata{i + n, 6} = valu;
   if numel(valu) == 1
      pardata{i + n, 7}=sprintf('%s = %g;',g_grind.pars{i},valu);
   else
      if max(max(valu)) - min(min(valu)) < eps
         pardata{i+n, 7}=sprintf('%s = zeros(%d,%d)+%g;',g_grind.pars{i},size(valu),mean(mean(valu)));
      else
         pardata{i + n, 7}=sprintf('%s = %s;',g_grind.pars{i},pardata{i,4});
      end

   end

end

n=n+length(g_grind.pars);
for i=1:length(g_grind.externvars)
   pardata{i + n, 1} = g_grind.externvars{i}.name;
   pardata{i + n, 2} = 'external vars';
   pardata{i + n, 3} = par('-d', g_grind.externvars{i}.name);
   pardata{i + n, 4} = g_grind.externvars{i}.default;
   pardata{i + n, 5} = par('-u', g_grind.externvars{i}.name);
   pardata{i + n, 6} = g_grind.externvars{i}.default;
   pardata{i + n, 7} = sprintf('%s = %g;',g_grind.externvars{i}.name,g_grind.externvars{i}.default);
end

parfilter={'all','all','all','','all'};
ud.Indices = [1, 1];
ud.debugunits=debugunits;
figs=findall(0,'type','figure');
tagndx=strcmp(get(figs,'tag'),'par.edit');
%make sure that not more than one dialog boxes are created 
if any(tagndx)
      hfig = figs(find(tagndx, 1));
      figure(hfig);
      set(hfig,'UserData', ud)
      h2=findobj(hfig,'Tag','FilterTable');
      set(h2,'Data',parfilter);
      h2=findobj(hfig,'Tag','ParTable');
      set(h2,'UserData', pardata);
      set(h2, 'Data', pardata(:,1:5));
      movegui(hfig,'onscreen');
else
hfig = figure( ...
   'Units','characters',...
   'PaperUnits',get(0,'defaultfigurePaperUnits'),...
   'Color', [0.94 0.94 0.94], ...
   'IntegerHandle','off',...
   'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
   'MenuBar','none',...
   'Name','Par. Edit',...
   'NumberTitle','off',...
   'PaperPosition',get(0,'defaultfigurePaperPosition'),...
   'PaperSize',get(0,'defaultfigurePaperSize'),...
   'PaperType',get(0,'defaultfigurePaperType'),...
   'Position', [103.8 29 121.4 32.4615384615385], ...
   'Resize','on',...
   'HandleVisibility','callback',...
   'UserData', ud, ...
   'Tag','par.edit',...
   'Visible','on' ,...
   'CreateFcn',@(h,evnt)movegui(h, 'center'));



uicontrol( ...
   'Parent', hfig, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @checkunitsbuttonclick, ...
   'FontSize', 10.6666666666667, ...
   'Position', [21.2 0.846 14 1.69230769230769], ...
   'String','Analyse units',...
   'TooltipString','Analyse and add units',...
   'Tag','CheckunitsBtn' );
uicontrol( ...
   'Parent', hfig, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @editparbuttonclick, ...
   'FontSize', 10.6666666666667, ...
   'Position', [39.8 0.846 13.8 1.69230769230769], ...
   'String','Edit par',...
   'TooltipString','Edit each parameter separately',...
   'Tag','EditparBtn' );
uicontrol( ...
   'Parent', hfig, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @okbuttonclick, ...
   'FontSize', 10.6666666666667, ...
   'Position', [58.4 0.846153846153846 14 1.69230769230769], ...
   'String','OK',...
   'TooltipString','OK',...
   'Tag','OKBtn' );

uicontrol( ...
   'Parent', hfig, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @cancelbuttonclick, ...
   'FontSize', 10.6666666666667, ...
   'Position', [77 0.846153846153846 14 1.69230769230769], ...
   'String',' Cancel',...
   'TooltipString','Cancel',...
   'Tag','CancelBtn' );

uicontrol( ...
   'Parent', hfig, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @helpbuttonclick, ...
   'FontSize', 10.6666666666667, ...
   'Position', [96 0.846153846153846 14 1.69230769230769], ...
   'String','Help',...
   'TooltipString','Help',...
   'Tag','HelpBtn' );
h2=uitable( ...
   'Parent', hfig, ...
   'Units','characters',...
   'BackgroundColor', [0.96 0.96 0.96], ...
   'ColumnFormat', [{getparfilter(pardata(:, 1))}, {getparfilter(pardata(:, 2))}, ...
   {getparfilter(pardata(:, 3))}, 'char', {getparfilter(pardata(:, 5))}], ...
   'ColumnEditable', true(1, 5), ...
   'ColumnName',{  'Parameter'; 'Group'; 'Description'; 'Value'; 'Unit' },...
   'ColumnWidth',  {  60 80 280 80 80 }, ...
   'CellEditCallback', @filtertableedit, ...
   'TooltipString','Select a filter here',...
   'Data', parfilter, ...
   'Position', [0 29.4615384615385 120.2 3.23076923076923], ...
   'RowName', [], ...
   'FontWeight','light',...
   'UserData', [], ...
   'Tag','FilterTable' );
set(h2,'units','normalized');
htable =  uitable('Parent', hfig, ...
   'Units','characters',...
   'BackgroundColor', [1 1 1; 0.96 0.96 0.96], ...
   'ColumnFormat', {  [] [] [] [] [] }, ...
   'ColumnEditable', [false true true true true], ...
   'KeyPressFcn', @onKeyPressData, ...
   'CellSelectionCallback', @ParTableCellSelection, ...
   'CellEditCallback', @partableedit, ...
   'ColumnName', [], ...
   'TooltipString','Press Alt-Enter to change parameter',...
   'ColumnWidth', {  60 80 280 80 80 }, ...
   'Data', pardata(:,1:5), ...
   'Position', [0 3 120.2 26.5384615384615], ...
   'RowName', [], ...
   'UserData', pardata, ...
   'Tag','ParTable' );
set(htable,'units','normalized');

%    uicontrol( ...
%       'Parent', hfig, ...
%       'Units','characters',...
%       'FontUnits','pixels',...
%       'Callback',@copybuttonclick,...
%       'FontSize', 10.6666666666667, ...
%       'Position', [39.8 0.846 13.8 1.69230769230769], ...
%       'String','Clipboard',...
%       'Tag','CopyBtn' );

hcmenu = uicontextmenu('Parent', hfig);
% Define the context menu items and install their callbacks
uimenu(hcmenu, 'Label', 'Copy Table', 'Callback',@copybuttonclick);
uimenu(hcmenu, 'Label', 'Clear units', 'Callback',@clearunitsclick);

%uimenu(hcmenu, 'Label', 'Paste', 'Callback',@copybuttonclick);
set(htable, 'uicontextmenu', hcmenu);
end

%  guide(hfig);
%   uiwait(hfig);
%    if ishandle(hfig)
%       htab=findobj(hfig,'Tag','ParTable');
%       allpardata = get(htab, 'userdata');
%       for i=1:size(allpardata,1)
%           apar=allpardata{i,1};
%           assignin('base',apar,allpardata{i, 6});
%           par('-set',apar, sprintf('%s,%s',allpardata{i, 3},allpardata{i, 5}));
%           g_grind.pargroups{i}=allpardata{i,2};
%       end

%       close(hfig);
%    end

% function resizetable(src,evnt)
% htab=findobj(hfig,'Tag','FilterTable');
% set(htab,'Units','characters'); 
% set(htab, 'Position', [0 3 120.2 26.5384615384615]);
% set(htab,'Units','Normalized');
% upos=get(htab, 'position');
% set(htab, 'position',[0,1-upos(4),1,upos(4)]);
% htab=findobj(hfig,'Tag','ParTable');
% set(htab,'Units','normalized'); 
% %left bottom width height
% set(htab,'Position', [0 upos(4) 1 0.8-upos(4)]);

function editparbuttonclick(src, ~)
hfig = getparentfig(src);
htab=findobj(hfig,'Tag','ParTable');
pardata = get(htab, 'data');
allpardata = get(htab, 'userdata');
ud = get(hfig, 'userdata');
if isstruct(ud) && isfield(ud, 'Indices') && ~isempty(ud.Indices)
   i = ud.Indices(1);
else
   i = 1;
end

j = strcmp(allpardata(:,1), pardata{i,1});
res = i_editpardlg(allpardata, find(j));
if ~isempty(res)
   allpardata = res;
end

set(htab, 'data', allpardata(:, 1:5));
set(htab, 'userdata', allpardata);
UpdateParfilter(hfig, allpardata)

function checkunitsbuttonclick(src, ~)
%function check the units
hfig = getparentfig(src);
htab=findobj(hfig,'Tag','ParTable');
allpardata = get(htab, 'userdata');
ud = get(hfig, 'userdata');
try
    [un, vars] = analunits({},allpardata(:, 1), allpardata(:, 5),ud.debugunits);
catch err
    h=i_errordlg(err.message,'par.edit','modal');
    set(findobj(h,'tag','MessageBox'),'fontname','Courier New', 'fontsize',8)
 %   set(h,'fontname','New courier');
    rethrow(err)
end

nundef=0;
ndef=0;
for i = 1:size(allpardata, 1)
   ndx = strcmp(allpardata{i, 1}, vars);
   if any(ndx)
      if un(ndx).isundefined
         allpardata{i, 5}  = '';
         nundef=nundef+1;
      else 
         if isempty(allpardata{i, 5})
             ndef=ndef+1;
         end

         allpardata{i, 5} = un(ndx).unit;
      end

   end

end
if ndef==0&&nundef>0
   msgbox(sprintf('No problems found, no units added but insufficient information to fill %d units.\nEnter some more units.',nundef),'analunits'); 
elseif ndef>0&&nundef>0    
   msgbox(sprintf('No problems found, %d units added but insufficient information to fill %d units.\nEnter some more units.',ndef,nundef),'analunits'); 
elseif ndef>0&&nundef==0    
   msgbox(sprintf('No problems found, %d units added, all units are now defined',ndef),'analunits'); 
elseif ndef==0&&nundef==0    
   msgbox('No problems found, all units consistent','analunits'); 
end

set(htab, 'data', allpardata(:, 1:5));
set(htab, 'userdata', allpardata);
UpdateParfilter(hfig, allpardata)

function clearunitsclick(src, ~)
hfig = getparentfig(src);
if strcmp('Yes',questdlg('Are you sure to clear all units?', 'par edit','Yes','No','No'))
    htab=findobj(hfig,'Tag','ParTable');
     allpardata = get(htab, 'userdata');
     allpardata(:, 5)={''};
     set(htab, 'data', allpardata(:, 1:5));
     set(htab, 'userdata', allpardata);
     UpdateParfilter(hfig, allpardata)
end



function onKeyPressData(src, evnt)
if any(strcmp(evnt.Key,'return')) && any(strcmp('alt',evnt.Modifier))
   editparbuttonclick(src, evnt)
end


function okbuttonclick(src, ~)
global g_grind;
hfig = getparentfig(src);
if ishandle(hfig)
   ud=get(hfig,'userdata');
   htab=findobj(hfig,'Tag','ParTable');
   allpardata = get(htab, 'userdata');
   for i = 1:size(allpardata, 1)
      apar = allpardata{i, 1};
      f=find(strcmp(apar,g_grind.pars),1);
      assignin('base',apar,allpardata{i, 6});
      if ~par('-set',apar, sprintf('%s,%s',allpardata{i, 3},allpardata{i, 5}))
          if ~(isempty(allpardata{i,3})&&isempty(allpardata{i,5}))
              g_grind.commands{end+1}=sprintf('%s=%s;   %%%s,%s',allpardata{i,1},allpardata{i,4},allpardata{i,3},allpardata{i,5});
          end

      end

      if ~isempty(f)
         g_grind.pargroups{f} = allpardata{i, 2};
      end

   end

   close(hfig);
   if ud.debugunits
      analunits('-c','-d');
   else
      analunits('-c');
   end
end

%uiresume(hfig);

function cancelbuttonclick(src, ~)
%uiresume(hfig);
hfig = getparentfig(src);
close(hfig);

function helpbuttonclick(~, ~)
commands('par');

function copybuttonclick(src, ~)
hfig = getparentfig(src);
htab=findobj(hfig,'Tag','ParTable');
pardata = get(htab, 'data');
hfilt=findobj(hfig,'Tag','FilterTable');
colnames = get(hfilt, 'ColumnName');
varcopy([transpose(colnames); pardata]);

function ParTableCellSelection(src, evnt)
hfig = getparentfig(src);
ud=get(hfig,'userdata');
ud.Indices=evnt.Indices;
set(hfig, 'Userdata', ud);

function partableedit(src, evnt)
hfig = getparentfig(src);
pardata = get(src, 'data');
par = pardata{evnt.Indices(1), 1};
allpardata = get(src, 'userdata');
i = strcmp(allpardata(:, 1), par);
allpardata{i, evnt.Indices(2)} = evnt.NewData;
if evnt.Indices(2) == 4
   allpardata{i, 6} = str2num(evnt.NewData);  %#ok<ST2NM>
end

set(src, 'userdata',allpardata);
UpdateParfilter(hfig, allpardata, evnt.Indices(2))

function filtertableedit(src, evnt)
hfig = getparentfig(src);
if any([1 2 3 5] == evnt.Indices(2))
   htab=findobj(hfig,'Tag','ParTable');
   allpardata = get(htab, 'userdata');
   currfilter = evnt.NewData;
   if strcmp(currfilter, 'all')
      set(htab, 'data', allpardata(:,1:5))
   else
      ndx = strcmp(currfilter, allpardata(:, evnt.Indices(2)));
      set(htab, 'data', allpardata(logical(ndx), 1:5))
   end

   filter = {'all'};
   filter = filter(ones(1, 5));
   filter{evnt.Indices(2)} = currfilter;
   set(src, 'data', filter);
end


function UpdateParfilter(h, allpardata, ipar)
hfilt=findobj(h,'Tag','FilterTable');
if nargin > 2
   if any(ipar == [1 2 3 5])
      filt = get(hfilt, 'ColumnFormat');
      filt{ipar}  = getparfilter(allpardata(:, ipar));
      set(hfilt, 'ColumnFormat', filt);
   end

else
   set(hfilt, 'ColumnFormat', [{getparfilter(allpardata(:, 1))}, {getparfilter(allpardata(:, 2))}, ...
      {getparfilter(allpardata(:, 3))}, 'char', {getparfilter(allpardata(:, 5))}]);
end


function f = getparfilter(pardata)
f=[{'all'},sort(unique(transpose(pardata)))];
ndx = false(size(f));
for i = 1:length(f)
   ndx(i) = ~isempty(f{i});
end

f = f(logical(ndx));


