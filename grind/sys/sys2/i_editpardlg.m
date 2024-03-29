function res = i_editpardlg(c, ic)
if nargin == 0
   avalue = rand(20);
   c={'A', 'group','descr A', mat2str(avalue),'m-1',avalue,'A = rand(20,20)';...
      'B', 'group','descr B', '10','m-1',10,'B = 10'};
end

if nargin < 2
   ic = 1;
end

apar = c{ic, 1};
agroup = c{ic, 2};
adescr = c{ic, 3};
aunit = c{ic, 5};
p = str2num(c{ic, 4}); 
if ~isempty(p)
   avalue = p;
else
   avalue = c{ic, 6};
end

if numel(avalue) > 1000
   acode=sprintf('%s = zeros(%d,%d)+%g',apar,size(avalue),mean(mean(avalue)));
else
   acode = c{ic, 7};
end
ud.data = c;
ud.ipar = ic;
h1 = figure( ...
   'Units','characters',...
   'PaperUnits',get(0,'defaultfigurePaperUnits'),...
   'Color', [0.941176470588235 0.941176470588235 0.941176470588235], ...
   'Colormap', [0 0 0.5625; 0 0 0.625; 0 0 0.6875; 0 0 0.75; 0 0 0.8125; 0 0 0.875; 0 0 0.9375; 0 0 1; 0 0.0625 1; 0 0.125 1; 0 0.1875 1; 0 0.25 1; 0 0.3125 1; 0 0.375 1; 0 0.4375 1; 0 0.5 1; 0 0.5625 1; 0 0.625 1; 0 0.6875 1; 0 0.75 1; 0 0.8125 1; 0 0.875 1; 0 0.9375 1; 0 1 1; 0.0625 1 1; 0.125 1 0.9375; 0.1875 1 0.875; 0.25 1 0.8125; 0.3125 1 0.75; 0.375 1 0.6875; 0.4375 1 0.625; 0.5 1 0.5625; 0.5625 1 0.5; 0.625 1 0.4375; 0.6875 1 0.375; 0.75 1 0.3125; 0.8125 1 0.25; 0.875 1 0.1875; 0.9375 1 0.125; 1 1 0.0625; 1 1 0; 1 0.9375 0; 1 0.875 0; 1 0.8125 0; 1 0.75 0; 1 0.6875 0; 1 0.625 0; 1 0.5625 0; 1 0.5 0; 1 0.4375 0; 1 0.375 0; 1 0.3125 0; 1 0.25 0; 1 0.1875 0; 1 0.125 0; 1 0.0625 0; 1 0 0; 0.9375 0 0; 0.875 0 0; 0.8125 0 0; 0.75 0 0; 0.6875 0 0; 0.625 0 0; 0.5625 0 0], ...
   'IntegerHandle','off',...
   'ResizeFcn',@resize,...
   'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
   'MenuBar','none',...
   'Name','Edit parameter',...
   'NumberTitle','off',...
   'PaperPosition',get(0,'defaultfigurePaperPosition'),...
   'PaperSize',get(0,'defaultfigurePaperSize'),...
   'PaperType',get(0,'defaultfigurePaperType'),...
   'Position', [52.8 10.8461538461538 99.2 33.7692307692308], ...
   'Resize','on',...
   'UserData', ud, ...
   'Tag','i_editpardlg',...
   'CreateFcn',@(h,evnt)movegui(h, 'center'));

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @Onparpopup, ...
   'BackgroundColor', [1 1 1], ...
   'FontSize', 10.6666666666667, ...
   'Position', [32.8 31.7692307692308 40.2 1.53846153846154], ...
   'String', c(:, 1), ...
    'Tooltipstring','Select parameter',...
   'Style','popupmenu',...
   'Value', ic, ...
   'Tag','parpopup' );

uitable( ...
   'Parent', h1, ...
   'Units','characters',...
   'ColumnEditable', true(1, size(avalue,2)), ...
   'Data', avalue, ...
    'Tooltipstring','Edit value(s) of parameter',...
   'Position', [9.8 1.61538461538462 65.4 16.7692307692308], ...
   'Tag','varTable');

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'BackgroundColor', [1 1 1], ...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [32.6 27 41.2 1.69230769230769], ...
    'Tooltipstring','Edit the description of the parameter',...
   'String', {  adescr }, ...
   'Style','edit',...
   'Tag','DescrEdit');

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'BackgroundColor', [1 1 1], ...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [32.6 29.3076923076923 41 1.69230769230769], ...
    'Tooltipstring','Edit the name of the group for this parameter (if any)',...
   'String', {  agroup }, ...
   'Style','edit',...
   'Tag','GroupEdit' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'BackgroundColor', [1 1 1], ...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [32.6 24.6923076923077 41.2 1.69230769230769], ...
    'Tooltipstring','Edit the unit of the parameter',...
   'String', {  aunit }, ...
   'Style','edit',...
   'Tag','UnitEdit');

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [9.8 27.2307692307692 17.4 1.07692307692308], ...
   'String','Description:',...
   'Style','text',...
   'Tag','text1' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [9.8 29.6153846153846 10.4 1.07692307692308], ...
   'String','Group:',...
   'Style','text',...
   'Tag','text2' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [9.8 24.8461538461538 10.4 1.07692307692308], ...
   'String','Unit:',...
   'Style','text',...
   'Tag','text4');

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'BackgroundColor', [1 1 1], ...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',... 
    'Tooltipstring','Edit MATLAB code for generating the parameter value(s), e.g. rand(10,1)',...
   'Max', 10, ...
   'Position', [9.8 19.2307692307692 65.2 2.92307692307692], ...
   'String', {  acode }, ...
   'Style','edit',...
   'Tag','CodeEdit');

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [9.8 22.2307692307692 60.6 1.07692307692308], ...
   'String','MATLAB code to set the parameter',...
   'Style','text',...
   'Tag','text5');

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback','commands(''par'')',...
   'FontSize', 10.6666666666667, ...
   'Position', [79.8 1.53846153846154 13.8 1.69230769230769], ...
   'String','Help',...
   'Tooltipstring','Open help system for par',...
   'Tag','Helpbutton' );

h = uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'Callback', @OnUpdateCClik, ...
   'FontUnits','pixels',...
   'FontSize', 10.6666666666667, ...
   'Position', [79.8 16.0769230769231 13.8 1.69230769230769], ...
   'String','Update code',...
    'Tooltipstring','Update the code based on the entered values',...
   'Tag','UpdateCbutton' );
if numel(avalue) > 1000
   set(h,'Visible','off');
end

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @OnUpdateTClik, ...
   'FontSize', 10.6666666666667, ...
   'Position', [79.8 19.8461538461538 13.8 1.69230769230769], ...
    'Tooltipstring','Update the value(s) of the parameter using the MATLAB code',...
   'String','Update table',...
   'Tag','UpdateTbutton' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @OnCancelClik, ...
   'FontSize', 10.6666666666667, ...
   'Position', [79.8 4.53846153846154 13.8 1.69230769230769], ...
    'Tooltipstring','Close without saving changes',...
   'String','Cancel',...
   'Tag','Cancelbutton' );
uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'Callback', @OnOKClik, ...
   'FontSize', 10.6666666666667, ...
   'Position', [79.8 7.53846153846154 13.8 1.69230769230769], ...
    'Tooltipstring','Save changes and exit',...
   'String','OK',...
   'Tag','OKbutton' );

uicontrol( ...
   'Parent', h1, ...
   'Units','characters',...
   'FontUnits','pixels',...
   'FontSize', 10.6666666666667, ...
   'HorizontalAlignment','left',...
   'Position', [9.8 31.7692307692308 20 1.07692307692308], ...
   'String','Parameter:',...
   'Style','text',...
   'Tag','Partext');

res = {};
%guide(h1);
uiwait(h1);
if ishandle(h1)
   ud = get(h1, 'userdata');
   ud = saveCurrentPar(h1,ud);
   res = ud.data;
   close(h1);
end

%
function res = CurrentPar(hfig)
h=findobj(hfig,'Tag','parpopup');
i = get(h, 'value');
pars = get(h, 'string');
res = pars{i};

function Onparpopup(hObject, ~)
hfig=getparentfig(hObject);
ud = get(hfig, 'userdata');
ud = saveCurrentPar(hfig,ud);
ipar = get(hObject, 'value');
ud.ipar = ipar;
set(hfig, 'userdata', ud);
h=findobj(hfig,'Tag','GroupEdit');
set(h, 'string', ud.data(ud.ipar,2));
h=findobj(hfig,'Tag','DescrEdit');
set(h, 'string',ud.data(ud.ipar,3));
h=findobj(hfig,'Tag','varTable');
set(h, 'data', ud.data{ud.ipar,6});
if numel(ud.data{ud.ipar, 6}) > 1000
   ud.data{ud.ipar, 4}=sprintf('%s = zeros(%d,%d)+%g',ud.data(ud.ipar,2),...
       size(ud.data{ud.ipar, 6}),mean(mean(ud.data{ud.ipar, 6})));
else
   ud.data{ud.ipar, 4} =  mat2str(ud.data{ud.ipar, 6});
end
h=findobj(hfig,'Tag','UnitEdit');
set(h, 'string',ud.data(ud.ipar,5));
h=findobj(hfig,'Tag','CodeEdit');
set(h, 'string',ud.data(ud.ipar,7));


function ud1 = saveCurrentPar(h1,ud)
h=findobj(h1,'Tag','GroupEdit');
ud.data{ud.ipar,2} = get(h, 'string');
h=findobj(h1,'Tag','DescrEdit');
ud.data{ud.ipar,3} = get(h, 'string');
h=findobj(h1,'Tag','varTable');
ud.data{ud.ipar,6} = get(h, 'data');
if numel(ud.data{ud.ipar, 6}) > 1000
   ud.data{ud.ipar, 4} = sprintf('[%dx%d %s]', size(ud.data{ud.ipar, 6}), class(ud.data{ud.ipar, 6}));
else
   ud.data{ud.ipar, 4} =  mat2str(ud.data{ud.ipar, 6});
end
h=findobj(h1,'Tag','UnitEdit');
ud.data{ud.ipar,5} = get(h, 'string');
h=findobj(h1,'Tag','CodeEdit');
ud.data{ud.ipar,7} = get(h, 'string');
for i = 1:7
   if iscell(ud.data{ud.ipar, i})
      ud.data{ud.ipar,i} = strtrim(sprintf('%s\n', ud.data{ud.ipar,i}{:}));
   end

end

ud1 = ud;

function OnUpdateTClik(hObject, ~)
hfig=getparentfig(hObject);
h=findobj(hfig,'Tag','CodeEdit');
g_code = get(h, 'string');
f=strfind(g_code, '=');
if length(f) == 1
   if iscell(g_code)
      g_code = g_code{1}(f{1} + 1:end);
   else
      g_code = g_code(f + 1:end);
   end

   l__v = eval(g_code);
else
   clear h f;
   eval(sprintf('%s;\n', g_code{:}));
   l__v = eval(CurrentPar(hfig));
end

h=findobj(hfig,'Tag','varTable');
set(h, 'data', l__v);
%

function OnUpdateCClik(hObject, ~)
hfig=getparentfig(hObject);
apar = CurrentPar(hfig);
h=findobj(hfig,'Tag','varTable');
valu = get(h, 'data');
if numel(valu) == 1
   acode=sprintf('%s = %g;',apar, valu);
else
   if max(max(valu)) - min(min(valu)) < eps
      acode=sprintf('%s = zeros(%d,%d)+%g;',apar,size(valu),mean(mean(valu)));
   elseif numel(valu) < 1000
      acode=sprintf('%s = %s;', apar, mat2str(valu));
   else
      i_errordlg('Variable too large to write as string');
   end

end

h=findobj(hfig,'Tag','CodeEdit');
set(h, 'string', {acode});
%
function OnCancelClik(hObject, ~)
hfig=getparentfig(hObject);
close(hfig);
%
function OnOKClik(hObject, ~)
hfig=getparentfig(hObject);
uiresume(hfig);

function resize(hobj,~)
hfig=getparentfig(hobj);
data={...
'parpopup',[1 1 1 0],[123 309.75 150.75 15];...
'varTable',[1 1 1 1],[36.75 15.75 245.25 163.5];...
'DescrEdit',[1 1 1 0],[122.25 263.25 154.5 16.5];...
'GroupEdit',[1 1 1 0],[122.25 285.75 153.75 16.5];...
'UnitEdit',[1 1 1 0],[122.25 240.75 154.5 16.5];...
'text1',[1 0 1 0],[36.75 265.5 65.25 10.5];...
'text2',[1 0 1 0],[36.75 288.75 39 10.5];...
'text4',[1 0 1 0],[36.75 242.25 39 10.5];...
'CodeEdit',[1 1 1 0],[36.75 187.5 244.5 28.5];...
'text5',[1 0 1 0],[36.75 216.75 227.25 10.5];...
'Helpbutton',[0 1 0 1],[299.25 15 51.75 16.5];...
'UpdateCbutton',[0 1 0 1],[299.25 156.75 51.75 16.5];...
'UpdateTbutton',[0 1 0 1],[299.25 193.5 51.75 16.5];...
'Cancelbutton',[0 1 0 1],[299.25 44.25 51.75 16.5];...
'OKbutton',[0 1 0 1],[299.25 73.5 51.75 16.5];...
'Partext',[1 0 1 0],[36.75 309.75 75 10.5];...
};
newpos=i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[443.25 219.75 372 329.25]);
