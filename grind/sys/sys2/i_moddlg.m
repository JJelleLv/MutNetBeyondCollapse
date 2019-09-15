function fig = i_moddlg(doclear, docheck, isvisible, doera)
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
global g_grind;
g=grindpath(2);
if isempty(g_grind)&&~strcmp(pwd, g)
    cd(g);
    fprintf('changed directory to %s\n',g);
end

fontsize=11;
font='Courier New';
if nargin == 0
    doclear = 0;
end

if nargin <2
    docheck = 1;
end

if nargin<3
    isvisible=true;
end
if nargin<4
    doera=false;
end
hfig=i_figno('dialog');
for i=1:10
    if ishandle(hfig+i)&&strcmp(get(hfig+i,'tag'),'GRIND Model')
        close(hfig+i);
    end

end

if ishandle(hfig)
    close(hfig);
    while ishandle(hfig)
        hfig=hfig+1;
    end

end

hfig= figure(hfig);
if ~isvisible
    set(hfig,'visible','off');
end


set(hfig, 'WindowButtonDown', '');
set(hfig, 'WindowButtonMotionFcn', '');
set(hfig,'Color',[0.941 0.941 0.941], ...
    'MenuBar','none', ...
    'Name','GRIND Model', ...
    'NumberTitle','off', ...
    'PaperPosition',[1296 12960 41472 31104], ...
    'PaperUnits','points', ...
    'ResizeFcn',@c_resize,...
    'Tag','GRIND Model', ...
    'ToolBar','none');
set(hfig,'KeypressFcn',@c_keydown)
set(hfig,	'Units','points');
pos=round(get(hfig,'position'));
pos(2)=pos(2)+25;
pos(3)=440;
pos(4)=320;
set(hfig,'position',pos);
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.941 0.941 0.941], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[11 281 106 16], ...
    'String','Name of model:', ...
    'Style','text', ...
    'Tag','namelabel');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'Callback',@c_updated, ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'FontName',font, ...
    'FontSize',fontsize, ...
    'Position',[130 284 181 16], ...
    'Style','edit', ...
    'Tag','FileName', ...
    'TooltipString','Enter here the name of the ini file');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_cd, ...
    'ListboxTop',0, ...
    'Position',[314 284 14 16], ...
    'String','...', ...
    'Tag','CDButton', ...
    'TooltipString','Select directory');
%328-14=314
if ~isempty(g_grind)&&(isfield(g_grind,'scheme') && length(g_grind.scheme)>10) && any(strncmp(g_grind.scheme,'%sym=',4))
    uicontrol('Parent',hfig, ...
        'Units','points', ...
        'BackgroundColor',[0.941 0.941 0.941], ...
        'HorizontalAlignment','left', ...
        'ListboxTop',0, ...
        'Position',[11 255 295 16], ...
        'String','Model equations (use vismod to change)', ...
        'Style','text', ...
        'Tag','modellabel');
    uicontrol('Parent',hfig, ...
        'Units','points', ...
        'Callback',@c_clearScheme, ...
        'ListboxTop',0, ...
        'Position',[200 260 70 16], ...
        'String','Clear scheme', ...
        'Tag','ClearSchemeButton', ...
        'TooltipString','Clear the vismod scheme');
    uicontrol('Parent',hfig, ...
        'Units','points', ...
        'BackgroundColor',[1 1 1], ...
        'Callback',@c_updated, ...
        'FontName',font, ...
        'FontSize',fontsize, ...
        'HorizontalAlignment','left', ...
        'ListboxTop',0, ...
        'Max',100, ...
        'Position',[11 143 317 116], ...
        'Enable','off', ...
        'Style','edit', ...
        'Tag','Model', ...
        'TooltipString','Enter here the differential equations (N''=..)');
else
    uicontrol('Parent',hfig, ...
        'Units','points', ...
        'BackgroundColor',[0.941 0.941 0.941], ...
        'HorizontalAlignment','left', ...
        'ListboxTop',0, ...
        'Position',[11 255 295 16], ...
        'String','Model equations (diffential/difference equations) ', ...
        'Style','text', ...
        'Tag','modellabel');
    uicontrol('Parent',hfig, ...
        'Units','points', ...
        'BackgroundColor',[1 1 1], ...
        'Callback',@c_updated, ...
        'FontName',font, ...
        'FontSize',fontsize, ...
        'HorizontalAlignment','left', ...
        'ListboxTop',0, ...
        'Max',100, ...
        'Position',[11 143 317 116], ...
        'Style','edit', ...
        'Tag','Model', ...
        'TooltipString','Enter here the differential equations (N''=..)');
end
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[0.941 0.941 0.941], ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Position',[11 120 295 16], ...
    'String','Parameters / Initial conditions / Commands', ...
    'Style','text', ...
    'Tag','parameterslabel');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'BackgroundColor',[1 1 1], ...
    'Callback',@c_updated, ...
    'FontName',font, ...
    'FontSize',fontsize, ...
    'HorizontalAlignment','left', ...
    'ListboxTop',0, ...
    'Max',100, ...
    'Position',[11 6 318 118], ...
    'Style','edit', ...
    'Tag','Parameters', ...
    'TooltipString','Assign default values to the parameters here');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_extract, ...
    'ListboxTop',0, ...
    'Position', [339 145 83 25], ...
    'String','Extract parameters', ...
    'Tag','ExtractButton', ...
    'TooltipString','Extract a list of parameters from the model eq.');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_cancel, ...
    'ListboxTop',0, ...
    'Position',[339 46 83 25], ...
    'String','Cancel', ...
    'Tag','CancelButton', ...
    'TooltipString','Cancel, don''t change/use model');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_ok, ...
    'ListboxTop',0, ...
    'Position', [339 13 83 25], ...
    'String','OK', ...
    'Tag','OKButton', ...
    'TooltipString','Save file and create model');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_openfile, ...
    'ListboxTop',0, ...
    'Position',[339 276 83 25], ...
    'String','Load model', ...
    'Tag','OpenFileButton', ...
    'TooltipString','Load an existing model');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'Callback',@c_newmodel, ...
    'ListboxTop',0, ...
    'Position',[339 111  83 25], ...
    'String','New model', ...
    'Tag','NewButton', ...
    'TooltipString','Clear all edit fields');
uicontrol('Parent',hfig, ...
    'Units','points', ...
    'ButtonDownFcn',')', ...
    'Callback','commands modelpanel', ...
    'ListboxTop',0, ...
    'Position',[339 78 83 25], ...
    'String','Help', ...
    'Tag','Helpbutton');
if nargout > 0, fig = hfig; end
c_init(hfig);
ud = get(hfig, 'Userdata');
ud.clear = doclear;
ud.check = docheck;
ud.era=doera;
ud.debug = 0;
set(hfig, 'Userdata', ud);
c_resize(hfig)

function c_makemodel(hobj,~)
hfig=getparentfig(hobj);
if ~isempty(hfig)
    h = get(hfig, 'Userdata');
    if ~isempty(h) && isfield(h, 'clear') && h.clear
        finishgrind;
        initgrind;
        evalin('base','clear global;');
        h = get(hfig, 'Userdata');
    else
        h.check = 1;
    end

else
    h.check  = 1;
end

if h.debug
    l = c_getmodel(hobj);
    i_makemodel(l);
    i_parcheck(1);
else
    try
        err = 0;
        s = '';
        try
            l = c_getmodel(hobj);
            i_makemodel(l);
            try
                [err, s] = i_parcheck(1);
                err = ~err;
            catch err1
                %       err1 = lasterror;
                s = striphtml(err1.message);
                err = 1;
            end

        catch err1
            %       err1 = lasterror;
            s = striphtml(err1.message);
            err = 1;
        end

        if ~err && h.check
            [err, s] = i_parcheck(1);
            err = ~err;
        end

        if err && h.check
            if i_myerrordlg(s, 'Error in model') == 2
                close(hfig);
                finishgrind;
                if ~isfield(h, 'era')
                    i_use(l.inifile, 1, 1, 0, h.era)
                else
                    i_use(l.inifile, 1, 1, 0, 1);
                end

            else
                model(l.inifile);
                return;
            end

        end

    catch err
        %err = lasterror;
        %    errordlg(striphtml([s err.message]));
        rethrow(err);
    end

end

function c_undoclearScheme(hobj,~)
global g_grind;
if isfield(g_grind, 'oldscheme')
    hfig=getparentfig(hobj);
    g_grind.scheme = g_grind.oldscheme;
    g_grind = rmfield(g_grind, 'oldscheme');
    h = get (hfig, 'Userdata');
    h.scheme = g_grind.scheme;
    set(hfig, 'Userdata', h);
    h=findobj(hfig,'Tag','ClearSchemeButton');
    set(h,'String','Clear scheme');
    set(h,'TooltipString','Clear the vismod scheme');
    set(h,'Callback',@c_clearScheme);
    h=findobj(hfig,'Tag','modellabel');
    set(h,'String','Model equations (use vismod to change)');
    h=findobj(hfig,'Tag','Model');
    set(h,'Enable','off');
end

function c_clearScheme(hobj,~)
global g_grind;
button = questdlg('Are you sure you want to delete the scheme?', ...
    'model','Yes', 'No','No');
if strcmp(button, 'Yes')
    hfig=getparentfig(hobj);
    g_grind.oldscheme = g_grind.scheme;
    h = get (hfig, 'Userdata');
    set(hfig,'Userdata',rmfield(h,'scheme'));
    h=findobj(hfig,'Tag','ClearSchemeButton');
    set(h,'String','Undo clear');
    set(h,'TooltipString','Undo clearing of the vismod scheme');
    set(h,'Callback',@c_undoclearScheme);
    h=findobj(hfig,'Tag','modellabel');
    set(h,'String','Model equations (diffential/difference equations) ');
    h=findobj(hfig,'Tag','Model');
    set(h,'Enable','on');
end

function c_resize(hobj,~)
hfig=getparentfig(hobj);
data={...
'namelabel',[1 0 1 0],[11 281.25 106 16];...
'FileName',[1 1 1 0],[130 284.25 180.75 16];...
'CDButton',[0 1  1 0],[313.75 284.25 14 16];...
'modellabel',[1 0 1 0],[11 255.25 295 16];...
'Model',[1 1 1 1],[11 142.125 316.75 117.125];...
'parameterslabel',[1 0  0 1],[11 119.125 295 16];...
'Parameters',[1 1 0 1],[11 6 317.75 117.125];...
'ExtractButton',[0 1 0 1],[338.75 145 83 25];...
'CancelButton',[0 1 0 1 ],[338.75 46 83 25];...
'OKButton',[0 1 0 1],[338.75 13 83 25];...
'OpenFileButton',[0 1 1 0],[338.75 276.25 83 25];...
'NewButton',[0 1 0 1],[338.75 111 83 25];...
'Helpbutton',[0 1 0 1],[338.75 78 83 25];...
};
i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[301.5 209.25 429.75 320.25]);
drawnow;
%distribute the height of the model and parameters panel
movetop=i_resizedlg(hfig,'vertdistrib',{'Model','Parameters'});
h=findobj(hfig,'tag','parameterslabel');
posparlab=get(h,'position');
if numel(posparlab)==4
   posparlab(2) = posparlab(2) + movetop(1);
   set(h,'position',posparlab);
end

function c_init(hobj,~)
global g_grind;
f = getparentfig(hobj);
ud.clear = 0;
ud.updated = 0;
ud.era = 1;
ud.scheme = {};
set(f, 'Userdata', ud);
if ~isempty(g_grind) && isfield(g_grind, 'model') && ~isempty(g_grind.model)
    %      global g_grind;
    if isfield(g_grind, 'scheme')  && ~isempty(g_grind.scheme)
        ud.scheme = g_grind.scheme;
        set(f, 'Userdata', ud);
    end

    h = findobj(f,'Tag','FileName');
    set(h, 'String', g_grind.inifile);
    h = findobj(f,'Tag','Model');
    for i = 1:length(g_grind.model)
        if ~isempty(g_grind.model{i})
            g_grind.model{i} = strtrim(g_grind.model{i});
        end

        if isempty(g_grind.model{i})
            g_grind.model{i} = ' ';
        else
        end

    end

    set(h, 'String', g_grind.model);
    h = findobj(f,'Tag','Parameters');
    for i = 1:length(g_grind.commands)
        if isempty(g_grind.commands{i})
            g_grind.commands{i} = ' ';
        end

    end

    set(h, 'String', g_grind.commands);
end

function res=c_getmodel(hobj,~)
global g_grind;
f=getparentfig(hobj);
h = findobj(f, 'Tag', 'Model');
r.model =  i_memo2cell(get(h, 'String'));
h = findobj(f,'Tag','Parameters');
r.commands = i_memo2cell(get(h, 'String'));
h = findobj(f, 'Tag', 'FileName');
r.inifile = strtrim(get(h, 'String'));
h = get(f, 'Userdata');
r.scheme = {};
if ~isempty(h) && isfield(h, 'scheme')
    r.scheme = h.scheme;
end

if ~isempty(h) && isfield(h, 'updated') && h.updated
    if isempty(r.inifile)
        s = inputdlg({'Enter file name'}, 'File name', 1, {'curr_mod.ini'});
        r.inifile = s{1};
    end

    g_grind.model = r.model;
    g_grind.commands = r.commands;
    g_grind.inifile = r.inifile;
    g_grind.scheme = r.scheme;
    savemodel(char(r.inifile));
end

if nargout == 1
    res = r;
else
    g_grind.model = r.model;
    g_grind.commands = r.commands;
    g_grind.inifile = r.inifile;
    g_grind.scheme = r.scheme;
end

function c_newmodel(hobj,~)
f=getparentfig(hobj);
savechanges(f);
setupdated(0,f);
s = ' ';
h = findobj(f, 'Tag', 'Model');
set(h, 'String', s);
h = findobj(f, 'Tag', 'Parameters');
set(h, 'String', s);
h = findobj(f, 'Tag', 'FileName');
set(h, 'String', s);
function c_updated(hobj,~)
setupdated(1,getparentfig(hobj))
function c_extract(hobj,~)
f=getparentfig(hobj);
%   global g_grind;
h = findobj(f, 'Tag', 'Model');
ggrind=i_analysemodel(i_memo2cell(get(h, 'String')));
h = findobj(f, 'Tag', 'Parameters');
s1 = i_memo2cell(get(h, 'String'));
hh = length(s1);
while (hh > 0) && isempty(strtrim(s1{hh}))
    hh = hh - 1;
end

s = cell(1,size(ggrind.pars, 2) + 10+hh);
s2 = cell(hh, 1);
k1=0;
for i = 1:hh
    f1 = strfind(s1{i},'=');
    if ~isempty(f1)
        k1 = k1 + 1;
        s2{k1} = strtrim(s1{i}(1:f1(1) - 1));
    end

    s{i} = s1{i};
end

s2 = s2(1:k1);
done = 0;
for i = 1:size(ggrind.pars, 2)
    f1 = 0;
    for j = 1:length(s2)
        if strcmp(ggrind.pars{i}, s2{j})
            f1 = 1;
            hh = hh - 1;
        end

    end

    if ~f1
        done = 1;
        s{i + hh} = [char(ggrind.pars{i}) '='];
    end

end

j = length(s);
while (j > 0) &&  isempty(s{j})
    j = j - 1;
end

s = s(1:j);
set(h, 'String', s);
if ~done
    msgbox('No new parameters found');
else
    setupdated(1,f)
end

function c_keydown(hobj,~)
f=getparentfig(hobj);
ch = get(f, 'CurrentCharacter');
if ~isempty(ch)
    if (ch == char(4)) %control - d sets in debug mode
        ud = get(f, 'userdata');
        ud.debug = ~ud.debug;
        set(f, 'userdata', ud);
        if ud.debug
            msgbox('Ctrl-D pressed: Debug mode','model')
        else
            msgbox('Ctrl-D pressed: Normal mode: error messages displayed','model');
        end

    elseif ch == char(8) %ctrl - h
        c_help(hobj);
        
    elseif ch == char(3) %ctrl - c
        if strcmp('Yes', questdlg('OK to cancel?','Model','Yes','No','No'))
            c_cancel(hobj);
        end

    elseif ch == char(13) %enter
        c_ok(hobj);
    end

end

function c_ok(hobj,~)
f=getparentfig(hobj);
set(f, 'DeleteFcn', @c_makemodel );
delete(f);
function c_cancel(hobj,~)
global g_grind;
f=getparentfig(hobj);
if isempty(g_grind) || (isfield(g_grind, 'model') && isempty(g_grind.model))
    clear g_grind;
end

delete(f);
function c_cd(hobj,~)
f=getparentfig(hobj);
setupdated(1,f)
h=findobj(f,'Tag','FileName');
curfil = get(h, 'String');
if isempty(curfil)
    curfil = '*.ini';
end

if ~strcontains(curfil, '.ini')
    curfil = [curfil '.ini'];
end

[filename, pathname] = uiputfile(curfil,'Select file name to save model');
if ~isempty(filename) && ischar(filename)
    h=findobj(f,'Tag','FileName');
    set(h, 'String', char(filename));
end

if ~isempty(pathname) && ischar(pathname)
    cd(pathname);
end

function c_openfile(hobj,~)
f=getparentfig(hobj);
savechanges(f);
[filename, pathname] = uigetfile('*.ini','Select file name');
h=findobj(f,'Tag','FileName');
if (filename ~= 0)
    set(h, 'String', char(filename));
    cd(pathname);
    c_loadfile(hobj);
end
function c_loadfile(hobj,~)
f=getparentfig(hobj);
setupdated(0,f);
h = findobj(f, 'Tag', 'FileName');
l_inifile = strtrim(get(h, 'String'));
[l_model, l_commands, l_scheme] = i_loadinifile(l_inifile);
if ~isempty(l_model)
    h = findobj(f, 'Tag','Model');
    set(h, 'String', l_model);
    h = findobj(f, 'Tag', 'Parameters');
    set(h, 'String', l_commands);
    ud = get(f, 'Userdata');
    ud.scheme = l_scheme;
    set(f, 'Userdata', ud);
else
    h = findobj(f, 'Tag', 'Model');
    set(h, 'String', {'', ''});
    h = findobj(f, 'Tag', 'Parameters');
    set(h, 'String', {'', ''});
    ud = get(f, 'Userdata');
    ud.scheme = {};
    set(f, 'Userdata', ud);
end



function setupdated(value,f)
ud = get(f, 'Userdata');
ud.updated = value;
set(f, 'Userdata', ud);
return;

function savechanges(f)
global g_grind;
ud = get(f, 'Userdata');
if ud.updated
    h = findobj(f, 'Tag', 'FileName');
    l_inifile = get(h, 'String');
    if ~isempty(l_inifile)&& strcmp('Yes',questdlg(['Do you want to save changes in ' l_inifile '?'], ...
            'Save changes','Yes','No','Yes'))
        oldgmodel = g_grind.model;
        oldgcommands = g_grind.commands;
        oldginifile = g_grind.inifile;
        h = findobj(f, 'Tag', 'Model');
        g_grind.model = transpose(get(h, 'String'));
        h=findobj(f, 'Tag', 'Parameters');
        g_grind.commands = transpose(get(h, 'String'));
        savemodel(l_inifile, 1);
        g_grind.model = oldgmodel;
        g_grind.commands = oldgcommands;
        g_grind.inifile = oldginifile;
    end

end


