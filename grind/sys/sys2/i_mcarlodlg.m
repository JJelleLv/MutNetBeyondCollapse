function  res=i_mcarlodlg(flag)
evalin('base','global g_mcarlo;');
global g_grind g_mcarlo;
if nargin < 1
    flag = 'init';
end

distrs={'Uniform','Normal','TruncNormal','LogNormal','paranal','hysteresis'};
switch flag
    case 'init'
        if ~isempty(g_mcarlo) && (~isfield(g_mcarlo, 'inifile') || strcmp(g_grind.inifile, g_mcarlo.inifile))
            if ~isfield(g_mcarlo, 'allpars')
                ud=i_mcarlodlg('initstruc');
            else
                ud = g_mcarlo.allpars;
            end

            for i = 1:length(ud)
                ud(i).value = evalin('base', ud(i).name);
                if ~isnan(ud(i).range)
                    ud(i).min = ud(i).value .* (1 - ud(i).range);
                    ud(i).max = ud(i).value .* (1 + ud(i).range);
                end

            end

        else
            g_mcarlo = [];
            ud=i_mcarlodlg('initstruc');
        end

        pars = cell(1, length(ud));
        for i = 1:length(ud)
            pars{i} = ud(i).name;
        end

        
        hfig = figure('Color',[0.941 0.941 0.941], ...
            'MenuBar','none', ...
            'Name','Monte-Carlo sensitivity analysis', ...
            'NumberTitle','off', ...
            'PaperPosition',[18 180 576 432], ...
            'PaperUnits','points', ...
            'Position',[336 313 567 403], ...
            'Tag','Fig1', ...
            'ResizeFcn',@resize,...
            'ToolBar','none',...
            'CreateFcn',@(h,evnt)movegui(h, 'center'));
        
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[1 1 1], ...
            'Callback',@c_listboxclick, ...
            'Position',[16.5 20.25 126 244.5], ...
            'String',pars, ...
            'Style','listbox', ...
            'Tooltipstring','Select parameters for sensitivity analysis',...
            'Tag','ParList', ...
            'Value', 1);
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[15 273.75 153 15], ...
            'String','Select parameters for sensitivity analysis:', ...
            'Style','text', ...
            'Tag','SelParText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'Callback',@c_selectedcheck, ...
            'ListboxTop',0, ...
            'Position',[157.5 230.75 72 21], ...
            'String','Selected?', ...
            'Tooltipstring','Check this box to select the current parameter',...
            'Style','checkbox', ...
            'Tag','SelectedCheck');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'Callback',@c_SelectAllclick, ...
            'ListboxTop',0, ...
            'Position',[342.5 230.75 71.25 26.25], ...
            'String','Select all', ...
            'Tooltipstring','Push button to select all parameters (using a certain distribution for the range)',...
            'Tag','SelectallButtom');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[1 1 1], ...
            'ListboxTop',0, ...
            'Position',[258.75 200 93.75 18], ...
            'String',' ', ...
            'Callback',@c_elemclick, ...
            'Style','popupmenu', ...
            'Visible','off',...
            'Tooltipstring','Select element of vector/matrix parameter',...
            'Tag','ElemList', ...
            'Value', 1);
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[157.5 200 91.5 18], ...
            'String','Element vector/matrix', ...
            'Style','text', ...
            'Visible','off',...
            'Tag','ElemText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[1 1 1], ...
            'Callback',@c_distrclick, ...
            'ListboxTop',0, ...
            'Position',[258.75 170.5 93.75 18], ...
            'String',distrs, ...
            'Tooltipstring','Select distribution for the Monte-Carlo analysis',...
            'Style','popupmenu', ...
            'Tag','DistrPopup', ...
            'Value', 1);
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[157.5 175.5 45 12], ...
            'String','Distribution', ...
            'Style','text', ...
            'Tag','DistrText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[1 1 1], ...
            'Callback',@c_minedited, ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[258.75 138.75 93.75 18], ...
            'String','', ...
            'Tooltipstring','Edit the minimum value of the distribution',...
            'Style','edit', ...
            'Tag','MinEdit');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[157.5 138.75 65 15], ...
            'String','Minimum:', ...
            'Style','text', ...
            'Tag','MinText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[157.5 102 65 15], ...
            'String','Maximum', ...
            'Style','text', ...
            'Tag','MaxText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[1 1 1], ...
            'Callback',@c_maxedited, ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[258.75 102 93.75 18], ...
            'String','', ...
            'Tooltipstring','Edit the maximum value of the distribution',...
            'Style','edit', ...
            'Tag','MaxEdit');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[1 1 1], ...
            'Callback',@c_rangeedited, ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[258.75 69 93.75 18], ...
            'String','', ...
            'Tooltipstring','Edit the range of the distribution',...
            'Style','edit', ...
            'Tag','RangeEdit');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[157.5 69 80 15], ...
            'String','Relative range:', ...
            'Style','text', ...
            'Tag','RangeText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'BackgroundColor',[0.941 0.941 0.941], ...
            'HorizontalAlignment','left', ...
            'ListboxTop',0, ...
            'Position',[158.25 260.25 200.25 15], ...
            'String','Parameter', ...
            'Style','text', ...
            'Tag','ParText');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'Callback',@c_OKclick, ...
            'ListboxTop',0, ...
            'Position',[176.25 17.25 71.25 26.25], ...
            'String','OK', ...
            'Tooltipstring','Save changes and close',...
            'Tag','OKButton');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'Callback',@c_Cancelclick, ...
            'ListboxTop',0, ...
            'Position',[252 17.25 71.25 26.25], ...
            'String','Cancel', ...
            'Tooltipstring','Close without saving changes',...
            'Tag','CancelButton');
        uicontrol('Parent',hfig, ...
            'Units','points', ...
            'Callback',@c_Helpclick, ...
            'ListboxTop',0, ...
            'Position',[327.75 17.25 71.25 26.25], ...
            'String','Help', ...
            'Tooltipstring','Display help on mcarlo',...
            'Tag','HelpButton');


        set(hfig, 'userdata', ud);
        c_listboxclick(hfig);
    case 'initstruc'
        res = struct('name',g_grind.pars,'descr','','value',0,'selected',...
        false,'range',0,'min',0,'max',0,'sd',0,'distr','Uniform');
        for i = 1:length(g_grind.pars)
            res(i).descr = par('-d', res(i).name);
            if isempty(res(i).descr)
                res(i).descr = 'Parameter';
            end

            res(i).value = evalin('base', res(i).name);
            res(i).selected = false(size(res(i).value));
            res(i).range = 0.1 + double(res(i).selected);
            res(i).min = res(i).value .* (1 - res(i).range);
            res(i).max = res(i).value .* (1 + res(i).range);
            res(i).sd = res(i).value .* res(i).range;
        end

        npar = length(res);
        if g_grind.statevars.vector
            for i = 1:length(g_grind.statevars.vectnames)
                p.name = g_grind.statevars.vectnames{i};
                p.descr = sprintf('Initial condition of %s', p.name);
                p.value = evalin('base', p.name);
                p.selected = false(size(p.value));
                p.range = 0.1 + double(p.selected);
                p.min = p.value .* (1 - p.range);
                p.max = p.value .* (1 + p.range);
                p.sd = p.value .* p.range;
                p.distr = 'Uniform';
                res(npar + i) = p;
            end
        else
            for i = 1:g_grind.statevars.dim
                p.name = g_grind.statevars.names{i};
                p.descr = sprintf('Initial condition of %s', p.name);
                p.value = evalin('base', p.name);
                p.selected = false;
                p.range = 0.1;
                p.min = p.value .* (1 - p.range);
                p.max = p.value .* (1 + p.range);
                p.sd = p.value .* p.range;
                p.distr = 'Uniform';
                res(npar + i) = p;
            end

        end
end

function resize(hobj,~)
hfig=getparentfig(hobj);
data={'ParList',[1 1 1 1],[16.5 20.25 126 244.5];...
'SelParText',[1 0 1 0],[15 273.75 153 15];...
'SelectedCheck',[0 1 1 0],[157.5 230.75 72 21];...
'SelectallButtom',[0 1 1 0],[342.5 230.75 71.25 26.25];...
'ElemList',[1 0 1 0],[258.75 200 93.75 18];...
'ElemText',[1 0 1 0],[157.5 200 91.5 18];...
'DistrPopup',[0 1 0 1],[258.75 170.5 93.75 18];...
'DistrText',[0 1 0 1],[157.5 175.5 45 12];...
'MinEdit',[0 1 0 1],[258.75 138.75 93.75 18];...
'MinText',[0 1 0 1],[157.5 138.75 65 15];...
'MaxText',[0 1 0 1],[157.5 102 65 15];...
'MaxEdit',[0 1 0 1],[258.75 102 93.75 18];...
'RangeEdit',[0 1 0 1],[258.75 69 93.75 18];...
'RangeText',[0 1 0 1],[157.5 69 80 15];...
'ParText',[0 1 1 0],[158.25 260.25 200.25 15];...
'OKButton',[0 1 0  1],[176.25 17.25 71.25 26.25];...
'CancelButton',[0 1 0 1],[252 17.25 71.25 26.25];...
'HelpButton',[0 1 0 1],[327.75 17.25 71.25 26.25];...
};
i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[633.75 192.75 425.25 302.25]);


function c_distrclick(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h=findobj(hfig,'tag','ParList');
parnr = get(h, 'value');
h=findobj(hfig,'tag','DistrPopup');
v = get(h, 'value');
distrs=get(h,'string');
ud(parnr).distr = distrs{v};
set(hfig, 'userdata', ud);
c_listboxclick(hfig);
function c_listboxclick(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h=findobj(hfig,'tag','ParList');
parnr = get(h, 'value');
if isempty(parnr)
    parnr = 1;
end

h3=findobj(hfig,'tag','ParText');
set(h3,'string',sprintf('%s:  %s = %g',ud(parnr).descr,ud(parnr).name,mean(ud(parnr).value(:))));
h=findobj(hfig,'tag', 'RangeEdit');
set(h, 'string', num2str(ud(parnr).range(1)));
h=findobj(hfig,'tag', 'MaxEdit');
set(h, 'string', num2str(ud(parnr).max(1)));
h=findobj(hfig,'tag','MinEdit');
set(h, 'string', num2str(ud(parnr).min(1)));
h=findobj(hfig,'tag', 'SelectedCheck');
set(h, 'value', double(ud(parnr).selected(1)));
h1=findobj(hfig,'Tag','ElemList');
h2=findobj(hfig,'Tag','ElemText');
if numel(ud(parnr).value) > 1
    set(h3,'string',sprintf('%s:  mean(%s) = %g',ud(parnr).descr,ud(parnr).name,mean(ud(parnr).value(:))));
    set(h1, 'visible', 'on');
    set(h2, 'visible', 'on');
    elems = cell(1, numel(ud(parnr).value) + 1);
    if size(ud(parnr).value, 2) > 1
        elems{1} = sprintf('%s(1:%d,1:%d)',ud(parnr).name,size(ud(parnr).value));
    else
        elems{1} = sprintf('%s(1:%d)', ud(parnr).name, length(ud(parnr).value));
    end

    k = 2;
    if size(ud(parnr).value, 2) > 1
        for i = 1:size(ud(parnr).value, 1)
            for j = 1:size(ud(parnr).value, 2)
                elems{k} = sprintf('%s(%d,%d)',ud(parnr).name,i,j);
                k = k + 1;
            end

        end

    else
        for k = 1:size(ud(parnr).value, 1)
            elems{k + 1} = sprintf('%s(%d)', ud(parnr).name, k);
        end

    end

    set(h1, 'value', 1);
    set(h1, 'string', elems);
else
    set(h1, 'visible', 'off');
    set(h1, 'value', 1);
    set(h2, 'visible', 'off');
end


h=findobj(hfig,'tag','DistrPopup');
distrs=get(h,'string');
for i = 1:length(distrs)
    if strcmpi(distrs{i}, ud(parnr).distr)
        set(h, 'value', i);
    end

end

if ~any(strcmpi({'uniform','paranal','hysteresis'},ud(parnr).distr))
    h=findobj(hfig,'tag', 'MaxEdit');
    set(h, 'visible', 'off');
    h=findobj(hfig,'tag','MinEdit');
    set(h, 'Tooltipstring','Edit the standard deviation of the distribution');
    h=findobj(hfig,'tag', 'MaxText');
    set(h, 'visible', 'off');
    h=findobj(hfig,'tag', 'MinEdit');
    set(h, 'string', num2str(ud(parnr).sd(1)));
    h=findobj(hfig,'tag', 'MinText');
    set(h,'string','Stand.dev.');
    h=findobj(hfig,'tag', 'RangeText');
    set(h,'string','CV (sd/mean)');
    h=findobj(hfig,'tag','RangeEdit');
    set(h, 'Tooltipstring','Edit the coefficient of variation (CV)');

elseif strcmpi('uniform', ud(parnr).distr)
    h=findobj(hfig,'tag','MinEdit');
    set(h, 'Tooltipstring','Edit the minimum value of the distribution');
    h=findobj(hfig,'tag', 'MaxEdit');
    set(h, 'visible', 'on');
    set(h, 'Tooltipstring','Edit the maximum value of the distribution');
    h=findobj(hfig,'tag', 'MaxText');
    set(h,'string','Maximum:');
    set(h, 'visible', 'on');
    h=findobj(hfig,'tag', 'MinText');
    set(h,'string','Minimum:');
    h=findobj(hfig,'tag', 'RangeText');
    set(h,'string','Relative range:')
    h=findobj(hfig,'tag','RangeEdit');
    set(h, 'Tooltipstring','Edit the relative range of the distribution');
else
    h=findobj(hfig,'tag', 'MinText');
    set(h,'string','[Start, End]:');
    h=findobj(hfig,'tag', 'MinEdit');
    set(h, 'Tooltipstring','Edit the range of paranal')
    if isfield(ud(parnr),'minmax')
        set(h,'string',sprintf('[%g %g]',ud(parnr).minmax));
    else
        set(h,'string',sprintf('[%g %g]',ud(parnr).min(1),ud(parnr).max(1)));
        ud(parnr).minmax=[ud(parnr).min(1),ud(parnr).max(1)];
    end

    h=findobj(hfig,'tag', 'MaxEdit');
    set(h, 'Tooltipstring','Edit the number of steps for paranal')
    set(h, 'visible', 'on');
    if isfield(ud(parnr),'steps')
        set(h,'string',sprintf('%g',ud(parnr).steps));
    else
        set(h,'string','40');
        ud(parnr).steps=40;
    end

    h=findobj(hfig,'tag', 'MaxText');
    set(h,'string','Number of steps:');
    set(h, 'visible', 'on');
    h=findobj(hfig,'tag', 'RangeText');
    set(h,'string','[Stabilizing, Writing]:')    
    set(h, 'Tooltipstring','Edit the number of steps of stabilizing and writing in paranal')

    h=findobj(hfig,'tag', 'RangeEdit');
    if isfield(ud(parnr),'nstabilwrite')
        set(h,'string',sprintf('[%g %g]',ud(parnr).nstabilwrite));
    else
        set(h,'string','[600 300]');
        ud(parnr).nstabilwrite=[600, 300];
    end

    set(hfig,'userdata',ud);
end

function c_selectedcheck(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h=findobj(hfig,'tag','ParList');
parnr = get(h, 'value');
h=findobj(hfig,'Tag','ElemList');
elemnr =  get(h, 'value');
h=findobj(hfig,'tag', 'SelectedCheck');
if elemnr == 1
    ud(parnr).selected =  false(size(ud(parnr).selected)) + get(h, 'value')>0;
else
    ud(parnr).selected(elemnr - 1) =  get(h, 'value')>0;
end

set(hfig, 'userdata', ud);
function c_rangeedited(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h=findobj(hfig,'tag','ParList');
parnr = get(h, 'value');
h=findobj(hfig,'Tag','ElemList');
elemnr =  get(h, 'value');
h=findobj(hfig,'tag', 'RangeEdit');
if any(strcmpi({'paranal','hysteresis'},ud(parnr).distr))
    if elemnr  > 1
        error('GRIND:mcarlo:notimplemented','Vector parameters not yet implemented for paranal');
    else
        ud(parnr).nstabilwrite = str2num(get(h, 'string'));  %#ok<ST2NM>
        set(hfig, 'userdata', ud);
    end

else
    %if ~isnan(ud(parnr).range)
    if elemnr == 1
        ud(parnr).range = zeros(size(ud(parnr).value)) + str2double(get(h, 'string'));
        ud(parnr).min = ud(parnr).value .* (1 - ud(parnr).range);
        ud(parnr).max = ud(parnr).value .* (1 + ud(parnr).range);
        ud(parnr).sd = ud(parnr).value .* ud(parnr).range;
    else
        elemnr = elemnr - 1;
        ud(parnr).range(elemnr) = str2double(get(h, 'string'));
        ud(parnr).min(elemnr) = ud(parnr).value(elemnr) .* (1 - ud(parnr).range(elemnr));
        ud(parnr).max(elemnr) = ud(parnr).value(elemnr) .* (1 + ud(parnr).range(elemnr));
        ud(parnr).sd(elemnr) = ud(parnr).value(elemnr) .* ud(parnr).range(elemnr);
    end

    %end

    set(hfig, 'userdata', ud);
    if strcmpi(ud(parnr).distr, 'uniform')
        h=findobj(hfig,'tag', 'MaxEdit');
        set(h, 'string', num2str(ud(parnr).max(elemnr)));
        h=findobj(hfig,'tag','MinEdit');
        set(h, 'string', num2str(ud(parnr).min(elemnr)));
    else
        h=findobj(hfig,'tag', 'MinEdit');
        set(h, 'string', num2str(ud(parnr).sd(elemnr)));
    end

end
function c_maxedited(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h=findobj(hfig,'tag','ParList');
parnr = get(h, 'value');
h=findobj(hfig,'Tag','ElemList');
elemnr =  get(h, 'value');
h=findobj(hfig,'tag', 'MaxEdit');
if any(strcmpi({'paranal','hysteresis'},ud(parnr).distr))
    if elemnr  > 1
        error('GRIND:mcarlo:notimplemented','Vector parameters not yet implemented for paranal');
    else
        ud(parnr).steps = str2num(get(h, 'string'));  %#ok<ST2NM>
        set(hfig, 'userdata', ud);
    end

else
    if elemnr == 1
        ud(parnr).max = zeros(size(ud(parnr).value)) + str2double(get(h, 'string'));
        if ud(parnr).value == 0
            ud(parnr).range = NaN;
        else
            r1 = iif(ud(parnr).value == 0, 0, 1 - ud(parnr).min ./ ud(parnr).value);
            r2 = iif(ud(parnr).value == 0, 0, ud(parnr).max ./ ud(parnr).value - 1);
            if r2 ~= r1
                ud(parnr).range = NaN;
            else
                ud(parnr).range = r1;
            end

        end

    else
        elemnr = elemnr - 1;
        ud(parnr).max(elemnr) = str2double(get(h, 'string'));
        if ud(parnr).value(elemnr) == 0
            ud(parnr).range(elemnr) = NaN;
        else
            r1 = iif(ud(parnr).value(elemnr) == 0, 0, 1 - ud(parnr).min(elemnr) ./ ud(parnr).value(elemnr));
            r2 = iif(ud(parnr).value(elemnr) == 0, 0, ud(parnr).max(elemnr) ./ ud(parnr).value(elemnr) - 1);
            if r2 ~= r1
                ud(parnr).range(elemnr) = NaN;
            else
                ud(parnr).range(elemnr) = r1;
            end

        end

    end

    set(hfig, 'userdata', ud);
    h=findobj(hfig,'tag', 'RangeEdit');
    set(h, 'string', num2str(ud(parnr).range(elemnr)));
end
function c_minedited(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h=findobj(hfig,'tag','ParList');
parnr = get(h, 'value');
h=findobj(hfig,'Tag','ElemList');
elemnr =  get(h, 'value');
h=findobj(hfig,'tag', 'MinEdit');
if any(strcmpi({'paranal','hysteresis'},ud(parnr).distr))
    if elemnr  > 1
        error('GRIND:mcarlo:notimplemented','Vector parameters not yet implemented for paranal');
    else
        ud(parnr).minmax = str2num(get(h, 'string'));  %#ok<ST2NM>
        set(hfig, 'userdata', ud);
    end

else
    if elemnr == 1
        if ~strcmpi(ud(parnr).distr, 'uniform')
            ud(parnr).sd = zeros(size(ud(parnr).value)) + str2double(get(h, 'string'));
            if ud(parnr).value > 0
                ud(parnr).range = iif(ud(parnr).value == 0, 0, ud(parnr).sd ./ ud(parnr).value);
            else
                ud(parnr).range = NaN;
            end

        else
            ud(parnr).min = zeros(size(ud(parnr).value)) + str2double(get(h, 'string'));
            if ud(parnr).value == 0
                ud(parnr).range = NaN;
            else
                r1 = iif(ud(parnr).value == 0, 0, 1 - ud(parnr).min ./ ud(parnr).value);
                r2 = iif(ud(parnr).value == 0, 0, ud(parnr).max ./ ud(parnr).value - 1);
                if r2 ~= r1
                    ud(parnr).range = NaN;
                else
                    ud(parnr).range = r1;
                end

            end

        end

    else
        elemnr = elemnr - 1;
        if ~strcmpi(ud(parnr).distr, 'uniform')
            ud(parnr).sd(elemnr) = str2double(get(h, 'string'));
            if ud(parnr).value{elemnr} > 0
                ud(parnr).range(elemnr) = iif(ud(parnr).value(elemnr) == 0, 0, ud(parnr).sd(elemnr) ./ ud(parnr).value(elemnr));
            else
                ud(parnr).range(elemnr) = NaN;
            end

        else
            ud(parnr).min(elemnr) = str2double(get(h, 'string'));
            if ud(parnr).value{elemnr} > 0
                r1 = iif(ud(parnr).value(elemnr) == 0, 0, 1 - ud(parnr).min(elemnr) ./ ud(parnr).value(elemnr));
                r2 = iif(ud(parnr).value(elemnr) == 0, 0, ud(parnr).max(elemnr) ./ ud(parnr).value(elemnr) - 1);
                if r2 ~= r1
                    ud(parnr).range(elemnr) = NaN;
                else
                    ud(parnr).range(elemnr) = r1;
                end

            else
                ud(parnr).range(elemnr) = NaN;
            end

        end

    end

    set(hfig, 'userdata', ud);
    h=findobj(hfig,'tag', 'RangeEdit');
    set(h, 'string', num2str(ud(parnr).range(elemnr)));
end

function c_elemclick(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
h1=findobj(hfig,'tag','ParList');
parnr = get(h1, 'value');
h=findobj(hfig,'tag','ElemList');
elemnr = get(h, 'value');
elems = get(h, 'string');
h=findobj(hfig,'tag','ParText');
if elemnr == 1
    set(h,'string',sprintf('%s:  mean(%s) = %g',ud(parnr).descr,ud(parnr).name,mean(ud(parnr).value(:))));
else
    set(h,'string',sprintf('%s:  %s = %g',ud(parnr).descr,elems{elemnr},ud(parnr).value(elemnr-1)));
end

if elemnr > 1
    elemnr = elemnr - 1;
end

h=findobj(hfig,'tag', 'RangeEdit');
set(h, 'string', num2str(ud(parnr).range(elemnr)));
h=findobj(hfig,'tag', 'MaxEdit');
set(h, 'string', num2str(ud(parnr).max(elemnr)));
h=findobj(hfig,'tag','MinEdit');
if strcmpi(ud(parnr).distr, 'uniform')
    set(h, 'string', num2str(ud(parnr).min(elemnr)));
else
    set(h, 'string', num2str(ud(parnr).sd(elemnr)));
end

h=findobj(hfig,'tag', 'SelectedCheck');
set(h, 'value', double(ud(parnr).selected(elemnr)));
function c_SelectAllclick(hobject,~)
hfig=getparentfig(hobject);
ud = get(hfig, 'userdata');
Ans=inputdlg({'Rel.range/CV for all parameters (empty=no update):', ...
    'Distribution (Normal,Uniform,LogNormal,empty= no update)'},'Select all',1,{'0.1','Uniform'});
if ~isempty(Ans)
    ran = str2double(Ans{1});
    if ~isempty(ran) && (ran > 0) && (ran < 1)
        for i = 1:length(ud)
            ud(i).range = ran + zeros(size(ud(i).value));
            ud(i).min = ud(i).value .* (1 - ud(i).range);
            ud(i).max = ud(i).value .* (1 + ud(i).range);
        end

    end

    distr = Ans{2};
    if ~isempty(distr)
        for i = 1:length(ud)
            ud(i).distr = distr;
        end

    end

    for i = 1:length(ud)
        ud(i).selected = true(size(ud(i).value));
    end

    set(hfig, 'userdata', ud);
    c_listboxclick(hfig);
end

function c_OKclick(hobject,~)
global g_mcarlo;
hfig=getparentfig(hobject);
g_mcarlo.allpars = get(hfig, 'userdata');
close(hfig);
function c_Helpclick(~,~)
disp('Help not yet implemented');
function c_Cancelclick(hobject,~)
hfig=getparentfig(hobject);
close(hfig);

