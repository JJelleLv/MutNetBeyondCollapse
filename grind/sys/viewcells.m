%VIEWCELLS   view a matrix variable
%  The values of a matrix variable are shown as a movie of 
%  changing colors. A dialog screen (<a href="matlab:help replayall">replayall</a>) is opened to move forward of backwards.
%  This command is commonly used for spatial explicit models and 
%  cellular automats.
%
%   Usage:
%   VIEWCELLS - shows all vector/matrix variables.    
%   VIEWCELLS NDAYS - runs the model for NDAYS days.
%   VIEWCELLS -out (-o) - change the variables to be plotted.
%   VIEWCELLS -out plotno var [siz1 siz2] - sets the output in a command line:
%   plotno  = number of plot, var = the state or auxiliary variable, [siz1 siz2] 
%   is the size of the variable.
%
%   VIEWCELLS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'ndays' [number>0] - number of days for running
%   VIEWCELLS('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - add: continue with the previous run.
%     '-d' - set back to default output
%     '-l' - list the output
%     '-n' - nocheck: do not check for changed option, just show last results.
%     '-o' - opens dialog box to select the plot.
%     '-o' plotno var [siz1 siz2] - sets the output in a command line:
%        plotno  = number of plot, var = the state or auxiliary variable, [siz1 siz2] 
%        is the size of the variable.
%     '-p' - makes it possible to replay a paranal
%     '-p1' or '-pa' - show the results of the last run of PARANAL.
%     '-p2' or '-paranal2' - show the results of the last two runs of PARANAL
%     '-r' - rerun the model always.
%     '-s' - silent (-s) do not plot the results.
%
%
%   See also vectplot, setcolormap, replayall, paranal
%
%   Reference page in Help browser:
%      <a href="matlab:commands('viewcells')">commands viewcells</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function viewcells(varargin)
global t g_grind g_Y g_t;
%opengl software
if (nargin > 0) && ischar(varargin{1}) && strncmpi(varargin{1}, '-o', 2)
    i_parcheck;
    if nargin == 1
        i_viewcellsdlg;
        return;
    end
    figno = varargin{2};
    if ischar(figno)
        figno=str2double(figno);
    end
    if isnan(figno)||numel(figno)~=1
        figno=1;
        h=1;
    else
        h=0;
    end
    fun = varargin{3-h};
    if length(varargin) > 3-h
        siz = i_checkstr(varargin{4-h});
    else
        siz=evalin('base', sprintf('size(%s);', fun));
    end
    g_grind.viewcells.currno = figno;
    g_grind.viewcells.vars{figno}.name = fun;
    g_grind.viewcells.vars{figno}.dims = siz;
    return;
end
fieldnams={'ndays', 'n>0', 'number of days for running',g_grind.ndays}';
res = i_time_options(i_time_options(i_parseargs(fieldnams,'ndays',...
    {'-r','-a','-n','-s','-p2|-paranal2','-p1|-pa','-d','-l','-p','-o'},varargin)));
i_parcheck;
if ~isfield(g_grind, 'viewcells')
    g_grind.viewcells.currno = 1;
    for i = 1:length(g_grind.statevars.vectnames)
        if g_grind.statevars.dims{i}.dim1 * g_grind.statevars.dims{i}.dim2 > 1
            viewcells('-out', i, g_grind.statevars.vectnames{i}, [g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2]);
        end
    end
end
i = 1;
while i <= length(res.opts)
    if (nargin == 0) || strncmpi(res.opts(i), '-d', 2) %-defaults
        if isfield(g_grind, 'viewcells')
            g_grind.viewcells.currno = 1;
            for i = 1:length(g_grind.statevars.vectnames)
                if g_grind.statevars.dims{i}.dim1 * g_grind.statevars.dims{i}.dim2 > 1
                    viewcells('-out',i, g_grind.statevars.vectnames{i}, [g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2]);
                end
            end
        end
        if nargin ~= 0
            return;
        end
    end
    if (nargin == 1) && strncmpi(res.opts(i), '-l', 2)
        displayout;
        return;
    elseif (nargin == 1) && strncmpi(res.opts(i), '-p', 2)
        updateparanalreplay(i_grindcolormap);
        return;
        
    end
end
if nargin < 2
    colmap = 'i_grindcolormap';
end
if ~g_grind.statevars.vector
    i_errordlg('The current model has no matrix state variables');
    error('GRIND:viewcells:NoMatrixVar','The current model has no matrix state variables');
end

if res.rerun
    %     disp('running');
    i_ru(t, res.ndays, res.N0, 1);
end
if ~res.silent
    hfig = makevarcontour(colmap);
    addreplaydata(hfig)
    replayall('variable','t');
end
if res.adding
    g_grind.solver.addmode = false;
end
if ~isempty(res.OldY)
    g_Y = res.OldY;
    g_t = res.Oldt;
end
function displayout
global g_grind;
if ~isfield(g_grind, 'viewcells')
    viewcells('-d');
end
for no = 1:length(g_grind.viewcells.vars)
    if ~isempty(g_grind.viewcells.vars{no}.name)
        fprintf('viewcells -out %d %s [%d %d]\n', no, g_grind.viewcells.vars{no}.name, ...
            g_grind.viewcells.vars{no}.dims);
    end
end
function h = viewcellsplot(hfig,data, ylim,colmap)
h = imagesc(data);
hax = get(h, 'parent');
%h=surf(data);
daspect([1 1 1]);
xlabel('Column');
ylabel('Row');
if ~isoctave&&verLessThan('matlab','8.4.0')
    set(hax, 'drawmode','fast','ydir','reverse');
else
    set(hax,'SortMethod','depth','ydir','reverse');
end
set(h,'CDataMapping','scaled');%,'BackFaceLighting','unlit');erasemode normal deleted
shading flat;
colormap(colmap);
h1 = colorbar;
ylim=ylim(~isnan(ylim)&~isinf(ylim));
if ~isempty(ylim)
    ylim = sort(ylim);
    if (length(ylim)==1)||(ylim(2) == ylim(1))
        ylim(2) = ylim(1) + 1;
    end
    set(hax,'CLimMode','manual','CLim',ylim);
    set(h1, 'ylim', ylim);
end
i_plotdefaults(hfig);
function hfig = makevarcontour(colmap)
global g_grind;
statevars = g_grind.viewcells.vars;
hfig = zeros(length(statevars), 1);
for i = 1:length(statevars)
    hfig(i) = i_makefig('varcontour', i);
    set(hfig(i),'doublebuffer','on')
    f1 =  find(strcmp(statevars{i}.name, g_grind.statevars.vectnames), 1);
    if ~isempty(f1)
        statevars{i}.dims = [g_grind.statevars.dims{f1}.dim1, g_grind.statevars.dims{f1}.dim2];
    end
    set(hfig(i), 'Name', ['Grid of ', statevars{i}.name]);
    yy  = outfun(statevars{i}.name);
    %  g_grind.viewcells.y{i} = outfun(statevars{i}.name);
    if any(~isreal(yy(:)))
        yy=real(yy);
        warning('grind:viewcells:complex','Ignored imaginary parts of complex results') 
    end
    maxY = max(max(yy));
    minY =  min(min(yy));
    lasti = size(yy, 1);
    A1 = reshape(yy(lasti,:), statevars{i}.dims(1),statevars{i}.dims(2));
    viewcellsplot(hfig(i),A1, [minY maxY], colmap);
end
function addreplaydata(hfig)
global g_t;
for i = 1:length(hfig)
    hax=findobj(hfig(i),'type','axes');
    tags = get(hax, 'tag');
    hax=hax(~(strcmp(tags,'legend')|strcmp(tags,'Colorbar')));
    ud = get(hax, 'userdata');
    ud.replay.callback = @replaycallback;
    ud.replay.onstart = @onstart;
    ud.replay.onend = [];
    ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
    ud.replay.no  = i;
    set(hax, 'userdata', ud);
end
function onstart(hax)
global  g_t;
if ishandle(hax)
    %    N0 = i_initvar;
    %    if i_settingschanged(N0, g_grind.ndays)
    %       i_ru(t, g_grind.ndays, N0, 1);
    %    end
    ud = get(hax, 'userdata');
    if ~isempty(g_t)
        ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
    end
    set(hax, 'userdata', ud);
    i_figure(get(hax, 'parent'));
end
function t = replaycallback(hax,avar, relt)
global g_grind;
t = [];
if ishandle(hax)&&isfield(g_grind,'viewcells')&&isempty(avar)||strcmp(avar,'t')
    ud = get(hax, 'userdata');
    t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
    ser = get(hax, 'children');
    if isfield(ud, 'replay')
        i = ud.replay.no;
        A1 = real(reshape(outfun(g_grind.viewcells.vars{i}.name, 'times', t),g_grind.viewcells.vars{i}.dims));
        if g_grind.viewcells.vars{i}.dims(2) == 1
            A1 = repmat(A1, 1, 3);
        end
        if g_grind.viewcells.vars{i}.dims(1) == 1
            A1 = repmat(A1, 3, 1);
        end
        if length(ser) == 1
            set(ser, 'CData', A1);
            if isprop(ser, 'ZData')
                set(ser, 'ZData', A1);
            end
        else
            set(ser{1}, 'CData', A1);
            if isprop(ser, 'ZData')
                set(ser(1), 'ZData', A1);
            end
        end
    end
end
%*******************Replay paranal**************************
%
%
function updateparanalreplay(colmap)
global g_paranal g_grind;
if isempty(g_paranal)||isempty(g_paranal.run)
    error('grind:null:paranal','Paranal shouldbe run before using this option');
end
if ~isfield(g_grind, 'viewcells')
    viewcells('-d');
end
for i = 1:length(g_grind.viewcells.vars)
    if ~isempty(g_grind.viewcells.vars{i})
        hfig = i_makefig('varcontour', i);
        if ishandle(hfig)
            hax=findobj(hfig,'type','axes');
            tags = get(hax, 'tag');
            if ~isempty(tags)&&any(~(strcmp(tags,'legend')|strcmp(tags,'Colorbar')))
                hax=hax(~(strcmp(tags,'legend')|strcmp(tags,'Colorbar')));
            end
            ud = get(hax, 'userdata');
        end
        set(hfig,'doublebuffer','on')
        set(hfig, 'Name', ['Grid of ', g_grind.viewcells.vars{i}.name]);
        Ys = outfun(g_grind.viewcells.vars{i}.name, '-p');
        maxY = max(max(Ys));
        minY =  min(min(Ys));
        lasti = size(Ys, 1);
        A1 = reshape(Ys(lasti,:), g_grind.viewcells.vars{i}.dims(1),g_grind.viewcells.vars{i}.dims(2));
        viewcellsplot(hfig,A1, [minY maxY], colmap);
        ud.replay.callback = @paranalreplaycallback;
        ud.replay.onstart = @paranalreplaystart;
        ud.replay.onend = @paranalreplayend;
        ud.replay.onturn = @i_replayparanalturn;
        ud.replay.pars = unique(g_paranal.run.parvalues);
        tstep = ud.replay.pars(2) - ud.replay.pars(1);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,:) g_paranal.run.parvalues(end)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',size(g_paranal.run.parvalues,1));
        ud.replay.no = i;
        hax=findobj(hfig,'type','axes');
        tags = get(hax, 'tag');
        hax=hax(~(strcmp(tags,'legend')|strcmp(tags,'Colorbar')));
        set(hax, 'userdata', ud)
    end
end
replayall('variable', g_paranal.run.pars{1});

function paranalreplaystart(hax)
global g_paranal;
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ud.pars = unique(g_paranal.run.parvalues);
        tstep = ud.replay.pars(2) - ud.replay.pars(1);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,:) g_paranal.run.parvalues(end)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',size(g_paranal.run.parvalues,1));
    end
    set(hax, 'userdata', ud);
end
function paranalreplayend(hax, closedlg)
if closedlg
    paranalreplaycallback(hax,'', 1);
end

function p = paranalreplaycallback(hax,avar, relt)
global g_paranal g_grind;
p = [];
%if ishandle(hax)&&isempty(avar)||strcmp(avar, g_paranal.run.pars{1});
if ishandle(hax)&&isempty(avar)||(~isempty(g_paranal.run)&&strcmp(avar, g_paranal.run.pars{1}))
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ndx = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
        [~, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx);
        p = g_paranal.run.parvalues(stepndx);
        %
        %       ndx = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
        %       ndx=find(g_paranal.p >= ndx, 1);
        if isempty(ndx)
            ndx = numel(g_paranal.run.t);
        end
        cdata = outfun(g_grind.viewcells.vars{ud.replay.no}.name, '-#ndxofparanal', ndx);
        ser = get(hax, 'children');
        set(ser, 'Cdata', cdata);
    end
end

