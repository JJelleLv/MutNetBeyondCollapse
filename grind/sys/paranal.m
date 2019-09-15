%PARANAL   1D parameter analyser
%   Change one parameter step-by-step and show the attractor by simulation. 
%   A 2D or 3D figure is created with the results (on X axis the parameter,
%   on Y axis the parameter/function of the x-axis of phase plane and on 
%   the Z-axis the parameter/function of the y-axis of the phase plane (if any).
%
%   Usage:
%   PARANAL - user is prompted for information
%   PARANAL PREVIOUS - replot the previous analysis (1 = replot, -1 = backwards)
%   PARANAL -1 - run backwards
%   PARANAL 1 - replot the previous analysis
%   PARANAL -out (-o) - change the default output in a dialog box.
%   PARANAL -out plotno [<param1> funy1 funy2 funz1] [minx maxx] [miny maxy] [minz maxz] 
%   - sets the output in a command line: plotno  = number of plot,  is first parameter
%   funy1 = 1st variable or function for y axis [minx maxx] range for xaxis, etc.
%   PARANAL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'lines' [integer>=0] - 0 = scatter plot, 1= line plot, 2 = contour plot (2D only), 3= surf plot (2D only)
%     'nend' [number and length(number)<=2] - The final value(s) of the parameter(s)
%     'outputtype' [integer>0] - The kind of output, name or number
%     'par' [parameter and length(parameter)<=2] - The parameter(s) that need to be changed (max 2)
%     'previous' [integer] - Plot the previous analysis or run backwards (1 = replot, -1 = backwards)
%     'silent' [logical] - Silent run, no figures
%     'stabil' [integer>=0] - Number of time units for stabilizing
%     'start' [number and length(number)<=2] - Start value(s) of the parameter(s)
%     'steps' [integer>0 and length(integer)<=2] - Number of steps for changing the parameter(s)
%     'writing' [integer>=0] - Number of time units for writing the result (0= one value)
%   PARANAL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-d' - default, reset the default output.
%     '-l' - list (-l) - list the outputs for paranal.
%     '-lo' - -lo=filename load the results of a previous paranal session
%     '-s' - -s=filename save the results of the last or current paranal to a mat 
%   file with name "filename.mat".
%
%   See also ax, conteq, paranal2d, forward_stabil, open_matcont_gui
%
%   Reference page in Help browser:
%      <a href="matlab:commands('paranal')">commands paranal</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function paranal(varargin)
global g_paranal g_grind;

if nargin>0&&strncmpi(varargin{1}, '-o', 2)
    if nargin == 1
        i_paranaloutdlg;
        paranal
    else
        p.xaxis = {''};
        p.yaxis = {''};
        p.zaxis = {''};
        p.xlim = [0 10];
        p.ylim = [0 10];
        p.zlim = [0 10];
        ip = 1;
        j = 0;
        for i = 2:nargin
            v = str2num(varargin{i}); 
            if isempty(v)
                %analyse axes
                s = varargin{i};
                if s(1) == '['
                    s = s(2:end - 1);
                end
                f=strfind(s,'''');
                if isempty(f)
                    f = strfind(s, ' ');
                    ff = repmat(f, 2, 1); %this is a trick to get double values
                    f = [0; ff(:); length(s) + 1];
                end
                if length(f) > 1
                    p.xaxis = mystr2cell(s(f(1) + 1:f(2) - 1));
                end
                if length(f) > 3
                    p.yaxis = mystr2cell(s(f(3) + 1:f(4) - 1));
                end
                if length(f) > 5
                    if f(5) < f(6) - 1
                        p.zaxis = mystr2cell(s(f(5) + 1:f(6) - 1));
                    end
                end
            elseif length(v) == 1
                ip = v;
            elseif length(v) == 2
                switch j
                    case 0
                        p.xlim = v;
                    case 1
                        p.ylim = v;
                    case 2
                        p.zlim = v;
                end
                j = j + 1;
            end
        end
        if (isempty(p.zaxis)||isempty(p.zaxis{1}))&&(isempty(p.yaxis)||isempty(p.yaxis{1}))
            error('GRIND:paranal:noYaxis','Cannot plot with only the x axis defined');
        end
        g_grind.paranal.plots{ip} = p;
        g_grind.paranal.defaultplots=0;
        g_grind.paranal.currno=1;
    end
    return;
end
  % 'outputtype', 'i>0#e[unchanged|mean|median|maxima|minima|minima+maxima|sd|cv|autoregr|skewness|sum|perc05|perc95|range90|min|max|range]', 'The kind of output, name or number',0;...
try
    apar=g_grind.pars{1};
catch
    apar='';
end
fieldnams={'previous', 'i', 'Plot the previous analysis or run backwards (1 = replot, -1 = backwards)',[];...
   'par', 'p&length(p)<=2', 'The parameter(s) that need to be changed (max 2)',apar;...
   'start', 'n&length(n)<=2', 'Start value(s) of the parameter(s)',0;...
   'nend', 'n&length(n)<=2', 'The final value(s) of the parameter(s)',10;...
   'steps', 'i>0&length(i)<=2', 'Number of steps for changing the parameter(s)',50;...
   'stabil', 'i>=0', 'Number of time units for stabilizing',600;...
   'silent', 'l', 'Silent run, no figures',false;...
   'writing', 'i>=0', 'Number of time units for writing the result (0= one value)',300;...
   'outputtype', 'i>0', 'The kind of output, name or number',1;...
   'lines', 'i>=0', '0 = scatter plot, 1= line plot, 2 = contour plot (2D only), 3= surf plot (2D only)',0}';
args=i_parseargs(fieldnams,'previous','-lo,-l,-d,-s',varargin);
if (nargin==0)|| ~ischar(varargin{1}) || ~strncmpi(varargin{1}, '-lo', 3) 
    i_parcheck;
end
if ~isfield(args,'par')&&~isfield(args,'previous')
   args.par={'',''};
end
if any(strcmp(args.opts,'-d'))
    g_grind.paranal.plots = {};
    avar=g_grind.xaxis.var;
    iX=i_getno(avar);
    if ~iX.isvar
        avar=i_statevars_names(1);
    end
    avar2=g_grind.yaxis.var;
    iY=i_getno(avar2);
    if ~iY.isvar||strcmp(avar,avar2)
        avar2=i_statevars_names(2);
    end
    g_grind.paranal.plots{1}.xaxis = {'<param1>'};
    g_grind.paranal.plots{1}.xlim = [0 10];
    g_grind.paranal.plots{1}.yaxis = {avar};
    g_grind.paranal.plots{1}.ylim = g_grind.xaxis.lim;
    g_grind.paranal.plots{1}.zaxis = {avar2};
    g_grind.paranal.plots{1}.zlim = g_grind.yaxis.lim;
    g_grind.paranal.plots2d = {};
    g_grind.paranal.plots2d{1}.xaxis = {'<param1>'};
    g_grind.paranal.plots2d{1}.xlim = [0 10];
    g_grind.paranal.plots2d{1}.yaxis = {'<param2>'};
    g_grind.paranal.plots2d{1}.ylim = [0 10];
    g_grind.paranal.plots2d{1}.zaxis = {avar};
    g_grind.paranal.plots2d{1}.zlim = g_grind.xaxis.lim;
    g_grind.paranal.defaultplots = 1;
    g_grind.paranal.currno = 1;
    return;
end
outfile = '';
if isfield(g_grind, 'pars')&&isempty(g_grind.pars)
    error('GRIND:paranal:NoPars','No parameters to analyse');
end
olddlg = [];
settingschanged=~(isfield(g_grind,'paranal')&&isfield(g_paranal,'settings'))||~struccmp(par('-v',0),g_paranal.settings);
if ~settingschanged && isfield(g_grind.paranal, 'dlg')
    olddlg = g_grind.paranal.dlg;
end
if ~isfield(args,'previous')
  args.previous = 0;
end

f1=strncmpi(args.opts, '-s',2);
if any(f1)
    outfile = 'paranal.mat';
    arg=args.opts{f1};
    f=strfind(arg, '=');
    if ~isempty(f)
        outfile = arg(f(1) + 1:end);
    else
        outfile =  fullfile(grindpath, outfile);
    end
    if ~isempty(g_paranal.run) && (nargin==1)
        saveparanal(outfile);
        return;
    end
end
f1=strncmpi(args.opts, '-lo',3);
if any(f1)
    arg=args.opts{f1};
    infile = 'paranal.mat';
    f=strfind(arg, '=');
    if ~isempty(f)
        infile = arg(f(1) + 1:end);
    else
        infile =  fullfile(grindpath, infile);
    end
    if ~strcontains(infile, '.')
        infile = [infile '.mat']; %(faster than sprintf)
    end
    if ~exist(infile, 'file')
        [infile,path]=uigetfile('*.mat','File for paranal results');
        if isempty(infile)
            disp('File not found');
            return;
        end
        cd(path);
    end
    loadparanal(infile);
    if nargin == 1
        paranal(1);
        return;
    end
end
if any(strcmp(args.opts, '-l'))
    if ~isfield(g_grind, 'paranal')
        paranal('-default');
    end
    for i = 1:length(g_grind.paranal.plots)
        fprintf('paranal -out %d [''%s'' ''%s'' ''%s''] %s %s %s\n',i, strtrim(sprintf('%s ',g_grind.paranal.plots{i}.xaxis{:})),...
            strtrim(sprintf('%s ',g_grind.paranal.plots{i}.yaxis{:})),strtrim(sprintf('%s ',g_grind.paranal.plots{i}.zaxis{:})),...
            mat2str(g_grind.paranal.plots{i}.xlim), mat2str(g_grind.paranal.plots{i}.ylim), mat2str(g_grind.paranal.plots{i}.zlim));
    end
    return;
end
%doreplay = replayall('close',1);
if any(args.previous==[1 -1 2])&&isempty(g_paranal.run)
    error('grind:paranal','There is no a previous analysis, first run paranal');
end
if args.previous == 2
    if ~isempty(g_paranal.prevrun)
        h=g_paranal.run;
        g_paranal.run=g_paranal.prevrun;
        g_paranal.prevrun=h;
        paranal(1);
        h=g_paranal.run;
        g_paranal.run=g_paranal.prevrun;
        g_paranal.prevrun=h;
    else
        warning('grind:paranal','There is no a previous analysis, plotting only the last run');
    end
    paranal(1);
    return;
end
answer=args;
if isfield(g_grind,'paranal')&&isfield(g_grind.paranal,'dlg')
    answer = mergestructs(g_grind.paranal.dlg,answer);
end
if isempty(args.previous)||(args.previous ==1) && isempty(answer)
    p.par = g_paranal.run.pars(1);
    p.start = min(g_paranal.run.parvalues(:,1));
    p.nend = max(g_paranal.run.parvalues(:,1));
    answer = i_paranaldialog('initstruct', p);
else
    answer = i_paranaldialog('initstruct', answer);
end
if isfield(args,'silent')
    g_grind.paranal.silent=args.silent;
end
if isempty(answer.par{1})   
    if isfield(g_grind,'paranal')&&isfield(g_grind.paranal,'dlg')
        answer.par=g_grind.paranal.dlg.par;
    end
    if args.previous == 0
        a1 = i_paranaldialog(answer);
        if isempty(a1)
            return
        elseif ~isempty(answer) && isfield(g_paranal.run, 'Y') && ~isempty(g_paranal.run.Y) && ~isempty(g_paranal.run.parvalues) && ...
                all(a1.start==answer.start) && strcmp(a1.par{1}, answer.par{1}) ...
                && strcmp(a1.par{2}, answer.par{2}) && all(a1.steps==answer.steps) ...
                && (a1.stabil==answer.stabil) && (a1.writing==answer.writing) && ...
                prod(double(a1.nend==answer.nend)) && (isfield(g_paranal,'settings')&&struccmp(g_paranal.settings, par('-v', 0))) && ...
                strcmp(questdlg('Do you want to use data of the previous paranal run?','paranal',...
                'Yes','No','No'),'Yes')
            args.previous = 1;
        end
        answer = a1;
    end
end
if ~isempty(answer)
    if ~isfield(g_grind, 'paranal')
        g_grind.paranal.defaultplots = 1;
    end
    answer.par{1} = checkmat(answer.par{1});
    answer.par{2} = checkmat(answer.par{2});
    if args.previous == -1
        s = answer.start;
        answer.start = answer.nend;
        answer.nend = s;
    end
    g_grind.paranal.dlg = answer;
    g_paranal.settings = par('-v', 0);
    if isfield(g_paranal,'sens')
        g_paranal=rmfield(g_paranal,'sens');
    end
    %    if ~isempty(answer.par{2})
    %       paranal2d(answer);
    %       clear answer;
    %       return;
    %    end
    if answer.lines == 2&&isempty(answer.par{2})
        i_warningdlg('GRIND:paranal:contour1D','Contour plot not supported in 1D paranal, making scatterplot instead');
    end
    if answer.lines == 3&&isempty(answer.par{2})
        i_warningdlg('GRIND:paranal:surface1D','Surface plot not supported in 1D paranal, making scatterplot instead');
    end
    if ~(args.previous==1)&&~(settingschanged||paranaldlgchanged(olddlg, g_grind.paranal.dlg))&&~isempty(g_paranal.run)&&(isempty(g_paranal.run.parvalues)||(g_paranal.run.parvalues(1,1)~=g_grind.paranal.dlg.start(1)))
        g_paranal.prevrun = g_paranal.run;
    else
        if settingschanged
            g_paranal.prevrun = [];
            if isfield(g_paranal,'forward_stabil')
                g_paranal=rmfield(g_paranal,'forward_stabil');
            end
            g_paranal.nulldata = [];
            if isfield(g_paranal,'potdata')
                g_paranal.potdata = [];
            end
        end
        if isempty(g_paranal.run)||(g_paranal.run.parvalues(1,1)~=g_grind.paranal.dlg.start(1))
            g_paranal.nulldata = [];
            if isfield(g_paranal,'forward_stabil')
                g_paranal=rmfield(g_paranal,'forward_stabil');
            end
            if isfield(g_paranal,'potdata')
                g_paranal.potdata = [];
            end
        end
    end
    %   function i_paranal(parname, start, nend, nsteps, nstabilizing, ndays, lines, args.previous, outputtype, disturb);
    i_paranal(args.previous);
    if ~isempty(outfile)
        saveparanal(outfile);
    end
    clear answer;
end
% if doreplay
%     replayall ('variable', g_paranal.run.pars{1});
% end
function saveparanal(outfile)
global g_paranal g_grind;  
path = pwd;              
inifile = g_grind.inifile; 
dlg = g_grind.paranal.dlg; 
save(outfile,'g_paranal','path','inifile','dlg');
fprintf('Saved paranal results to "%s"\n', outfile)
function loadparanal(infile)
global g_paranal g_grind; 
path = '';
inifile = '';
dlg = [];
load(infile,'g_paranal','path','inifile','dlg');
if isempty(path)
    i_errordlg('Not a valid file for paranal');
    error('GRIND:paranal:NoFile','Paranal: not a valid file');
end
if ~strcmp(path, pwd)
    cd(path);
end
if isempty(g_grind)||(~isfield('inifile',g_grind)&&~strcmp(inifile, g_grind.inifile))
    use(inifile);
    load(infile,'g_paranal','path','inifile','dlg');
end
g_grind.paranal.dlg = dlg;
function p = checkmat(apar)
global g_grind;
if ~isempty(apar) && ~strcontains(apar, '(')
    oldpar = '';
    s=evalin('base',sprintf('size(%s)',apar));
    if prod(s) > 1
        answer{1} = apar;
        if isfield(g_grind,'paranal')&&isfield(g_grind.paranal,'dlg')
            oldpar = g_grind.paranal.dlg.par{1};
            if ~strncmp(apar, oldpar, length(apar))
                oldpar = '';
            end
        end
        if ~isempty(oldpar)
            answer=inputdlg('Enter which element of parameter you want to use',sprintf('Parameter is %dx%d matrix',s),1,{oldpar});
        else
            if (max(s) > 1) && (min(s) == 1)
                answer=inputdlg('Enter which element of parameter you want to use',sprintf('Parameter is %dx1 vector',max(s)),1,{[apar '(1)']});
            end
            if (max(s) > 1) && (min(s) > 1)
                answer=inputdlg('Enter which element of parameter you want to use',sprintf('Parameter is %dx%d matrix',s),1,{[apar '(1,1)']});
            end
        end
        p = answer{1};
    else
        p = apar;
    end
else
    p = apar;
end
function A = mystr2cell(s)
s = outf('changeshortcut', s);
A=str2cell(strrep(s,' ',sprintf('\n'))); %#ok<SPRINTFN>
function res = paranaldlgchanged(dlg1, dlg2)
res=isempty(dlg1)||isempty(dlg2)||~strcmp(dlg1.par{1}, dlg2.par{1})||dlg1.steps(1)~=dlg2.steps(1);
if ~res
    min1 = min(dlg1.start(1), dlg1.nend(1));
    min2 = min(dlg2.start(1), dlg2.nend(1));
    max1 = max(dlg1.start(1), dlg1.nend(1));
    max2 = max(dlg2.start(1), dlg2.nend(1));
    res=min1~=min2||max1~=max2;
end




