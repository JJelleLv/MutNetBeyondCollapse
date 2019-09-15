%MCARLO   Monte Carlo uncertainty/sensitivity analysis
%   This function iterates a certain run, while changing the parameters and or initial conditions. 
%   Each parameter/initial condition can be selected separately and for each a distribution
%   for drawing each parameter independently can be selected. You can select either a Uniform distribution 
%   (=rand), a normal distribution (=randn), a truncated normal distribution (negative number are discarded) 
%   or a log-normal distribution. You can also use this function to run a stochastic model several 
%   times without changing parameters.
%   The data (in g_mcarlo) are saved to a file, and after the run the data can be analysed. If you select
%   an uncertainty analysis, time plots are generated with ranges (5%, 25% 50% 75% and 95% percentiles).
%   The sensitivity analysis performs a multivariate analysis of the parameter sensitivity (Klepper, 1989).
%   This analysis does not only find the overall effects of each parameter or initial condition, but it
%   also performs a cluster analysis, finding clusters of parameters with a similar effect on the results.
%   For this cluster analysis the statistics toolbox is <a href="matlab:commands toolboxes">required</a> (not for the sensitivity matrix).
%  
%  
%   Klepper, O. (1989) A model of carbon flows in relation to macrobenthic food supply in the Oosterschelde 
%   estuary (S.W. Netherlands), LUW Wageningen, The Netherlands. PhD thesis.
%  
%  
%   Usage:
%   MCARLO - all input in dialog boxes.
%   MCARLO UNCERTAIN - uncertainty analysis, all other input by user interface.
%   MCARLO g_mcarlo - All input is in this structure. Note that if g_mcarlo.i<g_mcarlo.n the previous run
%   is continued, else the previous data are used.
%   MCARLO UNCERTAIN ALLPARS N - ALLPARS is a structure with the parameter settings, N is the number of runs.
%   is continued, else the previous data are used.
%   MCARLO SENS -load=mcarlo - sensitivity analysis while loading the data from mcarlo.mat
%   MCARLO('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'allpars' [general] - cell with structures for each parameter (see -o setting)
%     'dendrpars' [cell] - list of parameters/variables that are used for the dendrogram.
%     'i' [integer>=0] - used internally (the number of runs done sofar)
%     'inifile' [string] - name of the ini-file (is assigned automatically)
%     'n' [integer>=0] - number of runs
%     'ndays' [number>0] - number of time units per run
%     'out' [cell] - list of variables/functions for output
%     'outfile' [string] - save the results to a mat file
%     'paranalpar' [logical] - parameters that are used for paranal
%     'run' [struct] - used internally for the results of the run
%     'silent' [logical] - silent run without output to the screen
%     'stabil' [integer>=0] - number of time units for stabilization
%     'tstep' [integer>0] - number of time units for output
%     'type' [uncertain | sens | ''] - type of run (unceartainty, sensitivity or unknown)
%   MCARLO('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - information about the status of MCARLO
%     '-l=mcarlo' - sensitivity analysis while loading the data from mcarlo.mat
%     '-o' - get a default structure for each parameter. Fields: name: name of parameter; descr: 'Parameter'
%     value: current value; selected: selected or not; range: relative range; min: minimum of range; max: maximum
%     of range; sd: standard deviation; distr: 'Uniform', 'Normal', 'TruncNormal','LogNormal', 'paranal','hysteresis'
%
%
%   See also paranal, out, time
%
%   Reference page in Help browser:
%      <a href="matlab:commands('mcarlo')">commands mcarlo</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res=mcarlo(varargin)
%(atype, opt)
evalin('base','global g_mcarlo;');
global g_mcarlo g_grind g_t g_Y t g_permanent g_paranal;
fieldnams={'type', 'e[u+ncertain|s+ens|]', 'type of run (unceartainty, sensitivity or unknown)','sens';...
   'allpars', '', 'cell with structures for each parameter (see -o setting)','';...
   'stabil', 'i>=0', 'number of time units for stabilization',0;...
   'ndays', 'n>0', 'number of time units per run',g_grind.ndays;...
   'tstep', 'i>0', 'number of time units for output',10;...
   'n', 'i>=0', 'number of runs',1000;...
   'i', 'i>=0', 'used internally (the number of runs done sofar)',0;...
   'out', 'c', 'list of variables/functions for output',{};...
   'inifile', 's', 'name of the ini-file (is assigned automatically)','';...
   'paranalpar', 'l', 'parameters that are used for paranal',[];...
   'run', 'r', 'used internally for the results of the run',[];...
   'outfile', 's', 'save the results to a mat file','';...
   'dendrpars', 'c', 'list of parameters/variables that are used for the dendrogram.','';...
   'silent', 'l', 'silent run without output to the screen',false}';
args=i_parseargs(fieldnams,'type,allpars,n','-?,-l,-o',varargin);
if any(strcmp(args.opts,'-o'))
   res= i_mcarlodlg('initstruc');
   return;
end
g_mcarlo=mergestructs(g_mcarlo,args);
%g_mcarlo = updatestruc(g_mcarlo);
if ~isfield(g_mcarlo,'silent')
    g_mcarlo.silent=false;
    f=find(strcmpi(g_grind.figopts,'visible'));
    if ~isempty(f)
        g_mcarlo.silent=strcmpi(g_grind.figopts{f+1},'off');
    end
end
if ~isfield(args,'type')
    g_mcarlo.type='';
else
    g_mcarlo.type = upper(args.type);
end
g_mcarlo.inifile = g_grind.inifile;
if ~isfield(g_mcarlo,'ndays')
    g_mcarlo.ndays=g_grind.ndays;
end
f1=strncmpi(args.opts,'-l',2);
if any(f1)
    %load from file and continue if necessary
    opt=args.opts{f1};
    f=strfind(opt, '=');
    if ~isempty(f)
        load(opt(f(1) + 1:end),'g_mcarlo');
        args=mergestructs(g_mcarlo,args);
        g_mcarlo=args;        
        if ~isfield(g_mcarlo,'silent')
           g_mcarlo.silent=false;
        end
%        g_mcarlo = updatestruc(g_mcarlo);
        if g_mcarlo.n == g_mcarlo.i
            disp('previous results loaded');
            if ~isempty(g_grind.permanent)
                time(1, '-s');
            end
        end
    end
    if ~strcmpi(g_mcarlo.inifile, g_grind.inifile)
        warning('GRIND:mcarlo:inifile','Data generated with another inifile (%s) than currently opened',g_mcarlo.inifile);
    end
else
    g_mcarlo.i=0;
end
if any(strcmp(args.opts, '-?'))
    fprintf('inifile: %s\n', g_mcarlo.inifile);
    fprintf('outfile: %s\n', g_mcarlo.outfile);
    disp('Parameters:')
    for i = 1:length(g_mcarlo.allpars)
        if g_mcarlo.allpars(i).selected
            if strcmp(g_mcarlo.allpars(i).distr, 'Uniform')
                fprintf('%s=%g;  %%selected: %s: [%g-%g]\n', g_mcarlo.allpars(i).name, g_mcarlo.allpars(i).value, ...
                    g_mcarlo.allpars(i).distr, g_mcarlo.allpars(i).min, g_mcarlo.allpars(i).max);
            elseif strcmp(g_mcarlo.allpars(i).distr, 'paranal')
                fprintf('%s=%g;  %%paranal: [%g-%g] %d steps\n', g_mcarlo.allpars(i).name, g_mcarlo.allpars(i).value, ...
                    g_mcarlo.allpars(i).minmax(1), g_mcarlo.allpars(i).minmax(2), g_mcarlo.allpars(i).steps);
            else
                fprintf('%s=%g;  %%selected: %s: sd=%g\n', g_mcarlo.allpars(i).name, g_mcarlo.allpars(i).value, ...
                    g_mcarlo.allpars(i).distr, g_mcarlo.allpars(i).sd);
            end
        else
            fprintf('%s=%g;\n', g_mcarlo.allpars(i).name, g_mcarlo.allpars(i).value);
        end
    end
    disp(' ');
    return;
end
i_parcheck;
%ndays = g_grind.ndays;
oldnsteps =  g_grind.tstep;
hasperm = ~isempty(g_grind.permanent);
if ~isfield(args,'allpars')||isempty(args.allpars)
    i_mcarlodlg;
    g_mcarlo.i = 0;
    uiwait;
end
nparanal = 0;
g_mcarlo.paranalpar  = [];
for j = 1:length(g_mcarlo.allpars)
    if g_mcarlo.allpars(j).selected(1) && any(strcmpi({'paranal','hysteresis'},g_mcarlo.allpars(j).distr))
        nparanal = nparanal + 1;
        g_mcarlo.paranalpar = j;
    end
end
if nparanal > 1
    i_errordlg(sprintf('There are now %d parameters for paranal/hysteresis selected, select one maximal.',nparanal),'mcarlo')
end
if ~isfield(args, 'n')
    if ~isempty(g_mcarlo.paranalpar)
        Ans=inputdlg({'Number of iterations (paranal):', 'File for results'}, 'mcarlo',1 ,{'1000','mcarlo.mat'});
        if ~isempty(Ans)
            p = g_mcarlo.allpars(g_mcarlo.paranalpar);
            g_mcarlo.paranal.par = {p.name  ''};
            g_mcarlo.paranal.start = [p.minmax(1) 0];
            g_mcarlo.paranal.nend = [p.minmax(2) 10];
            g_mcarlo.paranal.steps = [p.steps 50];
            g_mcarlo.paranal.stabil = p.nstabilwrite(1);
            g_mcarlo.paranal.writing = p.nstabilwrite(2);
            g_mcarlo.paranal.hysteresis = strcmpi(p.distr, 'hysteresis');
            if g_mcarlo.paranal.writing(1) == 0
                nsteps = g_mcarlo.paranal.steps(1) * (g_mcarlo.paranal.writing(1) + 1);
            else
                nsteps = (g_mcarlo.paranal.steps(1) + 1) * (g_mcarlo.paranal.writing(1) + 1) - 1;
            end
            Ans = [{'0'; num2str(ndays); num2str(nsteps)}; Ans];
        end
    else
        if ~isfield(g_mcarlo,'n')
            g_mcarlo.n=1000;
        end
        if ~isfield(g_mcarlo,'stabil')
            g_mcarlo.stabil=0;
        end
        if ~isfield(g_mcarlo,'tstep')
            g_mcarlo.tstep=10;
        end
        if ~isfield(g_mcarlo,'outfile')
            g_mcarlo.outfile='mcarlo.mat';
        end
        Ans=inputdlg({'Number of time units to stabilize before each run','Number of time units to simulate (each iteration):', 'Number of values of outcomes:',...
            'Number of iterations:', 'File for results'}, 'mcarlo',1 ,{int2str(g_mcarlo.stabil),int2str(g_mcarlo.ndays),...
            int2str(g_mcarlo.tstep),int2str(g_mcarlo.n),g_mcarlo.outfile});
    end
    if isempty(Ans)
        disp('Cancel pressed')
        return;
    end
    g_mcarlo.stabil = str2double(Ans{1});
    g_mcarlo.ndays = str2double(Ans{2});
    g_mcarlo.tstep = str2double(Ans{3});
    g_mcarlo.n = str2double(Ans{4});
    g_mcarlo.outfile = Ans{5};
    g_mcarlo.inifile = g_grind.inifile;
end
if isempty(g_grind.tstep)
    error('GRIND:mcarlo:notfixed','None of the outcomes should be fixed');
end

if ~isfield(g_mcarlo,'tstep')
    g_mcarlo.tstep=g_grind.tstep;
    if isnan(g_mcarlo.tstep)
        g_mcarlo.tstep=10;
    end
end
g_mcarlo.run.pars = cell(size(length(g_mcarlo.allpars)));
g_mcarlo.run.parndx = zeros(length(g_mcarlo.allpars), 3);
k = 1;
nnegl = 0;
for j = 1:length(g_mcarlo.allpars)
    if any(g_mcarlo.allpars(j).selected(:)) && ~any(strcmpi({'paranal','hysteresis'},g_mcarlo.allpars(j).distr))
        for m = 1:size(g_mcarlo.allpars(j).value, 1)
            for n = 1:size(g_mcarlo.allpars(j).value, 2)
                if g_mcarlo.allpars(j).selected(m, n)
                    pname = g_mcarlo.allpars(j).name;
                    if  (size(g_mcarlo.allpars(j).value, 2) > 1)
                        pname = sprintf('%s(%d,%d)',pname,m,n);
                    elseif (numel(g_mcarlo.allpars(j).value) > 1)
                        pname = sprintf('%s(%d)', pname, m);
                    end
                    psd = g_mcarlo.allpars(j).sd(m, n);
                    if ~isnan(psd) && (psd~=0)
                        g_mcarlo.run.pars{k} = pname;
                        g_mcarlo.run.parndx(k, :) = [j, m, n];
                        k = k + 1;
                    else
                        nnegl = nnegl + 1;
                    end
                end
            end
        end
    end
end
if k <= size(g_mcarlo.run.parndx, 1)
    g_mcarlo.run.parndx =  g_mcarlo.run.parndx(1:k - 1, :);
    g_mcarlo.run.pars = g_mcarlo.run.pars(1:k - 1);
end

if nnegl > 0
    warning('GRIND:mcarlo:varsneglected','%d parameters/variables were neglected as their range was zero.\n', nnegl);
end

npar = length(g_mcarlo.run.pars);
nsteps=g_mcarlo.n;
g_grind.tstep=g_mcarlo.tstep;
if ~isfield(g_mcarlo, 'i')||(g_mcarlo.i < g_mcarlo.n)
    g_mcarlo.i = 0;
    g_mcarlo.run.parvalues = zeros(nsteps, npar);
    if g_grind.tstep <= 1
        g_mcarlo.run.Y = zeros(1, g_grind.statevars.dim, nsteps);
    else
        g_mcarlo.run.Y = zeros(g_grind.tstep + 1, g_grind.statevars.dim, nsteps);
    end
    if hasperm
        time(1, '-s')
        g_mcarlo.run.perm = zeros(size(g_mcarlo.run.Y,1), g_grind.permanent{end}.to,nsteps);
    else
        g_mcarlo.run.perm = [];
    end
end

g_grind.solver.opt.OutputFcn = [];
try
    %Here starts the running
    for j=1:size(g_mcarlo.allpars)
        if any(g_mcarlo.allpars(j).selected)
            g_mcarlo.allpars(j).value=evalin('base',g_mcarlo.allpars(j).name);
        end
    end
    runned = 0;
    if ~isempty(g_mcarlo.paranalpar)
        g_grind.tstep =  max(1, g_mcarlo.paranal.writing);
    end
    i_waitbar(0, g_mcarlo.n, 'mcarlo','Running',2)
    for i = g_mcarlo.i + 1:g_mcarlo.n
        runned = 1;
        drawnow;
        % setting the parameter values
        if g_grind.statevars.vector
            k=1;
            for j=1:size(g_mcarlo.allpars)
                if any(g_mcarlo.allpars(j).selected(:))
                    ppar = getparvalue(g_mcarlo.allpars(j), size(g_mcarlo.allpars(j).value));
                    valu=g_mcarlo.allpars(j).value;
                    valu(g_mcarlo.allpars(j).selected)=ppar(g_mcarlo.allpars(j).selected);
                    %assign per parameter for speed
                    multassignin('base', g_mcarlo.allpars(j).name, valu);
                    ppar= ppar(g_mcarlo.allpars(j).selected);
                    g_mcarlo.run.parvalues(i, k:k+numel(ppar)-1) = transpose(ppar);
                    k=k+numel(ppar);
                end
            end
        else
            %can be very slow for vector models
            for j = 1:size(g_mcarlo.run.parndx,1)
                ppar = getparvalue(g_mcarlo.allpars(g_mcarlo.run.parndx(j, 1)), [1,1]);
                multassignin('base', char(g_mcarlo.run.pars{j}), ppar);
                g_mcarlo.run.parvalues(i, j) = ppar;
            end
        end
        %can also be state variables be changed
        N0 = i_initvar;
        if ~isempty(g_mcarlo.paranalpar)
            i_keep(N0);
            %paranal (only one parameter allowed)
            if ~g_mcarlo.paranal.hysteresis  %paranal
                era();
                g_grind.paranal.silent = 1;
                paranal(g_mcarlo.paranal)
                if i == 1
                    g_mcarlo.run.t = g_paranal.run.t(:);
                    g_mcarlo.paranal.parvalues = g_paranal.run.parvalues;
                end
                g_mcarlo.run.Y(:, :, i) = transpose(reshape(permute(g_paranal.run.Y, [2, 1, 3]), size(g_paranal.run.Y, 2), size(g_paranal.run.Y, 1) * size(g_paranal.run.Y, 3)));
                if hasperm
                    g_mcarlo.run.perm(:, :, i) =  g_paranal.perm;
                end
            else  %hysteresis - back and forth run of paranal
                era();
                paranal(g_mcarlo.paranal);
                paranal('-1')
                if i == 1
                    g_mcarlo.run.t = cat(3,g_paranal.prevrun.t, g_paranal.run.t);
                    g_mcarlo.paranal.parvalues = [g_paranal.prevrun.parvalues; g_paranal.run.parvalues];
                end
                g_mcarlo.run.Y(:, :, i) = cat(3,g_paranal.prevrun.Y, g_paranal.run.Y);
                if hasperm
                    g_mcarlo.run.perm(:, :, i) = cat(3,g_paranal.prevrun.perm, g_paranal.run.perm);
                end
            end
        else
            %running a normal model run
            if g_mcarlo.stabil > 0 %optionally stabilize
                i_keep(N0);
                stabil(g_mcarlo.stabil, '-s');
                N0 = i_initvar;
            end
            if g_grind.tstep <= 1
                stabil(g_mcarlo.ndays, '-s');
                g_Y = g_Y(end, :);
                g_t = g_t(end);
            else
                i_ru(t, g_mcarlo.ndays, N0, 0);
            end
            if size(g_Y,1)<g_mcarlo.tstep %in case there is an error and not the whole data is filled
                %complement with nan
                n=size(g_Y,1);
                g_Y(end+1:g_mcarlo.tstep+1,:)=nan(g_mcarlo.tstep-n+1,size(g_Y,2));
                g_t(end+1:g_mcarlo.tstep+1)=nan(g_mcarlo.tstep-n+1,1);
            end
            if i == 1
                g_mcarlo.run.t = g_t;
            end
            g_mcarlo.run.Y(:, :, i) = g_Y;
            if hasperm
                pperm = defpermanent('-g', []);
                g_mcarlo.run.perm(:, :, i) =  pperm;
            end
        end
        g_mcarlo.i = i;
        i_waitbar(1);
        if mod(i, 100) == 0
            save(g_mcarlo.outfile, 'g_mcarlo');
        end
    end
    if runned %save only if new data are generated
        save(g_mcarlo.outfile, 'g_mcarlo');
        fprintf('Saved data to %s\n', g_mcarlo.outfile);
    end
    i_waitbar([]);
    %set the parameters/statevars back to their original values
    g_grind.tstep = oldnsteps;
    for j = 1:length(g_mcarlo.allpars)
        if g_mcarlo.allpars(j).selected
            multassignin('base', char(g_mcarlo.allpars(j).name), g_mcarlo.allpars(j).value);
        end
    end
catch err
    %  err=lasterror;
    i_waitbar([]);
    g_grind.tstep = oldnsteps;
    save(g_mcarlo.outfile, 'g_mcarlo');
    for j = 1:length(g_mcarlo.allpars)
        if g_mcarlo.allpars(j).selected
            multassignin('base', char(g_mcarlo.allpars(j).name), g_mcarlo.allpars(j).value);
        end
    end
    rethrow(err);
end

%analyse the results either in an uncertainty analysis or sensitivity
%analysis
if ~isfield(args,'type')||isempty(args.type)
    args.type=questdlg('Select the analysis', 'Monte Carlo'...
        ,'Uncertainty analysis of the results','Sensitivity analysis of parameters','Sensitivity analysis of parameters');
end
if ~isempty(args.type) && (args.type(1) == 'U') %uncertainty analysis
    % the uncertainty analysis plots the timeplots but with the 5% 25% 50% 75% 95% percentiles
    for No = 1:size(g_grind.timevars, 2)
        if ~isempty(g_grind.timevars{No})
            hfig=i_makefig('time', No - 1);
            if g_mcarlo.silent
                set(hfig,'visible','off');
            end
            hold on;
            % plotedit('off');
            apen = g_grind.pen;
            apen.i = 1; %  do not cycle colors between runs
            apen.cycle = true;
            
            if strcmp(g_grind.outt{No}, 't')
                tt = g_mcarlo.run.t; %efficiency
            else
                %repeated calls to outfun are not efficient, maybe needs to
                %be reimplemented
                tt = mean(reshape(outfun(g_grind.outt{No}, '-m'), numel(g_mcarlo.run.t), g_mcarlo.n), 2);
            end
            %use the variables defined in 'out'
            for i = 1:length(g_grind.timevars{No})
                s = g_grind.timevars{No}{i};
                if ~strncmp(s, 'Observed ', 9)
                    xx = reshape(outfun(s, '-m'), numel(g_mcarlo.run.t), g_mcarlo.n);
                    xx =  i_makepercentiles(xx, [0.05, 0.25, 0.5, 0.75, 0.95]);
                    plot(tt, xx(:,[1 5]), apen.pen, 'Color', apen.color,'Linewidth',0.5, 'LineStyle',':');
                    plot(tt, xx(:,[2 4]), apen.pen, 'Color', apen.color,'Linewidth',1);
                    plot(tt, xx(:,3), apen.pen, 'Color', apen.color,'Linewidth',2);
                    apen = i_nextpen(apen);
                end
            end
            xlabel( g_grind.outt{No});
            s = i_disptext(g_grind.timevars{No}{1});
            ylabel(s);
        end
    end
elseif ~isempty(args.type)
    %sensitivity analysis
    %disp('sensitivity analysis not yet implemented')
    % here code to create (1) a sensitivity matrix
    % 1 matrix per kolom (voor alle SI verschillend)
    
    % - define parameters
    Nsteps = numel(g_mcarlo.run.t);
    Niter = g_mcarlo.i;
    Npars = size(g_mcarlo.run.parvalues, 2);
    
    %2 - get data
    %dit is alleen van belang voor vector variabelen,
    %bvzou ook relevante willen selecteren, bijv met out
    if isfield(args, 'out')&&~isempty(args.out)
        answer=args.out;
    else
        if isfield(g_mcarlo, 'out')
            s = strtrim(sprintf('%s\n', g_mcarlo.out{:}));
        elseif g_grind.statevars.vector
            s = strtrim(sprintf('%s\n', g_grind.statevars.vectnames{:}));
        else
            s = strtrim(sprintf('%s\n', g_grind.statevars.names{:}));
        end
        answer=inputdlg({sprintf('Enter variables or functions for sensitivity analysis\n(one per line)')},'mcarlo',[10,60],{s});
        if isempty(answer)
            disp('Cancel pressed');
            return;
        else
            answer=str2cell(answer{1});
        end
    end
    
    drawnow;
    vectout = answer;
    vectout = vectout(~strcmp(vectout, ''));
    %    s = strtrim(sprintf('%s ',s{:}));
    %    s=strrep(s,'  ',' '); %remove double spaces
    %    s=strrep(s,'  ',' ');
    %    f = [0 strfind(s, ' ') length(s) + 1];
    %    vectout = cell(length(f) - 1, 1);
    %    for i = 2:length(f)
    %       vectout{i - 1} = s(f(i - 1) + 1:f(i) - 1);
    %    end
    Nstate = 0;
    g_Y = g_mcarlo.run.Y(:, :, 1);  %#ok<NASGU>
    g_t = g_mcarlo.run.t;
    if hasperm
        g_permanent.Y = g_mcarlo.run.perm(:, :, 1);
        g_permanent.t = g_t;
    end
    i_waitbar(0, length(vectout), 'mcarlo','Analysing results',0.5)
    for i = 1:length(vectout)
        ss = outfun(vectout{i}); % the problem is that this may have a variable dimension
        Nstate = Nstate + size(ss, 2);
    end
    g_Y = [];
    g_t = [];
    mc_output = NaN + zeros(Niter, Nsteps * Nstate);
    g_mcarlo.out = cell(Nstate, 1);
    k = 1;
    l = 1;
    for i = 1:length(vectout)
        ss = outfun(vectout{i}, '-m');
        i_waitbar(1);
        nvar = size(ss, 2);
        for j = 1:nvar
            if nvar == 1
                g_mcarlo.out{l} = vectout{i};
            else
                g_mcarlo.out{l} = sprintf('%s(%d)', vectout{i}, j);
            end
            mc_output(:, k:k + Nsteps - 1) = transpose(reshape(ss(:,j), Nsteps, Niter));
            %I didn't manage to do this in one step
            k = k + Nsteps;
            l = l + 1;
        end
    end
    plabels = cell(Npars, 1);
    for j = 1:Npars
        plabels{j} = g_mcarlo.run.pars{j};
    end
    vlabels = cell(Nsteps * Nstate, 2);
    for i = 1:Nstate
        for j = 1:Nsteps
            vlabels{(i - 1) * Nsteps + j, 1} = g_mcarlo.out{i};
            vlabels{(i - 1) * Nsteps + j,2} = sprintf('%g', g_mcarlo.run.t(j));
        end
    end
    i_waitbar([]);
    %do the multivariate sensitivity analysis
    if ~isempty(g_mcarlo.outfile)
        save(g_mcarlo.outfile, 'g_mcarlo');
    end
    if isfield(args,'dendrpars')
       i_mcsensanal(mc_output, g_mcarlo.run.parvalues, vlabels, plabels, args.dendrpars, g_mcarlo.silent)
    else
       i_mcsensanal(mc_output, g_mcarlo.run.parvalues, vlabels, plabels,{}, g_mcarlo.silent)
    end
end
function ppar =  getparvalue(p, siz)
% if strfind(p.name, '(')||numel(n1)==1&&numel(n2)==1
%     siz = 1;
% else
%     siz = size(p.value);
%     n1 = 1:siz(1);
%     n2 = 1:siz(2);
% end
n1=1:siz(1);
n2=1:siz(2);
switch p.distr
    case 'Uniform'
        ppar = p.min(n1, n2) + rand(siz) * (p.max(n1, n2) - p.min(n1, n2));
    case 'Normal'
        ppar =  drawrandn(siz, p.value(n1,n2), p.sd(n1,n2));
    case 'TruncNormal'
        ppar = -1;
        while ppar < 0
            ppar =  drawrandn(siz, p.value(n1,n2), p.sd(n1,n2));
        end
    case 'LogNormal'
        ppar =  drawrandlogn(siz, p.value(n1,n2), p.sd(n1,n2));
    otherwise
        error('GRIND:mcarlo:UnknownDistr','Distribution %s not supported',p.distr);
end

function res = drawrandlogn(siz, mean, sd)
% draw logNormal distribution with expected mean of mean and expected variance of sd^2
s = log(1 + (sd / mean)^2)^0.5; %this is the sd needed to get a standard deviation of sd
mu = log(mean) - 0.5 * s; %mu is needed to get an expected value of mean
res = exp(randn(siz) .* s + mu);
function res = drawrandn(siz, mean, sd)
res = randn(siz) .* sd + mean;

% function astruc = updatestruc(astruc)
% %update old version of g_mcarlo if necessary
% if isfield(astruc, 'Ys')
%     astruc.run.Y = astruc.Ys;
%     astruc.run.t = astruc.t;
%     astruc.run.perm = astruc.perm;
%     astruc.run.parvalues = astruc.p;
%     astruc=rmfield(astruc,{'pars','p','t','Ys','perm'});
% end

