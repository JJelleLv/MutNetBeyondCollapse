%TIMESENS   Performs local sensitivity analysis of parameters/initial conditions in time
%  This algorithm does for each parameter two precise runs with a slightly different 
%  parameter or initial condition. It then computes the partial derivative (dV/dp) 
%  of the state variables V with respect to this parameter/initial condition p during the run. 
%  The results are shown with time or with timesens -t that implements a cluster
%  analysis of the parameters (Klepper 1998, see also <a href="matlab:help mcarlo">mcarlo</a>). Optionally the 
%  elasticities (dV/dp*p/V) instead of derivatives are used.
%
%  Usage:
%  TIMESENS - opens a dialog with all settings and runs also TIMESENS -T
%  TIMESENS PAR - adds PAR in the time plots
%  TIMESENS PAR PAR2 PAR3 etc - analyses PAR1-3
%  TIMESENS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'dendrpars' [equation] - List of parameters that are included in the dendrogram. (used with option '-t')
%     'disturbance' [number>0] - Disturbance for numerical derivatives (default=1E-8)
%     'elasticity' [logical] - Plot elasticities instead of derivatives
%     'maxdim' [integer>=0] - Limit the analysis to the first MAXDIM state variables (default=20)
%     'maxtime' [number] - Stop analysis after maxtime seconds (default=NaN)
%     'nsteps' [number>0] - Show the derivatives with nstep time units
%     'output' [string or cell] - Derivative or equation for output of timesens, for instance d[X]/d[p] (used internally)
%     'pars' [general] - List of parameters for calculating sensitivities
%     'symbolic' [logical] - Use symbolic toolbox
%     't' [number] - Time values for the '-t' option
%   TIMESENS('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - list the parameters used
%     '-all' - adds all parameters and state variables
%     '-c' - clears all sensitivities
%     '-e' - plots elasticities instead of derivatives.
%     '-e-' - plots derivatives (default)
%     '-o' VAR, PAR - returns the text used for <a href="matlab:help out">out</a> (timesens('d(VAR)/d(PAR')).
%     '-p' - use paranal output for getting the data
%     '-s' - stabilize before determining the sensitivities.
%     '-sym' - use the symbolic toolbox to derive semi analytical equations.
%     '-t' - analyse the sensitivities (table/figures/clustering).
%     '-t' STARTT ENDT NSTEPS - show the derivatives in the period STARTT-ENDT with NSTEPS time units.
%     '-t' [1:10:1000] - show the derivatives in the times [1:10:1000]
%     '-u' - update all sensitivities
% 
% 
%  See also time, simtime, out, mcarlo
%
%   Reference page in Help browser:
%      <a href="matlab:commands('timesens')">commands timesens</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [result, vtitles, ptitles] = timesens(varargin)
global g_Y g_t t g_grind g_sens g_paranal;
evalin('base','global g_sens');
defaultdisturb = 1E-8;
%i_parcheck;
if isempty(varargin)
    i_parcheck;
    hassymbolic= i_hastoolbox('symbolic');
    prompt = {'Parameter(s) (list with spaces):','Symbolic derivatives (Y/N)','Elasticity (Y/N)',...
        'Disturbance for numerical derivatives:','Start time:', ...
        'Number of time units to simulate', ...
        'Number of values for output'};
    step = g_grind.tstep;
    if isnan(step)
        step = g_grind.ndays;
    end
    answer={'',iif(hassymbolic,'Y','N/A'),'N',num2str(defaultdisturb),num2str(t),num2str(g_grind.ndays),num2str(step)};
    if isempty(g_sens)||~isfield(g_sens, 'pars')
        vars = [g_grind.pars g_grind.statevars.names];
        answer{1} = strtrim(sprintf('%s ', vars{:}));
    else
        answer{1} = strtrim(sprintf('%s ', g_sens.pars{:}));
        if ~isfield(g_sens,'symbolic')
            g_sens.symbolic=i_hastoolbox('symbolic');
        end
        if ~g_sens.symbolic
            answer{2}='N';
        end
        if g_sens.elasticity
            answer{3} = 'Y';
        end
        answer{4} = num2str(mean(g_sens.reldisturb));
    end
    answer = inputdlg(prompt, 'timesens', 1, answer);
    if isempty(answer) || isempty(answer{1})
        return;
    end
    answer{1}=str2cell(strrep([answer{1} ' '],' ',sprintf('\n'))); %#ok<SPRINTFN>
    if strncmpi(answer{3}, 'Y', 1)
        c=[{'-e'}, transpose(answer{1})];
        timesens(c{:});
    else
        timesens(answer{1}{:});
    end
    if strncmpi(answer{2}, 'Y', 1)
        %symbolic
        timesens('symbolic','on');
    else
        timesens('symbolic','off');
        timesens('disturbance', answer{4});
    end
    timesens('-t', [i_checkstr(answer{5}), i_checkstr(answer{6})], answer{7})
    return;
else
    fieldnams={'pars', '', 'List of parameters for calculating sensitivities','';...
   'nsteps', 'n>0', 'Show the derivatives with nstep time units',1000;...
   'output','s#c','Derivative or equation for output of timesens, for instance d[X]/d[p] (used internally)','';...
   't', 'n', 'Time values for the ''-t'' option',[0:100:1000];...
   'elasticity', 'l', 'Plot elasticities instead of derivatives',false;...
   'symbolic', 'l', 'Use symbolic toolbox', i_hastoolbox('symbolic');...
   'disturbance', 'n>0', 'Disturbance for numerical derivatives (default=1E-8)',1E-8;...
   'dendrpars', 'q', 'List of parameters that are included in the dendrogram. (used with option ''-t'')','';...
   'maxdim', 'i>=0', 'Limit the analysis to the first MAXDIM state variables (default=20)',20;...
   'maxtime', 'n', 'Stop analysis after maxtime seconds (default=NaN)',NaN}';
    args=i_parseargs(fieldnams,'if(hasoption(''-t'')),deffields=''t,nsteps'';elseif(hasoption(''-o'')),deffields=''pars(+)'';elseif argtype(1,''p+#v''),deffields=''pars(+)'';else,deffields=''output'';end',...
        '-all,-sym,-s,-c,-u,-t,-e-,-e,-o,-?,-p',varargin);
end
i_parcheck;
if any(strcmp(args.opts, '-o'))
    if ~isfield(args,'pars')||length(args.pars)==1
        error('grind:timesens:output','Expected a list of state variables and one parameter/variable for the option "-o"')
    end
    avar = args.pars(1:end-1);
    if ~iscell(avar)
        avar={avar};
    end
    result=cell(size(avar));
    apar = args.pars{end};
    ipar=i_getno(apar);
    if ipar.isvar
        apar = sprintf('%s(0)', apar);
    end
    if ~isfield(g_sens,'elasticity')||~g_sens.elasticity
        for i=1:length(avar)
            result{i}=sprintf('timesens(''d[%s]/d[%s]'')',avar{i}, apar);
        end
    else
        no=i_getno(args.pars{end});
        if no.isvar
            apar2= sprintf('val(''%s'')',args.pars{end});
        else
            apar2 = apar;
        end
        for i=1:length(avar)
            if strncmp(avar{i},'#',1)
                result{i}=sprintf('timesens(''d[%s]/d[%s]'')*%s/timesens(''statevar([%s@%s])'')', avar{i}, apar, apar2,avar{i},apar);
            else
                result{i}=sprintf('timesens(''d[%s]/d[%s]'')*%s/%s', avar{i}, apar, apar2,avar{i});
            end
        end
    end
    return;
end
if any(strcmp(args.opts,'-sym'))
    args.symbolic=true;
end
if any(strcmp(args.opts,'-e'))
    args.elasticity=true;
elseif any(strcmp(args.opts,'-e-'))
    args.elasticity=false;
end
if isfield(args,'maxdim')
    g_sens.maxdim=args.maxdim;
elseif ~isfield(g_sens,'maxdim')
    g_sens.maxdim=20;
end
if isfield(args,'maxtime')
    g_sens.maxtime=args.maxtime;
elseif ~isfield(g_sens,'maxtime')
    g_sens.maxtime=NaN;
end
if isfield(args,'elasticity')
    g_sens.elasticity=  args.elasticity;
    update_out;
    if g_sens.elasticity
        disp('elasticity on')
    else
        disp('elasticity off')
    end
elseif ~isfield(g_sens,'elasticity')
    g_sens.elasticity=false;
end


if any(strcmp(args.opts, '-?'))%display
    if ~isfield(g_sens,'pars')
        disp('No paremeters selected for timesens');
        return;
    end
    pars=sprintf('%s,',g_sens.pars{:});
    fprintf('Parameters: %s\nElasticity : %d\n',pars(1:end-1),g_sens.elasticity);
    if g_grind.solver.isdiffer
        sh='(t+1)';
    else
        sh='''';
    end
    if isfield(g_sens,'symbolic')&&g_sens.symbolic
        disp('Derivatives:')
        dispdiff=displaydiff(g_sens.symdiff);
        for i=1:size(g_sens.symdiff,2)
            for j=1:size(g_sens.symdiff,1)
                fprintf('Sens_%s_%s%s = %s\n',g_sens.pars{i},i_statevars_names(j),sh,dispdiff{j,i})
            end
        end
        [N1,odehan]=initsymbolic(i_initvar);
        symdiffs = odehan(t, N1);
        symdiffs = symdiffs(g_grind.statevars.dim+1:end);
        numdiffs=zeros(size(g_sens.symdiff));
        for i=1:length(g_sens.pars)
            apar=evalin('base',g_sens.pars{i});
            iX=i_getno(g_sens.pars{i});
            if iX.isvar
                numdiffs(:,i)=NaN; %I think I cannot test sensitivity to statevars
            else
                h1=i_getodehandle(8,sprintf(',%s',g_sens.pars{i}));
                numdiffs(:,i)=numjac(h1,t,apar,h1(t,apar),1E-6,[],false);
            end
        end
        numdiffs=numdiffs(:);
        mdiff= max(abs(symdiffs(~isnan(numdiffs))-numdiffs(~isnan(numdiffs))));
        if isempty(mdiff)
            mdiff=NaN;
        end
        if ~isnan(mdiff)&&mdiff>0.001
            disp([symdiffs,numdiffs])
            error('grind:timesens:symbolic','Numerical and analytical results differ by > %g',mdiff);
        else
            fprintf('OK max difference between numerical and analytical (single step)= %g\n',mdiff);
        end
    else
        dist=sprintf('%g, ',g_sens.reldisturb(:));
        fprintf('Disturbances: %s\n',dist);
    end
    return;
end
if isfield(args,'symbolic')
    if args.symbolic&&g_grind.statevars.vector
        warning('grind:timesens','Symbolic derivatives not supported for vector notation');
        g_sens.symbolic=false;
    elseif args.symbolic&&g_grind.solver.isimplicit
        warning('grind:timesens','Symbolic derivatives not supported for implicit models (DAE)');
        g_sens.symbolic=false;
    else
        g_sens.lastsettings=[];
        g_sens.symbolic = args.symbolic;
        if ~isfield(args,'pars')
            if isempty(g_sens)||~isfield(g_sens,'pars')
                args.pars=[g_grind.pars,g_grind.statevars.names];
            else
                args.pars=g_sens.pars;
            end
        end
        timesens(args.pars{:});
        if g_sens.symbolic
            i_model2mupad(true);
            if ~isempty(g_grind.syms.errormsg)
                g_sens.symbolic=false;
                error('grind:timesens','%s\n',g_grind.syms.errormsg{:});
            end
            g_sens.handles={};
            allpars=[g_grind.pars,g_grind.statevars.names];
            ndx=ismember(allpars,g_sens.pars);
            g_sens.pars=allpars(ndx);
            g_sens.symdiff=g_grind.syms.Sensitivp(:,ndx);
            if any(~ndx)
                fndx=find(ndx);
                for k=1:length(fndx)
                    for m=1:g_grind.statevars.dim
                        oldsensnr = g_grind.statevars.dim * fndx(k) + m;
                        newsensnr = g_grind.statevars.dim * k + m;
                        %replace statevar symbols with g_X1(i)
                        g_sens.symdiff=regexprep(g_sens.symdiff,sprintf('(?<![a-zA-Z_0-9])g_X1[(]%d,:[)]',oldsensnr),sprintf('g_X1(%d,:)',newsensnr));
                    end
                end
            end
            timesens('-?')
            return;
        end
    end
end
if any(strcmp(args.opts, '-s'))%stabilize and keep the last value
    if isfield(args,'nsteps')
        nsteps = args.nsteps;
    else
        nsteps = 1000;
    end
    oldstep = g_grind.tstep;
    oldndays = g_grind.ndays;
    try
        g_grind.tstep = 2;
        g_grind.ndays = nsteps;
        g_grind.checks.lastsettings = [];
        doupdate;
        ke;
        if isfield(g_sens, 'Ys')
            g_sens.init = g_sens.Ys;
        end
        fprintf('Simulated %d days.\n', nsteps);
        g_grind.tstep =  oldstep;
        g_grind.ndays =  oldndays;
    catch err
        g_grind.tstep =  oldstep;
        g_grind.ndays =  oldndays;
        rethrow(err);
    end
    return;
end
if any(strcmp(args.opts, '-t'))
    if ~isfield(g_sens,'pars')
        error('grind:timesens:nopars','No parameters selected')
    end
    if isfield(args,'t')&&numel(args.t)>2
        ts =  args.t;
    elseif isfield(args,'t')&&isfield(args,'npoints')
        ts = args.t(1):(args.t(2) -args.t(1))/ args.npoints:args.t(2);
    elseif ~isnan(g_grind.tstep)
        ts = (t:g_grind.ndays / g_grind.tstep:t + g_grind.ndays);
    else
        ts = (t:g_grind.ndays / 100:t + g_grind.ndays);
    end
    N0 = i_initvar;
    if  i_settingschanged(N0, ts(end)-ts(1))
        disp('running model');
        if ts(end)-ts(1)>g_grind.ndays
            g_grind.ndays=ts(end)-ts(1);
        end
        i_ru(t, g_grind.ndays, N0, 1);
    end
    if settingschanged
        doupdate;
    end
    dim = size(g_Y, 2);
    result1 = zeros(size(g_Y, 1), dim * length(g_sens.pars));
    titles1 = cell(1, dim * length(g_sens.pars));
    s_matrix = zeros(length(g_sens.pars), dim * length(ts));
    vlabels = cell(dim * length(ts), 2);
    plabels =  transpose(g_sens.pars);
    nt = length(ts);
    statevars_names=i_statevars_names;
    tnames=regexp(sprintf('%g,',ts),',','split');
    for j =  1:dim
        for n = 1:nt
            vlabels{n + nt * (j - 1), 1} = statevars_names{j};
            vlabels{n + nt * (j - 1), 2} = tnames{n};
        end
    end
    warning off MATLAB:interp1:NaNinY
    for k = 1:length(g_sens.pars)
        for j = 1:dim
            if g_sens.elasticity
                titles1{(k - 1) * dim + j} = sprintf('&partial;(%s)/&partial;(%s)*%s/%s', i_statevars_names(j), g_sens.pars{k}, g_sens.pars{k},i_statevars_names(j));
            else
                titles1{(k - 1) * dim + j} = sprintf('&partial;(%s)/&partial;(%s)', i_statevars_names(j), g_sens.pars{k});
            end
            result1(:,(k-1) * dim + j) = outfun(timesens('-o', 'pars',[{i_statevars_names(j)}, g_sens.pars(k)]));
            ttab = interp1(g_t, result1(:, (k - 1) * dim + j), ts);
            for n = 1:nt
                s_matrix(k, (j - 1)  * nt +  n) =  ttab(n);
            end
        end
    end
    warning on MATLAB:interp1:NaNinY
    if ~isfield(args,'dendrpars')
        args.dendrpars={};
    end
    f=find(strcmpi(g_grind.figopts,'visible'));
    silent= ~isempty(f)&&strcmpi(g_grind.figopts{f+1},'off');
    if g_sens.elasticity
        i_klepperen(s_matrix, vlabels, plabels, [], 'Full sensitivity matrix - elasticities d(x)/d(p)*p/x',args.dendrpars,silent)
    else
        i_klepperen(s_matrix, vlabels, plabels, [], 'Full sensitivity matrix - derivatives d(x)/d(p)',args.dendrpars,silent)
    end
    nmax = 20;
    [~,ndx] = sort(abs(result1(end, :)),'descend');
    if length(ndx) > nmax
        ndx = ndx(1:nmax);
        fprintf('The %d most significant sensitivities at the end of the run\n', nmax);
    else
        disp('The sensitivities at the end of the run');
    end
    for j = 1:length(ndx)
        fprintf('%s   %g\n', titles1{ndx(j)}, result1(end, ndx(j)));
    end
    if nargout > 0
        result = s_matrix;
        vtitles = vlabels;
        ptitles = plabels;
    end
    return;
end
if any(strcmp(args.opts, '-u'))
    doupdate;
    return;
end
if any(strcmp(args.opts, '-all'))
    args.pars = [g_grind.pars i_statevars_names];
end
if isfield(args,'disturbance')
    if ~isfield(args,'pars')||isempty(args.pars)
        for i = 1:length(g_sens.pars)
            g_sens.reldisturb(i) = args.disturbance;
        end
    else
        for j=1:length(args.pars)
            i = findadd(args.pars{j});
            if i <= length(g_sens.pars)
                g_sens.reldisturb(i) =  args.disturbance;
            end
        end
    end
end
if any(strcmp(args.opts, '-c'))
    if ~isempty(g_sens)&&isfield(g_sens, 'pars')
        for i = 1:length(g_grind.timevars)
            for k = 1:length(g_sens.pars)
                for j = 1:g_grind.statevars.dim
                    out('-silent',sprintf('-%d',i),'-remove',timesens('-o', 'pars',[{i_statevars_names(j)}, g_sens.pars(k)] ));
                end
            end
        end
    end
    g_sens = [];
    out('-silent','-cleantimevars')
    disp('Cleared all');
    return;
end
if ~(isfield(args,'pars')||isfield(args,'output'))
    return;
end
%****Start with adding parameters here ******
if isfield(args,'output')
par=args.output;
parvar=regexp(par,'[#A-Za-z0-9_\(\),]*(?![\[])','match');
if iscell(parvar)&&length(parvar)==4&&strcmp(parvar{1},'statevar')
    %statevar#1 is the most sensitive statevariable (maxdim is set)
    parvar= strrep(parvar,'(0)','');
    var = parvar{2};
    par1 = parvar{3};
    if var(1)=='#'
        ipar = findadd(par1);
        if settingschanged
            doupdate;
        end
        ivar.no=g_sens.maxndx(str2double(var(2:end)),ipar);
    end
    result = g_Y(:, ivar.no); %just return the most senstivite state variables
    return;
elseif iscell(parvar)&&length(parvar)==2   % d(par) / d(var) or  d(par(i)) / d(var(i))
    parvar= strrep(parvar,'(0)','');
    var = parvar{1};
    par1 = parvar{2};
    if var(1)~='#'
        ivar = i_getno(var);
        if ~ivar.isvar
            error('GRIND:timesens:NoStatevars','%s is not a state variable',var);
        end
    else
        ivar.rank=str2double(var(2:end));
    end
    symbolic = isfield(g_sens, 'symbolic')&&g_sens.symbolic;
    if ~any(strcmp(args.opts,'-p'))
        ipar = findadd(par1);
        if (~symbolic&&(~isfield(g_sens, 'Ys') || (length(g_sens.Ys) < ipar)))||(symbolic&&~isfield(g_sens,'Ysym'))
            if addpar(par1, defaultdisturb)
                g_sens.lastsettings = [];
            else
                return;
            end
        end
        if settingschanged
            doupdate;
        end
        if isfield(ivar,'rank')
            ivar.no=g_sens.maxndx(ivar.rank,ipar);
        end
        if g_sens.symbolic
            result = g_sens.Ysym(:, ipar * g_grind.statevars.dim + ivar.no);
            if size(g_Y, 1) ~= size(g_sens.Ysym, 1)&&size(g_sens.Ysym, 1)>1
                result = interp1(g_sens.tsym, result, g_t);
            end
        else
            yy = g_sens.Ys{ipar};
            if size(g_Y, 1) ~= size(yy, 1)
                result = (yy(:,ivar.no) - g_sens.origY(:, ivar.no)) / g_sens.disturb(ipar);
                result = interp1(g_sens.ts, result, g_t);
            else
                result = (yy(:,ivar.no) - g_sens.origY(:, ivar.no)) / g_sens.disturb(ipar);
            end
        end
    else
        if isempty(g_paranal)
            error('grind:timesens:paranal','First run paranal before using the timesens option "-p"')
        end
        if isfield(g_paranal, 'sens')
            ipar = find(strcmp(par1, g_paranal.sens.pars));
            if ~isempty(ipar)
                sens=g_paranal.sens.Ys{ipar}(:,ivar.no,:);
                result = sens(:);%transpose(reshape(permute(g_paranal.sens.Ys, [2, 1, 3]), size(g_paranal.sens.Ys, 2), size(g_paranal.sens.Ys, 1) * size(g_paranal.sens.Ys, 3)));
%                permute(g_paranal.sens.Ys(:, ivar.no,ipar),[2,1,3]);
                return;
            end
        else
            g_paranal.sens.pars = {};
            g_paranal.sens.reldisturb = [];
            g_paranal.sens.Ys = [];
        end
        g_paranal.sens.pars = [g_paranal.sens.pars, {par1}];
        ipar = length(g_paranal.sens.pars);
        if length(g_paranal.sens.reldisturb) < ipar
            g_paranal.sens.reldisturb(ipar) = 1E-8;
        end
        oldparanal = g_paranal;
        if isfield(g_grind, 'paranal')
            olddlg = g_grind.paranal;
        else
            olddlg = [];
        end
        upar = g_paranal.run.parvalues(:, 1);
        s.par =  {g_paranal.run.pars{1}  ''};
        s.start = [upar(1) 0];   %start value parameter 1 and 2
        s.nend = [upar(end) 10]; %end value parameter 1 and 2
        tt1=g_paranal.run.t(:, 1, g_paranal.run.parvalues == upar(1));
        tt2=g_paranal.run.t(:, 1, g_paranal.run.parvalues == upar(2));
        s.steps = [length(upar) - 1 50];      %number of steps parameter 1 and 2
        s.stabil =  tt2(1) - tt1(end);          %time of stabilizing per step
        s.writing =  tt1(end) - tt1(1);         %number of time units writing
        s.outputtype = 1;        %1 =  unchanged output; 2 = Mean 3 =  Median 4 = Maxima
        s.lines =  0;
        s.silent = true;
        if strcmp(g_paranal.run.pars{1}, par1)
            s.start = s.start * (1 + g_paranal.sens.reldisturb(ipar));
            s.nend = s.nend * (1 + g_paranal.sens.reldisturb(ipar));
            oldpar = [];
        else
            oldpar = evalin('base', par1);
            multassignin('base',  par1, oldpar * (1 +  g_paranal.sens.reldisturb(ipar)));
        end
        try
            paranal(s);
            if ~isempty(olddlg)
                g_grind.paranal = olddlg;
            end
            if ~isempty(oldpar)
                multassignin('base',  par1, oldpar);
                sens = (g_paranal.run.Y - oldparanal.run.Y) / (oldpar * (1 +  oldparanal.sens.reldisturb(ipar)));
            else
                sens = (g_paranal.run.Y - oldparanal.run.Y) ./ permute(repmat((g_paranal.run.parvalues - oldparanal.run.parvalues), size(g_paranal.run.Y,2),size(g_paranal.run.Y,1)),[2 3 1]);
                %repmat((g_paranal.run.parvalues - oldparanal.run.parvalues), 1, size(g_paranal.run.Y, 2));
            end
            g_paranal = oldparanal;
            g_paranal.sens.Ys{ipar} = sens;
            sens=g_paranal.sens.Ys{ipar}(:,ivar.no,:);
            result = sens(:);
        catch err
            if ~isempty(oldpar)
                multassignin('base',  par1, oldpar);
            end
            rethrow(err);
        end
    end
    return;
end
end
i_parcheck;
npar = 0;
if g_sens.elasticity
    s = 'Added elasticities of';
else
    s = 'Added sensitivities of';
end
if isfield(g_sens, 'pars')
    oldn = length(g_sens.pars);
else
    oldn = 0;
end
statevarlist=i_statevars_names;
% if ~g_grind.statevars.vector
%     statevarlist=g_grind.statevars.names;
% else
%     statevarlist={};
%     for i=1:length(g_grind.statevars.vectnames)
%         %very fast:
%         statevar=allelems(g_grind.statevars.vectnames{i},[g_grind.statevars.dims{i}.dim1 g_grind.statevars.dims{i}.dim2]);
%         statevarlist=[statevarlist; statevar(:)];  %#ok<AGROW>
%     end
% end
for j = 1:length(args.pars)
    if addpar(args.pars{j}, defaultdisturb,statevarlist)
        npar = npar + 1;
    end
end
if npar > 0
    if oldn == length(g_sens.pars)
        disp('Parameter(s) were already selected');
    else
        s=sprintf('%s ', s, g_sens.pars{oldn + 1:end});
        disp(s); %#ok<DSPS> % warning not correct sprintf is for getting more pars
    end
    disp('Run <a href="matlab: time">time</a> or <a href="matlab: timesens -t">timesens -t</a> to see the results');
    % updatesolver;
end
function updatesolver
global g_grind g_sens
if any(strcmp(solver('name'), {'c.ode45','ode45'}))&&(~isfield(g_sens,'symbolic')||~g_sens.symbolic)
    if g_grind.solver.opt.AbsTol > 1e-10 || g_grind.solver.opt.RelTol > 1e-10
        % disp('timesens has set solver tolerances to values 1E-11');
        g_sens.oldsolver.name=solver('name');
        g_sens.oldsolver.RelTol=solver('RelTol');
        g_sens.oldsolver.AbsTol=solver('AbsTol');
        solver(g_sens.oldsolver.name, 1E-11, 1E-11);
    end
% elseif isfield(g_sens,'symbolic')&&g_sens.symbolic&&isfield(g_sens,'oldsolver')
%     solver(g_sens.oldsolver.name,g_sens.oldsolver.AbsTol,g_sens.oldsolver.RelTol);
%     g_sens=rmfield(g_sens,'oldsolver') ;
end
if ~isfield(g_sens, 'symbolic')
    g_sens.symbolic=false;
end

function ok = addpar(par, disturb, varlist)
global g_grind g_sens;
if nargin==2
    varlist=i_statevars_names;
    %     if g_grind.statevars.vector
    %         varlist=cell(g_grind.statevars.dim,1);
    %         for i=1:length(g_grind.statevars.vectnames)
    %             elems=allelems(g_grind.statevars.vectnames{i},[g_grind.statevars.dims{i}.dim1,g_grind.statevars.dims{i}.dim2]);
    %             varlist(g_grind.statevars.dims{i}.from:g_grind.statevars.dims{i}.to)=elems;
    %         end
    %     else
    %         varlist=g_grind.statevars.names;
    %     end
end
if g_grind.statevars.dim>g_sens.maxdim
    varlist=regexp(strtrim(sprintf('#%03d ',1:g_sens.maxdim)),' ','split');
end
no = i_getno(par);
ok = 1;
if isempty(no.no)
    fprintf('unknown parameter %s\n', par);
    ok = 0;
    return;
end
if evalin('base',sprintf('isinteger(%s)',par))
    error('grind:timesens:integer','Timesens cannot handle integer parameters (%s)',par);
end
siz= evalin('base',sprintf('size(%s)',par));
if prod(siz) > 1
    if (siz(2) == 1) || (siz(1)==1)
        for j = 1:siz(1) * siz(2)
            ok = addpar(sprintf('%s(%d)', par, j),disturb);
        end
    else
        for i = 1:siz(1)
            for j = 1:siz(2)
                ok = addpar(sprintf('%s(%d,%d)',par,i,j),disturb);
            end
        end
    end
else
    i = findadd(par);
    parexist = ~isempty(g_sens) && isfield(g_sens, 'pars');
    if parexist
        parexist=(i <= length(g_sens.pars));
    end
    g_sens.pars{i} = par;
    g_sens.reldisturb(i) = disturb;
    if ~parexist
        g_sens.lastsettings = [];
        len = length(g_grind.timevars) + 1;
        vars=timesens('-o',varlist,par);
        out('-silent',sprintf('-%d',len),'-add','fun',vars);
    end
end
function res = settingschanged
global g_sens g_t;
if isempty(g_sens) || ~isfield(g_sens, 'lastsettings') || isempty(g_sens.lastsettings)
    res = 1;
    return;
end
settings = i_getsettings;
if ~isempty(g_t)
    settings.solver(1) = g_t(end);
end
res = ~struccmp(settings,g_sens.lastsettings);



% function res = isdifferent(A, B)
% res=~min(size(A) == size(B));
% if isempty(B) && isempty(A)
%    res = 0;
% elseif ~res
%    compar = A == B;
%    res = ~min(compar);
%    if res
%       res = ~min(compar + isnan(A) .* isnan(B));
%    end
% end
function update_out
global g_grind;
for i = 1:length(g_grind.timevars)
    s = g_grind.timevars{i};
    for j = 1:length(s)
        par = s{j};
        if strcontains(par, 'timesens(')
            parvar=regexp(par,'(?<=[\[])[#A-Za-z0-9_\(\),]*(?![\[])','match');
            if iscell(parvar)&&length(parvar)==2   % d(par) / d(var) or  d(par(i)) / d(var(i))
                parvar=strrep(parvar,'(0)','');
                avar = parvar{1};
                apar = parvar{2};
                g_grind.timevars{i}(j) = timesens(struct('opts','-o','pars',{{avar, apar}}));
            end
        end
    end
end
function [N1,odehan]=initsymbolic(N0)
global g_sens;
g_sens.handles{1} = i_getodehandle(1, '');
s = sprintf('%s;', g_sens.symdiff{:});
dim=size(g_sens.symdiff,1);
comm=sprintf('@(t,g_X1)[zeros(%d,size(g_X1,2));%s]',dim,s(1:end-1));
g_sens.handles{2}=evalin('base',comm);
odehan=@(t,g_X1)runsymbolic(t, g_X1, g_sens.handles{1},g_sens.handles{2},dim);
N1 = [N0; zeros(numel(g_sens.symdiff), 1)];
for i = 1:length(g_sens.pars)
    iX = i_getno(g_sens.pars{i});
    if iX.isvar
        N1(i * size(g_sens.symdiff, 1) + iX.no) = 1;
    end
end

function doupdate
global g_Y g_t t g_grind g_sens;
if ~isempty(g_sens)&&isfield(g_sens, 'pars')&&~isempty(g_sens.pars)
    if isempty(g_t)
        ndays=g_grind.ndays;
    else
        ndays=g_t(end)-g_t(1);
        t=g_t(1);
    end
    g_sens.lastsettings = i_getsettings;
    g_sens.lastsettings.solver(1)=ndays;
    if ~isfield(g_sens, 'symbolic')
       g_sens.symbolic=false;
    end
    N0 = i_initvar;
    oldY = g_Y;
    oldt = g_t;
    if g_sens.symbolic
        disp('Updating timesens (symbolic)...');
        [N1,odehan]=initsymbolic(N0);
        oldopt=g_grind.solver.opt;
        g_grind.solver.opt.Jacobian=[];
        g_grind.solver.opt.Vectorized='off';
        if ~isnan(g_sens.maxtime)
            loc_odespeedtest(g_sens.maxtime,[],'setup');
            tic
            g_grind.solver.opt.OutputFcn=@loc_odespeedtest;
        end
        [g_sens.tsym,g_sens.Ysym] = feval(str2func(g_grind.solver.name), odehan,g_t, N1 , g_grind.solver.opt);
        g_grind.solver.opt=oldopt;
        maxsens=max(abs(g_sens.Ysym),[],1);
        maxsens=reshape(maxsens,g_grind.statevars.dim,length(maxsens)/g_grind.statevars.dim);
        maxsens=maxsens(:,2:end);
        [~,g_sens.maxndx]=sort(maxsens,1,'descend');
    else
        try
            disp('Updating timesens (numeric)...');
            oldopt=g_grind.solver.opt;  
            oldtstep = g_grind.tstep;
            updatesolver;
            i_waitbar(0,length(g_sens.pars)+1,'timesens','Updating timesens',1);
            if isnan(g_grind.tstep)&&~g_grind.solver.isdiffer
                g_grind.tstep = max(100, g_grind.ndays);
            end
            i_ru(t, ndays, N0, 0);
            if size(g_Y,1)<ndays
                g_Y(end+1:ndays+1,:)=NaN;
                g_t=(0:ndays).';
            end
            g_sens.lastsettings.solver(1) = g_t(end);
            g_sens.origY = g_Y;
            g_sens.ts = g_t;
            if ~isnan(g_sens.maxtime)
                loc_odespeedtest(g_sens.maxtime,[],'setup');
                g_grind.solver.opt.OutputFcn=@loc_odespeedtest;
            end
%             if isnan(g_grind.tstep)&&~g_grind.solver.isdiffer
%                 g_grind.tstep = max(100, g_grind.ndays);
%                 tic;
%                 i_ru(t, g_grind.ndays, N0, 1);
%                 if size(g_Y,1)<g_grind.ndays
%                     g_Y(end+1:g_grind.ndays+1,:)=NaN;
%                     g_t=(0:g_grind.ndays).';
%                 end
%                 
%                 g_sens.origY = g_Y;
%                 g_sens.ts = g_t;
%                 i_waitbar(1);
%             else
%                 i_waitbar(0,length(g_sens.pars),'timesens','Updating timesens',1);
%             end
            for l_i = 1:length(g_sens.pars)
                oldpar = evalin('base', g_sens.pars{l_i});
                if oldpar==0
                    g_sens.disturb(l_i)=1E-8;
                else
                    g_sens.disturb(l_i) = g_sens.reldisturb(l_i) * oldpar;
                end
                multassignin('base',  g_sens.pars{l_i}, oldpar +  g_sens.disturb(l_i));
                if isempty(g_sens)||~isfield(g_sens, 'init')
                    N0 = i_initvar;
                else
                    N0 = transpose(g_sens.init{l_i}(end, :));
                end
                tic;
                i_ru(t, ndays, N0, 0);
                if size(g_Y,1)<ndays
                    g_Y(end+1:ndays+1,:)=NaN;
                    g_t=(0:ndays).';
                end
                
                i_waitbar(1);
                multassignin('base', g_sens.pars{l_i}, oldpar);
                g_sens.Ys{l_i} = g_Y; %interp1(g_t, g_Y, oldt);
            end
            g_grind.solver.opt=oldopt;
            g_sens.maxndx=zeros(g_grind.statevars.dim,length(g_sens.Ys));
            for i=1:length(g_sens.Ys)
                %what if the sizes are not consistent?
                yy=g_sens.Ys{i};
                if size(yy,1)~=size(g_sens.origY)
                    if size(yy,1)<size(g_sens.origY)
                        [~,ndx]=sort( max(abs(yy(size(yy,1),:)-g_sens.origY(size(yy,1),:))),2,'descend');
                    else
                        [~,ndx]=sort( max(abs(yy(size(g_sens.origY,1),:)-g_sens.origY(size(g_sens.origY,1),:))),2,'descend');
                    end
                else
                    [~,ndx]=sort( max(abs(yy-g_sens.origY)),2,'descend');
                end
                g_sens.maxndx(:,i)=transpose(ndx);
            end
            
            i_waitbar([])
            g_Y = oldY;
            g_t = oldt;
            g_grind.tstep = oldtstep;
            if isfield(g_sens,'oldsolver')
               solver(g_sens.oldsolver)
            end
        catch err
            if isfield(g_sens,'oldsolver')
               solver(g_sens.oldsolver)
            end
            i_waitbar([]);
            %   err=lasterror;
            g_Y = oldY;
            g_t = oldt;
            g_grind.tstep = oldtstep;
            rethrow(err);
        end
    end
end
function status = loc_odespeedtest(t, ~, flag)
persistent maxtime;
if strcmp(flag,'setup')
    maxtime=t;
    return;
end
status = toc>maxtime;
function j = closingbrack(s, i)
if s(i) == '('
    brack1 = '(';
    brack2 = ')';
elseif s(i) == '['
    brack1 = '[';
    brack2 = ']';
end
j = i + 1;
k = 1;
while (k~=1) || (s(j)~=brack2)
    if s(j) == brack1
        k = k + 1;
    end
    if s(j) == brack2
        k = k - 1;
    end
    j = j + 1;
end
function i = findadd(par)
global g_sens;
i = 1;
if ~isempty(g_sens) && isfield(g_sens, 'pars')
    i=find(strcmp(g_sens.pars,par),1);
    if isempty(i)
        i= length(g_sens.pars)+1;
    end
    %     while i <= length(g_sens.pars) && ~strcmp(g_sens.pars{i}, par)
    %         i = i + 1;
    %     end
end
% function g_X2 = runsymbolic(t, g_X1)
% global g_sens;
% g_X2 = g_sens.handles{2}(t, g_X1);
% g_X2(1:size(g_sens.symdiff, 1),:) = g_sens.handles{1}(t, g_X1(1:size(g_sens.symdiff, 1),:));
function g_X2 = runsymbolic(t, g_X1, han1,han2,dim)
g_X2 = han2(t, g_X1); %the jacobians
g_X2(1:dim,:) =han1(t, g_X1(1:dim,:)); %ode function

function s=displaydiff(differ)
global g_grind g_sens;
s=differ;
if g_grind.solver.isdiffer
    sh='(t)';
else
    sh='';
end
for i=1:g_grind.statevars.dim
    s=regexprep(s,sprintf('(?<![a-zA-Z_0-9#])g_X1[(]%d,:[)]',i),[i_statevars_names(i) sh]);
end
h=g_grind.statevars.dim;
for i=1:g_grind.statevars.dim
    for j=1:length(g_sens.pars)
        s=regexprep(s,sprintf('(?<![a-zA-Z_0-9#])g_X1[(]%d,:[)]',j*h+i),sprintf('Sens_%s_%s%s',g_sens.pars{j},i_statevars_names(i),sh));
    end
end
