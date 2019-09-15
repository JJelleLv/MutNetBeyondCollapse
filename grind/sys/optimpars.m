%OPTIMPARS   Fit data by optimizing parameters
%   Optimize parameters with the local simplex method or a global search (Shuffled Complex 
%   Evolution (SCE-UA)). First, a data set must be entered using <a href="matlab:help setdata">setdata</a> (a dialog box 
%   for entering data appears if this has not yet been done). Criterion is the sum 
%   of squares of normalized state variables.
%
%   Usage:
%   OPTIMPARS - Enter parameters and options in dialog boxes.
%   OPTIMPARS P1 P2 P3 .. - Optimize the list of parameters P1, P2, P3 using simplex optimization.
%   OPTIMPARS -local P1 P2 P3 .. - use simplex method (default).
%   OPTIMPARS -global P1 [L1 U1] P2 [L3 U3] P3 [L3 U3] .. - using global optimization and ranges for 
%   each parameter [L1=Lower P1, U1=Upper P1].
%   OPTIMPARS -reset - Reset the parameters as before the previous run.
%
%   OPTIMPARS('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'maxfunevals' [integer>0 or string] - maximum number of function evaluations allowed
%     'maxiter' [integer>0 or string] - maximum number of iterations allowed
%     'maxpcgiter' [integer>0] - number of complexes (sub-populations) (global method)
%     'maxtime' [number>=0|isnan(number)] - maximum time for optimization
%     'optmethod' [string] - method of optimization (fminsearch or fglobalmin or full name)
%     'parmaxs' [number] - minimum values of parameters (global method)
%     'parmins' [number] - maximum values of parameters (global method)
%     'pars' [general] - list of parameters and initial values of state variables to optimize, or with ranges
%     'penaltyfun' [normalized ss | ss or function_handle] - the function that determines the penalty (SS sum of squares, normalized is weight independent of scale)) (user function should take g_data)
%     'posonly' [logical] - restrict all parameters to positive values
%     'reseteachstep' [logical] - reset the simulation after each data point
%     'tolfun' [number>0] - termination tolerance on the function value (local method)
%     'tolx' [number>0] - termination tolerance on x
%   OPTIMPARS('-opt1','-opt2',...) - Valid command line options:
%     '-g' - use global method (fglobalmin)
%     '-l' - use global method (fminsearch)
%     '-r' - reset the parameters as before the previous run.
%       
%         
%See also setdata, fglobalmin, <a href="matlab:help fminsearch">fminsearch</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('optimpars')">commands optimpars</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function optimpars(varargin)
global g_data g_grind g_Y;
if ~isfield(g_data, 'options') || isempty(g_data.options)
    opts = optimset('fminsearch');
    g_data.options.Display = opts.Display;
    g_data.options.TolX = 1E-6;
    g_data.options.TolFun = 1E-6;
    g_data.options.MaxTime1=NaN;
    g_data.options.MaxFunEvals1 = opts.MaxFunEvals;
    g_data.options.MaxIter1 = opts.MaxIter;
    g_data.options.penaltyfun='normalized SS';
    g_data.options.PosOnly = 0;
    g_data.options.ResetEachStep  = 0;
end
if ~isfield(g_data, 'options2') || isempty(g_data.options2)
    opts =  optimset('fglobalmin');
    g_data.options2 = opts;
    g_data.options2.MaxTime1=NaN;
    g_data.options2.PosOnly = 0;
    g_data.options2.ResetEachStep = 0;
end
optmethods={'Nelder-Mead simplex (direct search) method', 'Shuffled Complex Evolution (SCE-UA) Duan, et al. 1992'};
if nargin==0
    i_parcheck;
    if ~i_optimparsdlg
        disp('Cancelled');
        return;
        %error('grind:optimpars','Cancelled by user')
    end
else
    fieldnams={'optmethod', 's', 'method of optimization (fminsearch or fglobalmin or full name)';...
   'pars', '', 'list of parameters and initial values of state variables to optimize, or with ranges';...
   'maxtime', 'n>=0|isnan(n)', 'maximum time for optimization';...
   'parmins', 'n', 'maximum values of parameters (global method)';...
   'parmaxs', 'n', 'minimum values of parameters (global method)';...
   'tolx', 'n>0', 'termination tolerance on x';...
   'tolfun', 'n>0', 'termination tolerance on the function value (local method)';...
   'maxfunevals', 'i>0#s', 'maximum number of function evaluations allowed';...
   'maxiter', 'i>0#s', 'maximum number of iterations allowed';...
   'penaltyfun', 'e[normalized SS|SS]#f', 'the function that determines the penalty (SS sum of squares, normalized is weight independent of scale)) (user function should take g_data)';...
   'posonly', 'l', 'restrict all parameters to positive values';...
   'reseteachstep', 'l', 'reset the simulation after each data point';...
   'maxpcgiter', 'i>0', 'number of complexes (sub-populations) (global method)'}';
    args=i_parseargs(fieldnams ,'pars(+)', '-g,-l,-r',varargin);
    i_parcheck;
    if any(strcmp(args.opts,'-g'))
        args.optmethod=optmethods{2};
    end
    if any(strcmp(args.opts,'-l'))
        args.optmethod=optmethods{1};
    end
    if isfield(args,'optmethod')
        if strncmpi(args.optmethod,'n',1)||strcmp(args.optmethod,'fminsearch')
            g_data.optmethod=optmethods{1};
        elseif strncmpi(args.optmethod,'s',1)||strcmp(args.optmethod,'fglobalmin')
            g_data.optmethod=optmethods{2};
        else
            error('grind:optimpars:unknownmethod','Unknown method "%s" for optimpars',args.optmethod)
        end
    else
        g_data.optmethod=optmethods{1};
    end
    if isfield(args,'parmins')
        g_data.parmins=args.parmins;
    end
    if isfield(args,'parmaxs')
        g_data.parmaxs=args.parmaxs;
    end
    allopts={'tolx','tolfun','maxtime','maxfunevals','maxiter','penaltyfun','posonly','reseteachstep','maxpcgiter'};
    allopts1={'TolX','TolFun','MaxTime1','MaxFunEvals1','MaxIter1','penaltyfun','PosOnly','ResetEachStep','MaxPCGIter'};
    f=fieldnames(args);
    ndx=ismember(allopts,f);
    if any(ndx)
        methodnr=find(strcmp(g_data.optmethod,optmethods));
        if methodnr==1
            optfield='options';
        else
            optfield='options2';
        end
        ndx=find(ndx);
        for i=1:length(ndx)
            g_data.(optfield).(allopts1{ndx(i)})=args.(allopts{ndx(i)});
        end
    end
    if ~isfield(args,'pars')
        if ~i_optimparsdlg
            error('grind:optimpars','Cancelled by user')
        end
    else
        if ~iscell(args.pars)
            args.pars={args.pars};
        end
        methodnr=find(strcmp(g_data.optmethod,optmethods));
        if methodnr==2
            if numel(args.pars)>1
                if isnumeric(args.pars{2})||~isnan(all(str2num(args.pars{2}))) %#ok<ST2NM>
                    pars={};
                    g_data.parmins=[];
                    g_data.parmaxs=[];
                    for i=1:floor(numel(args.pars)/2)
                        pars{i}=args.pars{2*(i-1)+1};
                        if isnumeric(args.pars{2*(i-1)+2})
                            vals=args.pars{2*(i-1)+2};
                        else
                            vals=str2num(args.pars{2*(i-1)+2}); %#ok<ST2NM>
                        end
                        if numel(vals)==1
                            vals(2)=vals(1)*2;
                        end
                        g_data.parmins(i)=vals(1);
                        g_data.parmaxs(i)=vals(2);
                    end
                    args.pars=pars;
                end
            end
        end
        g_data.pars=args.pars;
    end
end

%%%  start optimizing %%%%
g_grind.solver.reset_at_data=g_data.options.ResetEachStep;
if ~isfield(g_data,'obs')||isempty(g_data.obs)
    error('grind:optimpars:nodata','No data for optimizing')
end
if ~isfield(g_data,'pars')||isempty(g_data.pars)
    error('grind:optimpars:nodata','No parameter to optimize')
end
g_data.X = [];
g_data.fval = [];
g_data.iter = 0;
g_data.stopped = 0;
g_data.minobs = min(g_data.obs);
g_data.maxobs = max(g_data.obs);

par0 = zeros(1, size(g_data.pars, 2));
for i = 1:size(g_data.pars, 2)
    par0(i) = evalin('base', char(g_data.pars{i}));
end
oldtstep = g_grind.tstep;
g_data.X0 = par0;
g_grind.tstep = NaN;
try
    g_grind.solver.opt.OutputFcn = [];
    numberofvariables = length(g_data.pars); %#ok<NASGU> % needed
    %  numberOfVariables = length(g_data.pars); %version <7.1
    g_data.options.MaxFunEvals  = eval(lower(g_data.options.MaxFunEvals1));
    g_data.options.MaxIter  = eval(lower(g_data.options.MaxIter1));
    H0 = i_optimpardlg;
    if ~isempty(g_grind.figopts)
        set(H0,g_grind.figopts{:});
    end
    methodnr=find(strcmp(g_data.optmethod,optmethods));
    if methodnr==2
        if ~isnan(g_data.options2.MaxTime1)&&~isempty(g_data.options2.MaxTime1)
            tstart=tic;
            MaxTime=g_data.options2.MaxTime1;
            g_data.options2.OutputFcn=@(x, optimValues, state)timefun(MaxTime,tstart,x, optimValues, state);
        else
            g_data.options2.OutputFcn=[];
        end
        [X, fval, found] = fglobalmin('i_optimpars', g_data.parmins,g_data.parmaxs, g_data.options2);
    else
        if ~isnan(g_data.options.MaxTime1)&&~isempty(g_data.options.MaxTime1)
            tstart=tic;
            MaxTime=g_data.options.MaxTime1;
            g_data.options.OutputFcn=@(x, optimValues, state)timefun(MaxTime,tstart,x, optimValues, state);
        else
            g_data.options.OutputFcn=[];
        end
        [X, fval, found] = fminsearch('i_optimpars', par0, g_data.options);
    end
    assignpars(X);
    if ishandle(H0)
        close(H0);
    end
    g_grind.solver.opt.OutputFcn = str2func('i_odespeed');
    g_grind.tstep = oldtstep;
catch err
    %  err=lasterror;
    %if 0
    g_grind.solver.opt.OutputFcn = str2func('i_odespeed');
    g_grind.tstep = oldtstep;
    fval =    g_data.fval;
    if ishandle(H0)
        close(H0);
    end
    found=g_data.stopped == 1;
    if found
        assignpars(g_data.X);
        disp('Optimization stopped by user, using best solution so far');
    elseif g_data.stopped == 2
        assignpars(g_data.X0);
        g_Y = [];
        disp('Optimization cancelled by user, no change in parameters');
%        rethrow(err);
    else
        rethrow(err);
    end
end
if found
    fprintf('Method            : %s\n',g_data.optmethod );
    fprintf('Penalty function  : %s\n',char(g_data.options.penaltyfun));
    adjr2(length(g_data.pars));
    disp(['Sum of squares    : ' num2str(fval)]);
    disp(['No. of iterations : ' num2str(g_data.iter)]);
    disp(' ');
end
i = i_figno('time');
if ishandle(i)
    set(i, 'visible', 'off');
end

function multassignin(ws, name, V, fbrack)
%multassignin supports assignments to elements of arrays (e.g. name='A(1,2)')
%fbrack should be strfind(name, '(') (added for efficiency)
if isempty(fbrack)
    assignin(ws, name, V);
else
    temp = evalin(ws, name(1, fbrack(1) - 1));
    s = ['temp' name(fbrack(1):length(name)) ' = ' num2str(V) ';' ];
    eval(s);
    assignin(ws, name(1, fbrack(1) - 1), temp);
end
function assignpars(x)
global g_data;
g_data.pred = [];
for i = 1:size(g_data.pars, 2)
    fbrack = strfind(g_data.pars{i}, '(');
    multassignin('base', g_data.pars{i}, x(i),fbrack);
end

function stop=timefun(maxtime,tstart,~, ~, ~)
stop=toc(tstart)>maxtime;
