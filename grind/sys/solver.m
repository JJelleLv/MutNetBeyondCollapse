%SOLVER   Set the solver and settings.
%   For difference equations there is only one method available.
%   Enter the following information:
%   - Solver (rk4 Euler ode45 ode23 ode115 ode15S ode23S ode23T ode23TB)
%   - settings for the solver, see: MATLAB description or <a href="matlab:help rk4">rk4</a>/ <a href="matlab:help euler">euler</a>
%   It is also possible to define and use your own customized solver.
%
%   Usage:
%   solver - brings up a dialog box to select solver
%   solver name par1 par2 - select the solver "name" with  2 additional settings
%   solver -define name hasfixedstep usesode - adds a custom ode function with
%   the name "name" (must be a ode45 compatible m file in the search path), 
%   the two extra settings are used to determine whether the solver has a fixed time step 
%   (default 0) and if it uses the odefile (default 1) 
%   Solver can also be used in model definitions to get information about the current solver:
%   solver('ndays') - returns the number of time units for simulating (see <a href="matlab:help simtime">simtime</a>)
%   solver('step') - returns the fixed time step of the solver (<a href="matlab:help euler">euler</a>, <a href="matlab:help rk4">rk4</a>, <a href="matlab:help back_euler">back_euler</a>) or NaN if the solver
%   has not a fixed step.
%   solver('name') - returns the name of the solver.
%   SOLVER('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'abstol' [number>0] - absolute tolerance for the variable step solver
%     'addmode' [logical] - change the addmode of the model see also <a href="matlab:help addmode">addmode</a>
%     'backwards' [logical] - run backwards (see also <a href="matlab:help backw">backw</a>)
%     'fixedstep' [logical] - the defined solver has a fixed step
%     'name' [string] - name of the solver
%     'nonnegative' [logical] - non-negative option see also -n
%     'nonnegndx' [integer>0] - index of the variables that are non-negative
%     'reltol' [number>0] - relative tolerance for the variable step solver      
%     'stepsize' [number>0] - step for the fixed step solver
%     'usesode' [logical] - the defined solver uses the ode file or not?
%     'vectorized' [logical] - the defined solver can be run vectorized?
%   SOLVER('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - list the current solver
%     '-d' - brings up a dialog box to define a customized solver
%     '-defaultopt' - get the default options (used internally)
%     '-l' - get a list of all solvers
%     '-n' NONNEGATIVE - non-negative solutions only (not available in all solvers, see MATLAB Help).
%     '-n' off - sets the non-negative option off.
%     '-n' [2 3] - non-negative solutions for state variables 2 and 3 only.
%     '-properties' - get all properties of a solver
%
%   
%   See also rk4, euler, <a href="matlab:help ode45">ode45</a>, <a href="matlab:help ode23">ode23</a>, <a href="matlab:help ode113">ode113</a>, <a href="matlab:help ode15s">ode15s</a>, <a href="matlab:help ode23s">ode23s</a> 
%   <a href="matlab:help ode23t">ode23t</a>, <a href="matlab:help ode23tb">ode23tb</a>, <a href="matlab:help odeset">odeset</a>, <a href="matlab:help odeget">odeget</a>, c.ode45, c.euler, c.rk4, csolvers
%
%   Reference page in Help browser:
%      <a href="matlab:commands('solver')">commands solver</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [res] = solver(varargin)
%(solv, p1, p2, p3, p4)
global g_grind;
wrappers= {'dde23','ddesol';'ddesd','ddesolsd';'ode15i','i_ode15sol';'differ','i_differ';'c.','csolvers'};
if nargin==1 %get some information about the solver during run
    f=find(strcmpi(varargin{1}, {'ndays','name','step','reltol','abstol'}));
    if ~isempty(f)
        if f==1
            res=g_grind.ndays;
            return;
        elseif f==2
            i = find(strcmp(wrappers(:, 2), g_grind.solver.name));
            if isempty(i)
                res = g_grind.solver.name;
            else
                res = wrappers{i, 1};
                if strcmp(res, 'c.')
                    res = ['c.' g_grind.solver.csolver];
                end
            end
            return
        elseif f==3
            res = g_grind.solver.opt.StepSize;
            if isempty(res)
                res=nan;
            end
            return;
        elseif f==4
            res = g_grind.solver.opt.RelTol;
            if isempty(res)
                res=nan;
            end
            return;
        elseif f==5
            res = g_grind.solver.opt.AbsTol;
            if isempty(res)
                res=nan;
            end
            return;
        end
    end
end
if nargin==0
    i_parcheck;
    i_solverdlg;
    return;
end

if g_grind.solver.haslags
    if g_grind.dde.isvariable
        solverlist = {'ddesd'};
        fixedstep =  0;
        usesodefile =  1;
        usesJac= 0;
        implicit= 0;
        nonneg =  false;
    else
        solverlist = {'dde23','ddesd'};
        fixedstep = [0 0];
        usesodefile = [1 1];
        usesJac=[0 0];
        implicit=[0 0];
        nonneg = [false false];
    end
elseif g_grind.solver.isdiffer
    solverlist = {'differ','c.differ'};
    fixedstep = [1 1];
    usesodefile  = [1 0];
    usesJac=[0 0];
    implicit=[0 0];
    nonneg = [true true];
    g_grind.solver.opt.StepSize = 1;
elseif g_grind.solver.isimplicit
    solverlist = {'ode15i'};
    fixedstep = 0;
    usesodefile  = 1;
    usesJac=1;
    implicit=1;
    nonneg = true;
elseif g_grind.solver.eulerneeded&&isfield(g_grind.solver,'dwiener')
    solverlist = {'euler','c.euler','stochast_heun'};
    fixedstep =   [1 1 1];
    usesodefile = [1 0 1];
    usesJac=      [0,0,0];
    implicit=      [0,0,0];
    nonneg=         [ true, true,true];
elseif g_grind.solver.eulerneeded
    solverlist = {'euler','c.euler'};
    fixedstep =   [1 1];
    usesodefile = [1 0];
    usesJac=      [0,0];
    implicit=      [0,0];
    nonneg=         [ true, true];
else
    solverlist = {'rk4','euler', 'heun','back_euler','ode45','ode23','ode113','ode15s','ode15i', 'ode23s', 'ode23t', 'ode23tb','ode78','ode87','c.euler','c.rk4','c.ode45','bv78'};
    fixedstep =   [1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0];
    usesodefile = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1];
    usesJac=      [0,0,0,1,0,0,0,1,1, 1, 1, 1,0,0,0,0,0 0];
    implicit=      [0,0,0,0,0,0,0,0,1, 0, 0, 0,0,0,0,0,0 0];
    nonneg=         [ true,true,true,    false,         true,  true,   true,    true,   false,     false ,   true,     true,     false,false,true,true,true,false];
end
%wrappers= {'dde23','ddesol';'ddesd','ddesolsd';'ode15i','i_ode15sol';'differ','i_differ';'c.','csolvers'};
if isfield(g_grind.solver, 'customsolver')
    for i = 1:length(g_grind.solver.customsolver)
        solverlist = [solverlist, {g_grind.solver.customsolver{i}.name}];
        usesodefile = [usesodefile g_grind.solver.customsolver{i}.usesodefile];
        fixedstep = [fixedstep g_grind.solver.customsolver{i}.hasfixedstep];
        nonneg = [nonneg g_grind.solver.customsolver{i}.handles_nonnegative];
    end
end
fieldnams={'name', 's', 'name of the solver','';...
   'stepsize', 'n>0', 'step for the fixed step solver',0.1;...
   'reltol', 'n>0', 'relative tolerance for the variable step solver',5e-5;...
   'abstol', 'n>0', 'absolute tolerance for the variable step solver',1e-7;...
   'nonnegative', 'l', 'non-negative option see also -n',false;...
   'vectorized', 'l', 'the defined solver can be run vectorized?',true;...
   'nonnegndx', 'i>0', 'index of the variables that are non-negative',[];...
   'backwards', 'l', 'run backwards (see also <a href="backw.htm">backw</a>)',false;...
   'addmode', 'l', 'change the addmode of the model see also <a href="addmode.htm">addmode</a><br>',false;...
   'fixedstep', 'l', 'the defined solver has a fixed step',false;...
   'usesode', 'l', 'the defined solver uses the ode file or not?',true}';
args=i_parseargs(fieldnams,...
    ['if(hasoption(''-n'')),if(argtype(2,''l'')),deffields=''nonnegative'';else,deffields=''nonnegndx'';end;'...
    'elseif(hasoption(''-d'')),deffields=''name,fixedstep,usesode,nonegative'';else,deffields=''name,reltol,abstol'';end'],...
    '-defaultopt,-properties,-n,-d,-l,-?',varargin);
if ~isfield(args,'name')
    args.name=g_grind.solver.name;
end
if isfield(args,'vectorized')
    if args.vectorized
        g_grind.solver.opt.Vectorized='on';
    else
        g_grind.solver.opt.Vectorized='off';
    end
end
if strcmp(args.name,'list')||any(strcmp(args.opts, '-l'))
    res = solverlist;
    return
end
if any(strcmp(args.opts, '-properties'))
    if ~isfield(args,'name')
        solvername=g_grind.solver.name;
    else
        solvername=args.name;
    end
    i = find(strcmpi(wrappers(:, 2), solvername));
    res.name=solvername;
    if strcmp(solvername,'csolvers')
        res.name=['c.' g_grind.solver.csolver];
        f=strcmp(res.name,solverlist);
    elseif ~isempty(i)
        f=strcmp(wrappers(i,1),solverlist);
    else
        f=strcmpi(solvername,solverlist);
    end
    if ~any(f)
        error('grind:solver','solver not found');
    end
    res.hasfixedstep=fixedstep(f);
    res.usesodefile=usesodefile(f);
    res.usesJac=usesJac(f);
    res.implicit=implicit(f);
    res.handles_nonnegative=nonneg(f);
    if strcmp(solvername,g_grind.solver.name)
        if res.hasfixedstep
            res.StepSize=g_grind.solver.opt.StepSize;
        else
            res.RelTol=g_grind.solver.opt.RelTol;
            res.AbsTol=g_grind.solver.opt.AbsTol;
        end
        if res.handles_nonnegative
            res.NonNegative=g_grind.solver.opt.NonNegative;
        end
        res.backwards=g_grind.solver.backwards;
        res.addmode=g_grind.solver.addmode;
        if g_grind.solver.isdiffer
            res.iters=g_grind.solver.iters;
        end
    end
    return;
elseif any(strcmp(args.opts, '-defaultopt'))
    res=odeset;
    res.NonNegative=[];
    %     i = find(strcmpi(wrappers(:, 2), g_grind.solver.name));
    %     if ~isempty(i)
    %         f=strcmp(wrappers(i,1),solverlist);
    %     else
    %         f=strcmpi(g_grind.solver.name,solverlist);
    %     end
    %   if fixedstep(f)
    if g_grind.solver.isdiffer
        res.StepSize=1;
    else
        res.StepSize=0.1;
    end
    %   end
    res.RelTol=5e-5;
    res.AbsTol=1e-7;
    return
end

i_parcheck;
if any(strcmp(args.opts, '-s'))
    %use a parameter to switch stochasticity off
    if iscell(args.name)
        g_grind.solver.switch_stochast=args.name;
    else
        g_grind.solver.switch_stochast={args.name};
    end
    return;
end
if any(strcmp(args.opts, '-n'))
    if ~isfield(args,'nonnegndx')
        args.nonnegndx=transpose(1:g_grind.statevars.dim);
    end
    if isfield(args,'nonnegative')
        if ~args.nonnegative
            args.nonnegndx=[];
        elseif numel(args.nonnegative)>1
            args.nonnegndx=find(args.nonnegative);
        end
    end
    if max(args.nonnegndx)>g_grind.statevars.dim
        error('grind:solver:nonNegative','Index of nonNegative larger than the number of state variables');
    end
    if min(args.nonnegndx)<1
        error('grind:solver:nonNegative','Index of nonNegative must be larger than 0');
    end
    g_grind.solver.opt.NonNegative=args.nonnegndx;
    
    if length(g_grind.solver.opt.NonNegative)==g_grind.statevars.dim
        disp('All state variables are forced non-negative');
    elseif isempty(g_grind.solver.opt.NonNegative)
        disp('None of the state variables are forced non-negative');
    elseif ~g_grind.statevars.vector
        s=sprintf('%s,',g_grind.statevars.names{g_grind.solver.opt.NonNegative});
        fprintf('Some state variables (%s) are forced non-negative\n',s(1:end-1));
    end
    return
end
if any(strcmp(args.opts, '-d'))
    if nargin == 1
        answer={'','0','1','0'};
        prompt={'Name of the solver file (must be *.m file)','Has a fixed time step?','Uses the ode file?','Can handle option NonNegative?'};
        answer = inputdlg(prompt, 'solver', 1, answer);
        if ~isempty(answer)
            solver('-d', answer{1}, answer{2}, answer{3}, answer{4});
        end
    else
        f = which(args.name);
        if isempty(f)
            error('GRIND:solver:filenotfound','Cannot find a file named "%s.m"\n', args.name);
        end
        f =  strfind(args.name, '.m');
        if ~isempty(f)
            args.name = args.name(1:f(1) - 1);
        end
        if ~isfield(args,'fixedstep')
            args.fixedstep= false;
        end
        if ~isfield(args,'usesode')
            args.usesode= false;
        end
        if ~isfield(args,'nonnegative')
            args.nonnegative= false;
        end
        if ~isfield(g_grind.solver, 'customsolver')
            p = 1;
        else
            p = length(g_grind.solver.customsolver) + 1;
            for i = 1:length(g_grind.solver.customsolver)
                if strcmp(args.name, g_grind.solver.customsolver{i}.name)
                    p = i;
                end
            end
        end
        g_grind.solver.customsolver{p}.name = args.name;
        g_grind.solver.customsolver{p}.hasfixedstep = args.fixedstep;
        g_grind.solver.customsolver{p}.usesodefile = args.usesode;
        g_grind.solver.customsolver{p}.usesJac=0;
        g_grind.solver.customsolver{p}.implicit=0;
        g_grind.solver.customsolver{p}.handles_nonnegative = args.nonnegative;
        fprintf('Solver "%s" registered\n', args.name)
    end
    return;
end
if any(strcmp(args.opts,'-?'))||strcmp(args.name, '?')
    i = find(strcmpi(wrappers(:, 2), g_grind.solver.name));
    if isempty(i)
        args.name = g_grind.solver.name;
        if strcmp(args.name, 'csolvers')
            args.name = ['c.' g_grind.solver.csolver];
        end
    else
        args.name = wrappers{i, 1};
    end
    n = strcmpi(solverlist, args.name);
    if ~fixedstep(n)
        fprintf('solver %s %g %g\n', args.name, g_grind.solver.opt.RelTol, g_grind.solver.opt.AbsTol);
    elseif ~g_grind.solver.isdiffer
        fprintf('solver %s %g\n', args.name, g_grind.solver.opt.StepSize);
    else
        fprintf('solver %s\n', args.name);
    end
    if length(g_grind.solver.opt.NonNegative)==g_grind.statevars.dim
        disp('All state variables are forced non-negative');
    elseif ~isempty(g_grind.solver.opt.NonNegative)&& ~g_grind.statevars.vector
        s=sprintf('%s,',g_grind.statevars.names{g_grind.solver.opt.NonNegative});
        fprintf('Some state variables (%s) are forced non-negative\n',s(1:end-1));
    end
    return;
end

n = find(strcmpi(solverlist, args.name));
if ~isempty(n)
    %     if ~(g_grind.solver.isdiffer || g_grind.solver.haslags)
    %     oldname=g_grind.solver.name;
    g_grind.solver.name = solverlist{n};
    if usesJac(n)
        if g_grind.statevars.dim>100
            warning('grind:solver:jacobian','This solver needs an Jacobian: is probably very slow for the current high dimensional model')
        end
        enterjac('-update',implicit(n));
    else
        g_grind.solver.opt.Jacobian=[]; %for efficiency
    end
    i = find(strcmpi(wrappers(:, 1), g_grind.solver.name));
    if ~isempty(i)
        g_grind.solver.name  = wrappers{i, 2};
    end
    if ~isempty(strfind(g_grind.solver.name, 'c.'))
        g_grind.solver.csolver = g_grind.solver.name(3:end);
        g_grind.solver.name = 'csolvers';
        if ~csolvers('-isok')
            csolvers('-c','-s');
        end
    elseif  ~usesodefile(n)
        warning('GRIND:solver:sync','The solver "%s" uses its own model definition, be sure the GRIND model is in sync',g_grind.solver.name);
    end
    %   end
    if ~fixedstep(n)
        g_grind.solver.opt.StepSize = [];
    end
    if isfield(args,'reltol')
        if fixedstep(n)
            g_grind.solver.opt.StepSize = args.reltol;
        else
            g_grind.solver.opt.RelTol = args.reltol;
        end
    end
    if isfield(args,'abstol')
        g_grind.solver.opt.AbsTol = args.abstol;
    end
    
elseif g_grind.solver.eulerneeded&&g_grind.solver.isstochastic
    s=solver('-l');
    s=sprintf('%s, ',s{:});
    error('grind:solver:eulerneeded','Solver "%s" not available for a stochastic model, but a fixed time step is needed.\nAvailable solvers: %s',args.name,s(1:end-2))
else
    s=solver('-l');
    s=sprintf('%s, ',s{:});
    error('grind:solver:notfound','Solver "%s" not found for this model\nAvailable solvers: %s',args.name,s(1:end-2))
end
