%DJUMP   Adds a jump process to a stochastic differential equation
%   This is a versatile implementation of a compound Poisson process. You
%   have to define a function that returns the waiting time between two
%   jumps (typically modelled as Poisson process, where intervals are generated 
%   by drawing from an exponential distrbution (exprnd(lambda), where lambda is 
%   the average interval between jumps)) and a function that generates the
%   size of the jumps (can be a normal distribution normrnd(mu,sd)).
%   You can use this function in combination with <a href="matlab:help dwiener">dwiener</a> to
%   create a jump-diffusion process. The solver is by default set to euler integration.
%      
%   Usage:
%   DJUMP(timing,size) where timing is an equation (including arguments) that defines
%    intervals and size is an equation that defines the jumps.
%   DJUMP('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'corr' [number>=0 and number<=1] - correlation in the jump size between the state
%    variables (single value is for the correlation within a vector state variable)
%     'no' [integer>0] - number of the djump function (not used when defining a djump)
%     'size' [string] - equation that generates the size of a jump
%     'synchronize' [integer>=0] - use this to synchronize the jumps of different variables,
%    the djumps with the same number will synchronize (0 means independent). The timing equation of 
%    the factors should be exactly the same (default=0)
%     'timing' [string] - equation that generates the duration of the intervals
%   DJUMP('-opt1','-opt2',...) - Valid command line options:
%     '-i' - used internally to initialize the djump functions
%     '-l' - list the defined jumps
% 
%  Example:
%  DJUMP(exprnd(meantiming),normrnd(meansize,sdsize)) - this creates a compound Poisson
%    process with Gaussian distributed jump sizes. Note that meantiming,meansize and sdsize
%    are added as parameters to the model and should be initialized.
%  DJUMP(exprnd(meantiming),meansize,'synchronize',1) - Poisson process with fixed jumps. To synchronize
%   the timing between state variable 'synchronize' 1 is used (terms with the same number>0
%   will be synchronized).
%     
%
%  See also model, dwiener, setevent
%
%   Reference page in Help browser:
%      <a href="matlab:commands('djump')">commands djump</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function args1=djump(varargin)
global g_grind
fieldnams={'timing', 's', 'equation that generates the duration of the intervals',100;...
   'size', 's', 'equation that generates the size of a jump','';...
   'synchronize', 'i>=0', 'use this to synchronize the jumps of different variables,',0;...
   'corr', 'n>=0&n<=1', 'correlation in the jump size between the state',0;...
   'no', 'i>0', 'number of the djump function (not used when defining a djump)',0}';
args=i_parseargs(fieldnams,'timing,size','-i,-l',varargin);
if strcmp(args.opts,'-i')
    maxsynch=0;
    if isfield(g_grind.solver,'djump')
        for i=1:length(g_grind.solver.djump.args)
            if g_grind.solver.djump.args(i).synchronize>maxsynch
                maxsynch=g_grind.solver.djump.args(i).synchronize;
            end
            g_grind.solver.djump.args(i).nextt=[];
            g_grind.solver.djump.args(i).timing_pars=getvars(g_grind.solver.djump.args(i).timing);
            g_grind.solver.djump.args(i).size_pars=getvars(g_grind.solver.djump.args(i).size);
            g_grind.solver.djump.args(i).timing_fun=makefun(g_grind.solver.djump.args(i).timing,g_grind.solver.djump.args(i).timing_pars,false);
            g_grind.solver.djump.args(i).size_fun=makefun(g_grind.solver.djump.args(i).size,g_grind.solver.djump.args(i).size_pars,g_grind.statevars.vector);
        end
        if maxsynch>0
            g_grind.solver.djump.nextt=cell(maxsynch,1);
            sync_groups=[g_grind.solver.djump.args(:).synchronize];
            ugroups=unique(sync_groups);
            for i=1:length(ugroups)
                if ugroups(i)>0
                    aa1={g_grind.solver.djump.args(sync_groups==ugroups(i)).timing};
                    aa=unique(aa1);
                    if length(aa)>1
                        error('grind:djump:sync','Error in synchronizing jumps, the timing functions within the sync group %d are not the same: \n%s',ugroups(i),sprintf('"%s" ',aa{:}));
                    end
                    if length(aa1)==1
                        warning('grind:djum:sync','Cannot synchronize sync group %d, because the group has only one member',ugroups(i));
                    end
                end
            end
        end
    end
    return;
end
if strcmp(args.opts,'-l')
    for i=1:length(g_grind.solver.djump.args)
        fprintf('Jump process no: %d\n',i)
        fprintf('  Timing of jumps (timing): %s\n',char(g_grind.solver.djump.args(i).timing));
        fprintf('  Jump sizes (size): %s\n',char(g_grind.solver.djump.args(i).size));
    end
    return
end
if nargout>0
    %only if used for the first time, use defaults end give errors
    if ~isfield(args,'synchronize')
        args.synchronize=0;
    end
    if ~isfield(args,'timing')
        %default exponential distribution of the next time of Poisson process
        args.timing='exprnd(100)';
    end
    if ~isfield(args,'size')
        %default normal distribution of the size of Poisson process
        args.size='normrnd(0,10)';
    end
    if ~isfield(args,'corr')
        args.corr=0;
    end
end
if isfield(args,'timing')
    args.timing_pars=getvars(args.timing);
    args.timing_fun=makefun(args.timing,args.timing_pars,false);
end
if isfield(args,'size')
    args.size_pars=getvars(args.size);
    args.size_fun=makefun(args.size,args.size_pars,false);
end
if nargout>0
    args1=args;
else
    if size(g_grind.solver.djump.args)==1
        g_grind.solver.djump.args=mergestructs(g_grind.solver.djump.args,args);
    elseif isfield(args,'no')
        g_grind.solver.djump.args(args.no)=mergestructs(g_grind.solver.djump.args(args.no),args);
    else
        error('grind:djump','More than one djump fields please indicate the number (field "no")')
    end
end
function res=exprnd(lambda,siz) %#ok<DEFNU>
if nargin<2
    siz=1;
end
res=-log(rand(siz))*lambda;

function res=normrnd(mu,sigma,siz) %#ok<DEFNU>
if nargin<3
    siz=1;
end
res=randn(siz)*sigma+mu;

function res=makefun(fun,vars,makesize)
addsize=false;
if makesize
    sizeable_funcs= {'betarnd','binornd','chi2rnd','evrnd','exprnd','frnd','gamrnd',...
        'geornd','gevrnd','gprnd','hygernd','lognrnd','nbinrnd','ncfrnd','nctrnd',...
        'ncx2rnd','normrnd','poissrnd','rand','randg','randn','random','raylrnd','trnd','unifrnd','wblrnd'};
    fun1=parsed_equation(fun);
    f=find(ismember(fun1.fs,sizeable_funcs));
    
    if ~isempty(f)
        for i=1:length(f)
            [pars,j]=fun1.getfunctionpars(f(i));
            if  j(end,end)+1<=length(fun1.fs)&& strcmp(fun1.fs{j(end,end)+1},')')
                if isempty(pars)
                    fun1.fs{j(end,end)+1}='g_size_var)';
                else
                    fun1.fs{j(end,end)+1}=',g_size_var)';
                end
                addsize=true;
            end
        end
        fun=sprintf('%s',fun1.fs{:});
    end
end
s=sprintf('%s,',vars{:});
if length(s)>1
    s=s(1:end-1);
end
if addsize
    res=eval(sprintf('@(%s,g_size_var)%s',s,fun));
else
    res=eval(sprintf('@(%s)%s',s,fun));
end

function res=getvars(fun)
pfun=parsed_equation(fun);
res=symvar(pfun);


