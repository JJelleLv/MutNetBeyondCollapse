%STOCHAST_HEUN   Heun integration of stochastic differential equation
%   Solve a stochastic differential equation using the two-step improved Euler integration (Heun).
%   This simple integration method is much more accurate than the Euler-Maruyama integration.
%  
%   For a detailed description of the method see:
%   Roberts, A.J. (2012) Modify the improved Euler scheme to integrate stochastic differential equations. <a href="https://arxiv.org/abs/1210.0933v1">arXiv:1210.0933v1</a>
%  
%   The odefunction needs to be split in drift and diffusion terms (each a function handle). This is 
%   done automatically in grind if <a href="matlab:help dwiener">dwiener</a> is used (drift and diffusion terms should add to the result). Alternatively a cell 
%   list of both functions needs to be given (first: drift, second: diffusion (can return a matrix if there are more than 
%   one terms of dwiener per equation))
%
%   This method can solve a Stratonovich or Ito SDE. The default is Ito, if
%   options has a field SDEmethod named 'Stratonovich', the SDE will be solved in the Stratonovich
%   style. If the field is empty or 'Ito' the Ito solution will be used (small difference in practice).
%
%   
%   See also dwiener, solver, euler, heun
%
%   Reference page in Help browser:
%      <a href="matlab:commands('stochast_heun')">commands stochast_heun</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [tout, yout] = stochast_heun(funcs, tspan, y0, options)
global g_grind;
if ~isa(funcs,'cell')
    funcs={i_getodehandle('drift'),i_getodehandle('diffusion'),i_getodehandle('jumps')};
end
if (nargin < 4)
    %if there are no valid options use default stepsize
    delta = 0.3;
    options.StepSize =0.3;
else
    %use the option StepSize as step size
    if ~isfield(options,'StepSize')
        options.StepSize=options.MaxStep;
    end
    if isempty(options.StepSize)
        options.StepSize=0.1;
    end
    
    delta = options.StepSize;
end
if nargin < 4
    nonNegative=[];
    haveOutputFcn = false;
    isIto=true;
else
    nonNegative=odeget(options,'NonNegative',[],'fast');
    haveOutputFcn=~isempty(odeget(options,'OutputFcn',[],'fast'));
    outputFcn=odeget(options,'OutputFcn',[],'fast');
    SDEmethod=odeget(options,'SDEmethod',[],'fast');
    if isempty(SDEmethod)
        isIto=true;
    else
        isIto=strcmpi('ito',SDEmethod);
    end
end
driftfun=funcs{1};
diffusionfun=funcs{2};
if length(funcs)==3
    jumpsfun=funcs{3};
else
    jumpsfun=[];
end
anyNonNegative=~isempty(nonNegative);

if haveOutputFcn
    feval(outputFcn,tspan,y0,'init');
end

% Test that tspan is internally consistent.
tspan = tspan(:);
ntspan = length(tspan);
if ntspan == 1
    t0 = 0;
    next = 1;
else
    t0 = tspan(1);
    next = 2;
end
tfinal = tspan(ntspan);
if t0 == tfinal
    error('GRIND:stochast_heun:tpan','The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t0);
if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan-1)) <= 0)
    error('GRIND:stochast_heun:tspan','The entries in tspan must strictly increase or decrease.');
end
t = t0;
y = y0(:);
neq = length(y);
%adapt delta if there should be given more output
step_tspan=median(diff(tspan));
delta=min(delta,step_tspan);
g_grind.solver.opt.StepSize=delta;
% Set the output flag.

outflag = ntspan > 2;                          % output only at tspan points

% Allocate memory if we're generating output.
delta = delta * tdir;

if nargout > 0
    if outflag                        % output only at tspan points
        tout = tspan;
        outflag=delta~=step_tspan;
        yout = zeros(ntspan, neq);
    else
        tout = transpose(t0:delta:tfinal);
        if tout(end)<tfinal %if tfinal cannot divided in delta's
            tout(end+1)=tfinal;
        end
        
        yout = zeros(size(tout, 1),neq);
    end
    nout = 1;
    tout(nout) = t;
    yout(nout, :) = transpose(y);
end
%
%MAIN LOOP
%evaluate the odefunction for the next time steps
%fold=[];
running=1;
sqrt_delta=sqrt(delta);
while running
    if isIto
        Sk=sign(rand(1)-0.5);
    else
        Sk=0;
    end
    drift=feval(driftfun, t, y);
    diffus=feval(diffusionfun, t, y);
    deltaW=randn(size(diffus))*sqrt_delta;
    if isfield(g_grind.solver.dwiener,'L')
        for i=1:size(deltaW,3)
            deltaW(:,i)=deltaW(:,i).'*g_grind.solver.dwiener.L(:,:,i);
        end
    end
    f1=delta*drift+sum((deltaW - Sk*sqrt_delta).*diffus,2);
    tnew = t + delta;
    drift=feval(driftfun, tnew, y+f1);
    diffus=feval(diffusionfun, tnew, y+f1);
    f2=delta*drift+sum((deltaW + Sk*sqrt_delta).*diffus,2);
    ynew = y + (f1+f2)/2;
    if ~isempty(jumpsfun)
        ynew=ynew+sum(feval(jumpsfun,t,y)*delta,2);
    end
    if anyNonNegative
        ynew(nonNegative) = max(ynew(nonNegative),0);
    end
    if tnew>=tfinal
        running=0;
    end
    if ~outflag                 % computed points, no refinement only the last value
        nout = nout + 1;
        if ~running
            t1=tout(nout);
            yout(nout, :)=transpose((y+(ynew-y)./(tnew-t).*(t1-t)));
        else
            yout(nout, :) = transpose(ynew);
        end
        tnew=tout(nout);
        if haveOutputFcn
            if feval(outputFcn,tnew,transpose(ynew),'')
                running = false;
            end
        end
    elseif (tdir * (tnew - tspan(next)) >= 0) % at tspan, tspan assumed to be larger than delta
        nout = nout + 1;
        t1=tout(nout);
        y1=transpose((y+(ynew-y)./(tnew-t).*(t1-t)));
        yout(nout, :) = y1;
        next = next + 1;
        if haveOutputFcn
            if feval(outputFcn,t1,y,'')
                running = false;
            end
        end
    end
    y = ynew;
    t = tnew;
end
g_grind.solver.opt.StepSize=options.StepSize;
% if nout<length(tout)
%    tout=tout(1:nout);
%    yout=yout(1:nout,:);
% end

