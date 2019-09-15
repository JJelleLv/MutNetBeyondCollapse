function [tout, yout, varargout] = dde_euler(odefile, tspan, y0, options, varargin)
global g_grind;
%DDE23('F',LAGS,HISTORY,TSPAN)
if nargin >= 4
   nonNegative = ~isempty(odeget(options,'NonNegative',[],'fast'));
   if nonNegative
      warning('GRIND:DDE23:NonNegativeIgnored','DDE23 does not constrain solution to be non-negative. Option ''NonNegative'' will be ignored.');
   end
end
lags = zeros(size(g_grind.dde.lags));
varargout = varargin;
for i = 1:length(g_grind.dde.lags)
   lags(i) = evalin('base', g_grind.dde.lags{i});
end
if any(lags<0)
    error('grind:ddesol','All time lags should be constant and positive');
end
if ~isempty(g_grind.solver.history)&&isstruct(g_grind.solver.history)
   g_grind.solver.history.x = g_grind.solver.history.x - g_grind.solver.history.x(end) + tspan(1);
   if ~all(transpose(y0) == ddehist(tspan(1)))
      g_grind.solver.history = y0;
   end
else
   g_grind.solver.history = y0;
end
sol = dde_euler1(odefile, lags, @ddehist, [tspan(1), tspan(length(tspan))], options);
if ~any(isnan(sol.y))
   g_grind.solver.newhist = sol;
else
   g_grind.solver.newhist = [];
end
if length(tspan) > 2
    if tspan(end)>sol.x(end)
        %this can be due to an event, don't want an error then
        j=length(tspan);
        k=length(sol.x);
        while tspan(j)>sol.x(k)&&j>0&&k>0
            tspan(j)=sol.x(k);
            j=j-1;
            k=k-1;
        end
    end
   [yout] = deval(sol, transpose(tspan));
   tout = transpose(tspan);
   yout = transpose(yout);
else
   yout = transpose(sol.y);
   tout = transpose(sol.x);
end

function y = ddehist(t)
global g_grind;
if isstruct(g_grind.solver.history)
    if t<g_grind.solver.history.x(1)
      y=g_grind.solver.history.y(1);
    else
      y = transpose(deval(g_grind.solver.history, t));
    end
else
   y = g_grind.solver.history;
end

function sol = dde_euler1(odefile, lags, history, tspan, options)
%EULER   Euler integration
%   Solve an ordinary differential equation using simple Euler integration.
%   Use this integration method only with very small
%   time steps or to study numerical problems with differential equations.
%
%   See also solver, rk4, heun, <a href="matlab:help ode45">ode45</a>, <a href="matlab:help ode23">ode23</a>, <a href="matlab:help ode113">ode113</a>, <a href="matlab:help ode15s">ode15s</a>, <a href="matlab:help ode23s">ode23s</a> 
%   <a href="matlab:help ode23t">ode23t</a>, <a href="matlab:help ode23tb">ode23tb</a>, <a href="matlab:help odeset">odeset</a>, <a href="matlab:help odeget">odeget</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('euler')">commands euler</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 22-Mar-2019 14:31:43 $
%function [tout, yout] = euler(odefunct, tspan, y0, options)
global g_grind;
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
else
    nonNegative=odeget(options,'NonNegative',[],'fast');
end
anyNonNegative=~isempty(nonNegative);
haveOutputFcn=~isempty(options.OutputFcn);
outputFcn=options.OutputFcn;
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
    error('GRIND:euler:tpan','The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t0);
if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan-1)) <= 0)
    error('GRIND:euler:tspan','The entries in tspan must strictly increase or decrease.');
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
        yout = nan(ntspan, neq);
    else
        tout = transpose(t0:delta:tfinal);
        if tout(end)<tfinal %if tfinal cannot divided in delta's
            tout(end+1)=tfinal;
        end
        
        yout = nan(size(tout, 1),neq);
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
while running
    f=feval(odefunct, t, y);
    %  simple Euler
    ynew = y + f .*  delta;
    if anyNonNegative
        ynew(nonNegative) = max(ynew(nonNegative),0);
    end
    
    
    tnew = t + delta;
    if tnew>=tfinal
        running=0;
    end
    
    if ~outflag                 % computed points, no refinement only the last value
        nout = nout + 1;
        if ~running
            if nout>length(tout)
                nout=length(tout);
            end
            t1=tout(nout);
            yout(nout, :)=transpose(y+(ynew-y)./(tnew-t).*(t1-t));
        else
            yout(nout, :) = transpose(ynew);
        end
        tnew=tout(nout);
        if haveOutputFcn
            if feval(outputFcn,tnew,transpose(ynew),'')
                ndx=~isnan(yout(:,1));
                yout=yout(ndx,:);
                tout=tout(ndx);
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
                ndx=~isnan(yout(:,1));
                yout=yout(ndx,:);
                tout=tout(ndx);
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


