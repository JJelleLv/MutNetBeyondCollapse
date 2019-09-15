%HEUN   Heun integration
%   Solve an ordinary differential equation using improved Euler integration (Heun).
%   This simple integration method is much more accurate than Euler integration.
%
%   See also solver, rk4, stochast_heun, euler <a href="matlab:help ode45">ode45</a>, <a href="matlab:help ode23">ode23</a>, <a href="matlab:help ode113">ode113</a>, <a href="matlab:help ode15s">ode15s</a>, <a href="matlab:help ode23s">ode23s</a> 
%   <a href="matlab:help ode23t">ode23t</a>, <a href="matlab:help ode23tb">ode23tb</a>, <a href="matlab:help odeset">odeset</a>, <a href="matlab:help odeget">odeget</a>
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('heun')">commands heun</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [tout, yout] = heun(odefunct, tspan, y0, options)
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
while running
   f1=feval(odefunct, t, y);
   ypred = y + f1 .*  delta; %predictor  
   tnew = t + delta;
   f2=feval(odefunct, tnew, ypred); %corrector
   ynew = y + (f1+f2) .*  delta/2;
   if anyNonNegative
       ndx = nonNegative( ynew(nonNegative) < 0 );
       ynew(ndx) = max(ynew(ndx),0);
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

