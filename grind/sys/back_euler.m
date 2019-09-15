%BACK_EULER   Backwards Euler integration
%   Solve an ordinary differential equation using the implicit backwards Euler integration.
%   This method is much more efficient for stiff equations, however it
%   cannot be used in stochastic models. The Jacobian is evaluated
%   numerically.
%
%   See also solver, rk4, <a href="matlab:help ode45">ode45</a>, <a href="matlab:help ode23">ode23</a>, <a href="matlab:help ode113">ode113</a>, <a href="matlab:help ode15s">ode15s</a>, <a href="matlab:help ode23s">ode23s</a> 
%   <a href="matlab:help ode23t">ode23t</a>, <a href="matlab:help ode23tb">ode23tb</a>, <a href="matlab:help odeset">odeset</a>, <a href="matlab:help odeget">odeget</a>

%   Copyright 2014 WUR
%   Revision: 1.1.8 $ $Date: 12-Mar-2014 09:06:31 $

function [tout, yout] = back_euler(odefunct, tspan, y0, options)
Jfun=[];
if (nargin < 4)
   %if there are no valid options use default stepsize
   delta = 0.1;
else
   %use the option StepSize as step size
   if ~isfield(options,'StepSize')
      options.StepSize=options.MaxStep;
   end
   delta = options.StepSize;
   if isempty(delta)
       delta=0.1;
   end
end
if nargin==4 && ~isempty(options.Jacobian)
   Jfun=options.Jacobian;
end
if nargin>=4
   nonNegative = ~isempty(odeget(options,'NonNegative',[],'fast'));
   if nonNegative
      warning('GRIND:back_euler:NonNegativeIgnored','BACK_EULER does not constrain solution to be non-negative. Option ''NonNegative'' will be ignored.');   
   end
end
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

% Set the output flag.

outflag = ntspan > 2;                          % output only at tspan points

% Allocate memory if we're generating output.
delta = delta * tdir;

if nargout > 0
   if ntspan > 2                         % output only at tspan points
      tout = zeros(ntspan, 1);
      yout = zeros(ntspan, neq);
%      nsteps=round(tdir*(tfinal-t0)/delta)+1;
   else
      tout = transpose(t0:delta:tfinal);
%      nsteps=size(tout,1);
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
   f=feval(odefunct, t, y);
   %simple euler as first guess
   ynew = y + f .*  delta;
   %optimize the result
   tnew = t + delta; 
   %need fast newton-raphson method
   ynew=newton4euler(odefunct,tnew,y,ynew,delta,Jfun);
   if tnew+0.01*delta>=tfinal
      running=0;
      tnew=tfinal;
   end
   if ~outflag                 % computed points, no refinement
      nout = nout + 1;
      tout(nout) = tnew;
      yout(nout, :) = transpose(ynew);
      if haveOutputFcn
        if feval(outputFcn,tnew,transpose(ynew),'')
          running = false;
        end  
      end 
   elseif (tdir * (tnew - tspan(next)) >= 0) % at tspan, tspan assumed to be larger than delta
       nout = nout + 1;
       tout(nout) = tnew;
       yout(nout, :) = transpose(ynew);
       next = next + 1;
       if haveOutputFcn
        if feval(outputFcn,tnew,transpose(ynew),'')
          running = false;
        end  
      end
   end
   y = ynew;
   t = tnew;
end   
if nout<length(tout)
   tout=tout(1:nout);
   yout=yout(1:nout,:);
end

function [Y,isConverged]=newton4euler(f,t,yold,Y,h,Jfun)
% [Y,isConverged]=newton4euler(f,x,ytranspose,Y,h)
% special function to evaluate Newton's method for back_euler
TOL = 1.e-8;
MAXITS = 100;
 
isConverged= false;  % starts out FALSE
for n=1:MAXITS
  fValue = feval( f, t, Y);
  if isempty(Jfun)
     fPartial = i_calcjac(1, 1, Y);
  else
     fPartial=Jfun(t,Y);
  end
  F = yold + h * fValue - Y;
  dFdY = h * fPartial - eye(length(Y));
  %solve the linear equations
  increment=dFdY\F;
  Y = Y - increment;
  if norm(increment,inf) < TOL*norm(Y,inf)
    isConverged=true;  % turns TRUE here
    return
  end
end


