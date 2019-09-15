%RK4   Runga-Kutta 4 solver with fixed time step
%   Integrates a system of ordinary differential equations
%   using the fourth order Runge-Kutta  method.  See also <a href="matlab:help ode45">ode45</a>
%   and <a href="matlab:help odedemo">odedemo</a>.
%   [t,y] = rk4('yprime', tspan, y0) integrates the system
%   of ordinary differential equations described by the M-file
%   yprime.m over the interval tspan=[t0,tfinal] and using
%   initial conditions y0.
%   [t, y] = rk4(F, tspan, y0, options) options.StepSize= step
%   size to be used
%
%   See also solver, c.rk4, euler, <a href="matlab:help ode45">ode45</a>, <a href="matlab:help ode23">ode23</a>, <a href="matlab:help ode113">ode113</a>, <a href="matlab:help ode15s">ode15s</a>, <a href="matlab:help ode23s">ode23s</a> 
%   <a href="matlab:help ode23t">ode23t</a>, <a href="matlab:help ode23tb">ode23tb</a>, <a href="matlab:help odeset">odeset</a>, <a href="matlab:help odeget">odeget</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('rk4')">commands rk4</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [tout, yout] = rk4(FunFcn, tspan, y0, options)

% INPUT:
% F     - String containing name of user-supplied problem description.
%         Call: yprime = fun(t,y) where F = 'fun'.
%         t      - Time (scalar).
%         y      - Solution column-vector.
%         yprime - Returned derivative column-vector; yprime(i) = dy(i)/dt.
% tspan = [t0, tfinal], where t0 is the initial value of t, and tfinal is
%         the final value of t.
% y0    - Initial value column-vector.
% ssize - The step size to be used. (Default: ssize = (tfinal - t0)/100).
%
% OUTPUT:
% t  - Returned integration time points (column-vector).
% y  - Returned solution, one solution column-vector per tout-value.
%
% The result can be displayed by: plot(t,y).


% Initialization

t0 = tspan(1);
tfinal = tspan(end);
pm = sign(tfinal - t0);  % Which way are we computing?
if nargin < 4 || isempty(options.StepSize)
   ssize = 0.3;   
else
   ssize = options.StepSize;
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

if ssize < 0, ssize = -ssize; end
%if ssize > 0.3, ssize = 0.3; end
h = pm * ssize;
t = t0;
y = y0(:);

% We need to compute the number of steps.

dt = abs(tfinal - t0);
N = floor(dt / ssize) + 1;
if (N - 1) * ssize < dt
   N = N + 1;
end

% Initialize the output.

tout = zeros(N, 1);
tout(1) = t;
yout = zeros(N, size(y, 1));
yout(1, :) = transpose(y);
k = 1;
running=true;
% The main loop
while (k < N)&&running
   if pm * (t + h - tfinal) > 0
      h = tfinal - t;
      tout(k + 1) = tfinal;
   else
      tout(k + 1) = t0  + k * h;
   end
   k = k + 1;
   % Compute the slopes
   s1 = feval(FunFcn, t, y);
   s2 = feval(FunFcn, t + h / 2, y + h * s1 / 2);
   s3 = feval(FunFcn, t + h / 2, y + h * s2 / 2); 
   s4 = feval(FunFcn, t + h, y + h * s3);
   y = y + h * (s1 + 2 * s2 + 2 * s3 +s4) / 6;
   t = tout(k);
   if anyNonNegative
       %truncate negative values to zero
       ndx = nonNegative( y(nonNegative) < 0 );
       y(ndx) = max(y(ndx),0);
    end
    if haveOutputFcn
        if feval(outputFcn,t,y,'')
          running = false;
          yout=yout(1:k,:);
          tout=tout(1:k,:);
        end  
    end
   yout(k, :) = transpose(y);
end

if length(tspan)>2 % not very efficient
   yout=interp1(tout,yout,tspan);
   tout=tspan;
end

