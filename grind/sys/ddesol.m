function [tout, yout, varargout] = ddesol(odefile, tspan, y0, options, varargin)
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
sol = dde23(odefile, lags, @ddehist, [tspan(1), tspan(length(tspan))], options);
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
