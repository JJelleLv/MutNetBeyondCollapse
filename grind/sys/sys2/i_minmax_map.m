function g_X2=i_minmax_map(~,g_X1,varargin)
global g_grind;

opt=g_grind.solver.map.opt;
if g_grind.solver.map.haslags
    opt.Events=@lag_minevents;
else
    opt.Events=@minevents;
end

g_grind.solver.map.h=i_getodehandle(7,'');
[~, Y] = feval(str2func(g_grind.solver.map.solver),g_grind.solver.map.h, [0,5000,10000], g_X1,opt);
g_X2=transpose(Y(end,:));
function [value,isterminal,direction] = minevents(x,y)
% we define two events.
% at minimum, dydx = 0, direction = 1 event function is increasing after
% the minimum.
%
% at maximum, dydx = 0, direction = -1 event function is decreasing after
% the maximum
%
% We create a vector for each output.
global g_grind;
%myode=str2func(g_grind.solver.map.odefile);
value = g_grind.solver.map.h(x,y); % this is the derivative
value=value(g_grind.solver.map.varno);
isterminal = 1;   % do not stop the integration
direction = g_grind.solver.map.dir;
function [value,isterminal,direction] = lag_minevents(x,y,lags)
% we define two events.
% at minimum, dydx = 0, direction = 1 event function is increasing after
% the minimum.
%
% at maximum, dydx = 0, direction = -1 event function is decreasing after
% the maximum
%
% We create a vector for each output.
global g_grind;
%myode=str2func(g_grind.solver.map.odefile);
value = g_grind.solver.map.h(x,y,lags); % this is the derivative
value=value(g_grind.solver.map.varno);
isterminal = 1;   % do not stop the integration
direction = g_grind.solver.map.dir;
