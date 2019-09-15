function [tout,yout] = i_ode15sol(odefile,tspan,y0,options,varargin)
%yp0=zeros(size(y0));
global g_grind;
if g_grind.solver.isimplicit
   yp0 =rand(size(y0));
   [y0,yp0] = decic(odefile,tspan(1),y0,ones(size(y0)),yp0,zeros(size(y0)));
   [tout,yout] = ode15i(odefile,tspan,y0,yp0,options);
else
   %normal ode
   yp0 = feval(odefile,tspan(1),y0);
%  decic not needed as we have consistent yp0
   [tout,yout] = ode15i(i_getodehandle(3,''),tspan,y0,yp0,options);
end


