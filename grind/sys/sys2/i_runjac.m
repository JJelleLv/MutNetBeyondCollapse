function res=i_runjac(t,g_A)
global g_grind;
%if t>g_grind.lyapspect.tJac
%   g_grind.lyapspect.tJac=t+g_grind.lyapspect.dtJac;
   g_grind.lyapspect.J=i_calcjac(g_grind.lyapspect.num,1,g_A(1:g_grind.lyapspect.d));
%end

res=g_A;
res(1:g_grind.lyapspect.d)=feval( g_grind.lyapspect.odehandle,t,g_A(1:g_grind.lyapspect.d));
res(g_grind.lyapspect.d+1:g_grind.lyapspect.d+g_grind.lyapspect.L)=g_grind.lyapspect.J*reshape(g_A(g_grind.lyapspect.d+...
   1:g_grind.lyapspect.d+g_grind.lyapspect.L),g_grind.lyapspect.d,g_grind.lyapspect.d);



