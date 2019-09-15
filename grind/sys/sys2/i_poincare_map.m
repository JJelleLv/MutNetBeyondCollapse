function g_X2=i_poincare_map(~,g_X1,varargin)
%g_grind.solver.map.plane = formula of a plane
%g_grind.solver.map.tperiod = period of simulation (keep small for speed,
%but you get error if too small).
%npoints = number of points for output
%increasing/decreasing not possible to change yet
global g_grind;
npoints=1000;
opt=g_grind.solver.map.opt;
ts=linspace(0,g_grind.solver.map.tperiod,npoints);
[ts, Y] = feval(str2func(g_grind.solver.map.solver), i_getodehandle(7,''), ts, g_X1,opt);
g_X2=findfirstpoincare(ts,Y);
if isempty(g_X2)
    error('grind:makemap','poincare plane not found: period for simulation may be too small');
end


function [y_poincare] = findfirstpoincare(t,Y)
%determines the interpolated y values where the plane is crossed
global g_grind;
data.pars={};
data.parvalues=[];
data.perm=[];
data.Y=zeros(size(Y,1),size(Y,2),1);
data.Y(:,:,1)=Y;
data.t=t;
y_poincare=[];
yy = outfun(g_grind.solver.map.plane,data);
crossing=yy(2:end-3,:)<0&yy(3:end-2,:)>=0; %increasing trajectory
%crossing=yy(2:end-3,:)>0&yy(3:end-2,:)<=0; %decreasing
i=find(crossing,1);
if ~isempty(i)
    i=i+3;
    y_poincare = interp1(yy(i - 2:i + 2), Y(i - 2:i + 2,:), 0);
end



