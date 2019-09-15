function ydot=i_polar2cart(t,ycart,varargin)
%transforms a model with polar coordinates to cartesian
%is vectorizable if the odehandle is vectorizable
%(used in null for instance and is very fast)
%At least 2 dimensions are requiered and coordinates of radius and theta
%should be supplied + handle to odefile
%Using g_grind is slow but complete binding to anonymus function is more complex
%to implement
global g_grind
x=ycart(g_grind.solver.polar.radius,:);
y=ycart(g_grind.solver.polar.theta,:);
theta = atan2(y,x);
radius =sqrt(x.^2+y.^2);
ycart(g_grind.solver.polar.radius,:)=radius;
ycart(g_grind.solver.polar.theta,:)=theta;
ypoldot=g_grind.solver.polar.odehandle(t,ycart,varargin{:});
ydot=ypoldot;
ydot(g_grind.solver.polar.radius,:)=x./radius.*ypoldot(g_grind.solver.polar.radius,:)-y.*ypoldot(g_grind.solver.polar.theta,:);
ydot(g_grind.solver.polar.theta,:)=y./radius.*ypoldot(g_grind.solver.polar.radius,:)+x.*ypoldot(g_grind.solver.polar.theta,:);


