%%example of how to use basic MATCONT commands with a GRIND model
%use the model
use MLfast

%create a grind_matcont class

g_cont=grind_matcont;
%add the matcont paths
g_cont.open

%% here the MATCONT commands start
% note that @g_cont.handles is a MATCONT style handle to the model
% the grind command par('-vector') returns the current parameters
% val('-vector') returns the current state variables as vector
%
%fold1 example, see MatCont manual
p=par('-vector');ap1=[1];

[x0,v0]=init_EP_EP(@g_cont.handles,[0.047222;0.32564],p,ap1);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',65);
opt=contset(opt,'MinStepsize',0.00001);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'Backward',1);

[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);

x1=x(1:2,s(2).index);p=[x(end,s(2).index);0.1];
[x0,v0]=init_H_LC(@g_cont.handles,x1,p,ap1,0.0001,30,4);

opt=contset;
opt=contset(opt,'IgnoreSingularity',1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
[x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);

plotcycle(x2,v2,s2,[1 2]);
%% here the MATCONT commands stop
% remove matcont from the current path (to avoid problems)
g_cont.close;
