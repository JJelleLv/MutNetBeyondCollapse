%Example MatCont
%autonomous electronic circuit
%Ref:
%Freire, E., Rodriguez-Luis, A., Gamero E. and Ponce, E., A case study for homo-
%clinic chaos in an autonomous electronic circuit: A trip form Takens-Bogdanov to
%Hopf- Shilnikov, Physica D 62 (1993) 230ÿ253.
%
%open the model
use('torBPC.ini');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[0.32858;0.93358;0.5;0.001;-0.6;-0.9;0.6];
g_cont.add_points('id','EP1','msg','unstable spiral (--): 0.5612 -0.001 0.4464 ','p0',g_set{1},...
  'x0',[0.5612161775;-0.001;0.4464138828]);
g_cont.add_points('id','EP2','msg','unstable spiral (+/0): 0.0013 -0.001 0.0005 ','p0',g_set{1},...
  'x0',[0.001250028186;-0.0009999999959;0.0005250247257]);
g_cont.add_points('id','EP3','msg','unstable spiral (--): -0.5647 -0.001 -0.4497 ','p0',g_set{1},...
  'x0',[-0.5646852112;-0.0009999999952;-0.4496523034]);
g_cont.add_points('id','P1','msg','initial point: 0.0013 -0.001 0.0005 ','p0',g_set{1},...
  'x0',[0.00125003;-0.001;0.000525025]);

%change the matcont/GRIND settings
g_cont.set('Userhandles',@(t,x,a3,b3,beta,epsilon,gamma,nu,r)epsilon,...
    'InitStepsize',[],...
    'MaxStepsize',[],...
    'MinStepsize',[],...
    'UserfunctionsInfo',struct('name','epsilon','state',1,'label','E0'),...
    'MaxNumPoints',10,...
    'par1','nu');

%Continue Equilibrium to a EP - Equilibrium points curve
g_cont.cont('EP2','EP');

%change the matcont/GRIND settings
g_cont.set('Backward',false,...
    'MaxNumPoints',50);

%Continue Hopf bifurcation to a LC - Limit cycle curve
g_cont.cont('H1','LC');


g_cont.set('par2','epsilon',...
    'NTST',25,...
    'NCOL',4,...
    'VarTolerance',1e-4,...
    'FunTolerance',1e-4,...
    'Backward',true,...
    'MaxNumPoints',15);

g_cont.cont('NS+1','NS');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

