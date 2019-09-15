%Example MatContM

%open the model
use('tnfmap.ini');

%create object "g_cont" that uses the matcontm engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[-1;0;1;1];
g_cont.add_points('id','EP1','msg','unstable node (++): 1.291 2.582 ','p0',g_set{1},...
  'x0',[1.290994449;2.581988897]);
g_cont.add_points('id','EP2','msg','unstable spiral (0): -0 0 ','p0',g_set{1},...
  'x0',[-3.06674568e-012;0]);
g_cont.add_points('id','EP3','msg','unstable node (--): -1.291 -2.582 ','p0',g_set{1},...
  'x0',[-1.290994449;-2.581988897]);
g_cont.add_points('id','P1','msg','initial point: 0.001 0.001 ','p0',g_set{1},...
  'x0',[0.001;0.001]);

%change the matcontm/GRIND settings
g_cont.set('par1','beta2');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP2','FP');

%change the matcontm/GRIND settings
g_cont.set('Backward',false,...
    'MaxNumPoints',50,...
    'Userfunctions',1,...
    'UserfunctionsInfo',struct('label',{'B2_','B3_'},'name',{'beta2=2','beta2=0.5'},'state',{1,1}),...
    'Userhandles',{@(t,g_X,beta1,beta2,C,D)beta2-2,@(t,g_X,beta1,beta2,C,D)beta2-0.5});

%Continue Transcritical (branch point) bifurcation to a FP - Fixed point (maps) curve
g_cont.cont('T1','FP');

%change the matcontm/GRIND settings
g_cont.set('Backward',true);

%Continue Transcritical (branch point) bifurcation to a FP - Fixed point (maps) curve
g_cont.cont('T1','FP');

%change the matcontm/GRIND settings
g_cont.set('Backward',[],...
    'MaxNumPoints',300,...
    'Userfunctions',0,...
    'UserfunctionsInfo',[],...
    'Userhandles',[]);

%Continue Transcritical (branch point) bifurcation to a FP - Fixed point (maps) curve
g_cont.cont('T1','FP');

%Continue Period Doubling (m) to a FP2 - Fixed point 2 (maps) curve
g_cont.cont('PD1','FP2');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

