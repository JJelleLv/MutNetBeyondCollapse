% See CodStock5-11 examples MatcontM
%%Parameters:
%[55;0.5;0.5;0.5;3;  2; 1]
%(F, P,  m1, m2, B1,B2,B3)
%

%open the model
use('NAFstock.ini');

%create object "g_cont" that uses the matcontm engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[3;2;1;55;0.5;0.5;0.5];
g_cont.add_points('id','EP1','msg','unstable spiral (++): 2.762 1.008 ','p0',g_set{1},...
  'x0',[2.761950592;1.007984437]);
g_cont.add_points('id','EP2','msg','unstable node (0): 0 0 ','p0',g_set{1},...
  'x0',[2.598612855e-012;3.312720418e-014]);
g_cont.add_points('id','P1','msg','initial point: 0.001 0.001 ','p0',g_set{1},...
  'x0',[0.001;0.001]);

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',560,...
    'Backward',true,...
    'par1','F');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('Backward',false,...
    'MaxNumPoints',1560,...
    'par1','B1');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',50,...
    'par1','B2');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',500,...
    'par1','B3');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('par2','B3',...
    'MaxNumPoints',500,...
    'par1','F');

%Continue Neimark-Sacker(m) to a NS - Neimark sacker curve curve
g_cont.cont('NS1','NS');

%change the matcontm/GRIND settings
g_cont.set('amp',1e-005,...
    'Backward',true,...
    'MaxNumPoints',250);

%Continue 1:4 Resonance to a LP4m1 - Limit point map curve
g_cont.cont('R4_1','LP4m1');

%change the matcontm/GRIND settings
g_cont.set('Backward',false,...
    'amp',1e-010);

%Continue 1:4 Resonance to a LP4m2 - Limit point map curve
g_cont.cont('R4_1','LP4m2');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',400,...
    'par1','m2');

%Continue Neimark-Sacker(m) to a NS - Neimark sacker curve curve
g_cont.cont('NS1','NS');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

