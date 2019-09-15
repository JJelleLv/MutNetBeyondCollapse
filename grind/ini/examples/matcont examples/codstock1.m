% See CodStock1-5 examples MatcontM

%open the model
use('rfish.ini');

%create object "g_cont" that uses the matcontm engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points (can also be done as
%g_cont.add_points(findeqs))
g_set{1}=[1;120;0.5;0.9;0.5];
g_cont.add_points('id','EP1','msg','stable spiral (++)','p0',g_set{1},...
  'x0',[32.28471687;2.130466379]);
g_cont.add_points('id','EP2','msg','unstable node (0)','p0',g_set{1},...
  'x0',[-9.721803317e-013;1.04293643e-014]);
g_cont.add_points('id','P1','msg','initial point','p0',g_set{1},...
  'x0',[0.001;0.001]);

%change the matcontm/GRIND settings
g_cont.set('Backward',false,...
    'par1','F');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',1250,...
    'Backward',true);

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',300);

%Continue Transcritical (branch point) bifurcation to a FP - Fixed point (maps) curve
g_cont.cont('T1','FP');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',2750,...
    'par2','m2',...
    'Backward',false);

%Continue Neimark-Sacker(m) to a NS - Neimark sacker curve curve
g_cont.cont('NS1','NS');

%change the matcontm/GRIND settings
g_cont.set('Backward',true,...
    'MaxNumPoints',150);

%Continue Neimark-Sacker(m) to a NS - Neimark sacker curve curve
g_cont.cont('NS1','NS');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

