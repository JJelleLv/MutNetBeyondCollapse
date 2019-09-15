%example of MatCont
%open the model
use('cstr');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[3;0;0;0];
g_cont.add_points('id','EP1','msg','stable node (--): -0.9 ','p0',g_set{1},...
    'x0',-0.9);
g_cont.add_points('id','P1','msg','initial point: 0.001 ','p0',g_set{1},...
    'x0',0.001);

%change the matcont/GRIND settings
g_cont.set('FunTolerance',0.001,...
    'VarTolerance',0.001,...
    'MaxNumPoints',50,...
    'MaxStepsize',1,...
    'par1','lambda',...
    'symbolic',true,...
    'Backward',false);

%Continue Equilibrium to a EP - Equilibrium points curve
g_cont.cont('EP1','EP');

%change the matcont/GRIND settings
g_cont.set('BranchingPars',[1:4],...
    'MaxNumPoints',300,...
    'par2','beta');

%Continue Fold bifurcation (limit point) to a F - Fold bifurcation curve
g_cont.cont('F2','F');

%change the matcont/GRIND settings
g_cont.set('parranges1',[-1 3],'par3','gamma','Backward',true);

%Continue Branch Point of Fold curves to a BP - Branchpoint of fold bifurcation curve
g_cont.cont('BP1','BP');

g_cont.cont('BP2','BP');

g_cont.show;

g_cont.plot('curvendx',1);

g_cont.plot('fun',{'lambda','gamma','beta'},'curvendx',[2 3 4]);

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

