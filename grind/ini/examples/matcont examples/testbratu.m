%Example MatCont

%open the model
use('bratu.ini');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=0;
g_cont.add_points('id','EP1','msg','stable node (0): -0 0 ','p0',g_set{1},...
  'x0',[-1.647654073e-012;0]);
g_cont.add_points('id','P1','msg','initial point: 0.001 0.001 ','p0',g_set{1},...
  'x0',[0.001;0.001]);

g_cont.select_point('P1',true); %select the last parameter values in GRIND

%change the matcont/GRIND settings
g_cont.set('par1','a',...
    'MaxNumPoints',50,...
    'Backward',false,...
    'Singularities',1,...
    'Userfunctions',1,...
    'UserfunctionsInfo',struct('name','userf1','state',1,'label','u1'),...
    'Userhandles',@(t,x,a)a-0.2)

%Continue Equilibrium to a EP - Equilibrium points curve
g_cont.cont('EP1','EP');

g_cont.cont('EP1','EP','expand');
g_cont.set('Backward',[]);
g_cont.cont('T1','EP')
g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

