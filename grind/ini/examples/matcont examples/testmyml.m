%Example MatCont

%open the model
use('myml');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%first search all equilibria automatically (this model has only one)
g_cont.add_points(findeqs);

%change the matcont/GRIND settings
g_cont.set('MaxStepsize',10,...
    'MaxNumPoints',2500,...
    'par1','Inp');

%Continue Equilibrium to a EP - Equilibrium points curve
g_cont.cont('EP1','EP');

%change the matcont/GRIND settings
g_cont.set('Backward',false,...
    'Multipliers',false,...
    'MaxNumPoints',500,...
    'MaxStepsize',5);

%Continue Hopf bifurcation to a LC - Limit cycle curve
g_cont.cont('H+1','LC');

%change the matcont/GRIND settings
g_cont.set('Adapt',0,...
    'par2','V3',...
    'MaxNumPoints',20,...
    'MaxStepsize',1);

%Continue Limit cycle to a Hom - Homoclinic bifurcation curve
g_cont.cont('LC1','Hom');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')


