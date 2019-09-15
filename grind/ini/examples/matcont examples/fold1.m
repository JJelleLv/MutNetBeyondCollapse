%fold1 example MatCont

%open the model
use('mlfast.ini');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%first add initial point
g_cont.add_points([0.047222;0.32564]);

%change the matcont/GRIND settings
g_cont.set('MinStepsize',1e-05,...
    'MaxStepsize',0.01,...
    'MaxNumPoints',65,...
    'par1','y');

%Continue Equilibrium to a EP - Equilibrium points curve
g_cont.cont('EP1','EP');

%change the matcont/GRIND settings
g_cont.set('par2','z',...
    'parranges2',[-0.1 0.3],...
    'MaxNumPoints',200);

%Continue Hopf bifurcation to a H - Hopf bifurcation curve
g_cont.cont('H+1','H');

%Expand the previous curve
g_cont.cont('H+1','H','expand');

%Continue Fold bifurcation (limit point) to a F - Fold bifurcation curve
g_cont.cont('F1','F');

%Continue Fold bifurcation (limit point) to a F - Fold bifurcation curve
g_cont.cont('F2','F');

%change the matcont/GRIND settings
g_cont.set('IgnoreSingularity',1,...
    'par2','',...
    'MaxNumPoints',50,...
    'MaxStepsize',1);

%Continue Hopf bifurcation to a LC - Limit cycle curve
g_cont.cont('H+1','LC');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

