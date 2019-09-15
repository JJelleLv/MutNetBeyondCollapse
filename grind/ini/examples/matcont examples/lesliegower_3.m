%Example MatContM
%Leslie-Gower model

%open the model
use('lesliegower.ini');

%create object "g_cont" that uses the matcontm engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[20;18;0.55;0.36;0.3;0.23;0.08;0.3;0.18;0.26;0.72;0.23;0.29;0.98];
g_cont.add_points('id','EP1','msg','stable node (++):  [11.17052012;11.94625068;20.99523834;15.21083594]','p0',g_set{1},...
  'x0',[11.17052012;11.94625068;20.99523834;15.21083594]);
g_cont.add_points('id','EP2','msg','saddle (+/0):  -0 -0 32.687 23.6814 ','p0',g_set{1},...
  'x0',[-1.827858793e-011;-1.531380443e-011;32.68698061;23.68138391]);
g_cont.add_points('id','EP3','msg','saddle (+/0):  21.5029 22.9961 0 0 ','p0',g_set{1},...
  'x0',[21.50285631;22.99611022;6.288387266e-012;3.937669331e-012]);
g_cont.add_points('id','EP4','msg','unstable node (0):  0 -0 0 0 ','p0',g_set{1},...
  'x0',[1.696888526e-011;-2.019649754e-013;1.36963757e-011;0]);
g_cont.add_points('id','P1','msg','initial point:  0.001 0.001 0.001 0.001 ','p0',g_set{1},...
  'x0',[0.001;0.001;0.001;0.001]);

g_cont.select_point('P1',true); %select the last parameter values in GRIND

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',900,...
    'par1','b1',...
    'parranges1',[10 40]);

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP1','FP');

%change the matcontm/GRIND settings
g_cont.set('par2','cjy',...
    'parranges2',[0 1.4]);

%Continue Period Doubling (m) to a PD - Period doubling curve curve
g_cont.cont('PD1','PD');

%change the matcontm/GRIND settings
g_cont.set('Backward',false,...
    'amp',0.001,...
    'MaxNumPoints',200);

%Continue Generalized Period Doubling to a LP2 - Limit point map curve
g_cont.cont('GPD1','LP2');

%Continue Generalized Period Doubling to a LP2 - Limit point map curve
g_cont.cont('GPD2','LP2');

%Continue Generalized Period Doubling to a LP2 - Limit point map curve
g_cont.cont('GPD2','LP2');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

