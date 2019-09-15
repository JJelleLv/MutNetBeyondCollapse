%Example MatContM
%Leslie-Gower model

%open the model
use('lesliegower.ini');

%create object "g_cont" that uses the matcontm engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[20;18;0.55;0.36;0;0.23;0.08;0;0.18;0.26;0.72;0.23;0.29;0.98];
g_cont.add_points('id','EP1','msg','stable node (++): 16.4291 17.57 28.8712 20.9169 ','p0',g_set{1},...
  'x0',[16.42912109;17.57003227;28.87121737;20.91690238]);
g_cont.add_points('id','EP2','msg','saddle (+/0): -0 0 32.687 23.6814 ','p0',g_set{1},...
  'x0',[0;0;32.68698061;23.68138391]);
g_cont.add_points('id','EP3','msg','saddle (+/0): 21.5029 22.9961 0 0 ','p0',g_set{1},...
  'x0',[21.50285631;22.99611022;0;0]);
g_cont.add_points('id','EP4','msg','unstable node (0): 0 -0 0 -0 ','p0',g_set{1},...
  'x0',[9.52056708e-012;-1.014365058e-013;2.019505878e-011;-5.825033874e-015]);
g_cont.add_points('id','P1','msg','initial point: 0.001 0.001 0.001 0.001 ','p0',g_set{1},...
  'x0',[0.001;0.001;0.001;0.001]);

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',50,...
    'Backward',0,...
    'par1','cyj');

%Continue Equilibrium to a FP - Fixed point (maps) curve
g_cont.cont('EP3','FP');

%change the matcontm/GRIND settings
g_cont.set('MaxNumPoints',250,...
    'Backward',[]);

g_cont.cont('T1','FP');

g_cont.set('MaxNumPoints',100,...
    'Backward',1);

g_cont.cont('PD1','FP2');

g_cont.set('MaxNumPoints',200,...
    'Backward',[],...
    'par1','cjy');

g_cont.cont('EP2','FP');

g_cont.cont('PD3','FP2');

g_cont.set('MaxNumPoints',200,...
    'Backward',1)
g_cont.cont('T3','FP');

g_cont.set('MaxNumPoints',350,...
    'Backward',[])

g_cont.cont('EP1','FP');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

