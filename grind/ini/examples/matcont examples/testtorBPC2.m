%Example MatCont

%open the model
use('torBPC.ini');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%initialize parameters
beta=0.5;
gamma=-0.6;
r=0.6;
a3=0.32858;
b3=0.93358;
nu=-0.9;
epsilon=0;
%add [0;0;0]

pnt=g_cont.add_points([0;0;0]);

g_cont.set('MaxNumPoints',150,...
    'par1','nu',...
    'par2','',...
    'InitStepsize',[],...
    'MaxStepsize',[],...
    'MinStepsize',[],...
    'Backward',0,...
    'Adapt',5);

g_cont.cont(pnt,'EP');

g_cont.cont('H1','LC');

g_cont.set('amp',1e-2);

g_cont.cont('BPC1','LC');

g_cont.set('MaxNumPoints',200,...
    'BranchingPars','epsilon',...
    'par1','nu',...
    'par2','epsilon',...
    'par3','r');

g_cont.cont('BPC1','BPC');

g_cont.show;

g_cont.plot;
disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

