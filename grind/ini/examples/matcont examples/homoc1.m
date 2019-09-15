%Example MatCont: a modification of the fast subsystem of the Morris-Lecar equations
%the Morris-Lecar equations were introduced  as a model for the electrical activity in
%the barnacle giant muscle fiber. Example from matcont
%
%select the model
use mlfast

%create an instance of the class grind_matcont (g_cont facilitates the use
%of conteq later)
global g_cont;
g_cont=grind_matcont;

%select first parameter and set various Matcont settings
g_cont.set('par1','y',...
    'Singularities',1,...
    'MaxNumPoints',65,...
    'MinStepsize',0.00001,...
    'MaxStepsize',0.01,...
    'Backward',1);

%add an initial point (is at an equilibrium)
pnt=g_cont.add_points([0.047222;0.32564]);

%continue the equilibrium
g_cont.cont(pnt,'EP');

%change some settings and set the specific settings amp (initial amplitude)
%NTST number of mesh intervals and NCOL number of collocation points
g_cont.set('MaxNumPoints',200,...
    'MaxStepsize',1,...
    'IgnoreSingularity',1,...
    'amp',0.0001,...
    'NTST',30,...
    'NCOL',4);

%start form the subcritical Hopf bifurcation and continue the unstable
%limitcycle
g_cont.cont('H+1','LC');

%%set second parameter and specific settings for continuing Homoclinic
%bifurcation
g_cont.set('par2','z',...
    'MaxNumPoints',15,...
    'MaxStepsize',1,...
    'MinStepsize',[],...
    'InitStepsize',[],...
    'extravec',[0 1 1],...
    'Backward',0,...
    'eps0',0.01,...
    'eps1',0.01);

%continue homoclinic bifurcation
g_cont.cont('LC1','Hom');

%show the special points
g_cont.show

%plot all curves
g_cont.plot

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')
