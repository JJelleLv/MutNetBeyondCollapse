%%Example MatCont: Model of Vandermeer 2004 describing a coupled predator-prey system
%   Bifurcations include: T BPC CH GH H HH LPC NSNS R1-4
%
%
%Vandermeer, J. (2004) "Coupled oscillations in food webs: Balancing competition 
%   and mutualism in simple ecological models" American Naturalist 163:857-867. 
%   http://dx.doi.org/10.1086/420776
%
%open the model
use('vandermeer.ini');

%create object "g_cont" that uses the matcont engine
global g_cont;
g_cont=grind_matcont;

%first add all initial points
g_set{1}=[2;0.8;1.3;0.18;1;0.1;1;1];
g_set{2}=[2;0.8;1.3;0.18;1.021904264;0.1;0.2670567244;1];
g_cont.add_points('id','EP1','msg','unstable spiral (++)','p0',g_set{1},...
  'x0',[0.4162179546;0.4162179546;0.04531858968;0.04531858969]);
g_cont.add_points('id','EP10','msg','saddle (0)','p0',g_set{2},...
  'x0',[0;0;0;0]);
g_cont.add_points('id','EP11','msg','saddle (--)','p0',g_set{2},...
  'x0',[0.03753712709;0;-0.1496994339;1.128752054]);
g_cont.add_points('id','EP2','msg','unstable spiral (++)','p0',g_set{2},...
  'x0',[0.04427262855;0.4841031437;0.0453185897;0.04531858969]);
g_cont.add_points('id','EP3','msg','unstable spiral (+/0)','p0',g_set{2},...
  'x0',[2.107188017;0;0;0.2970885324]);
g_cont.add_points('id','EP4','msg','saddle (+/0)','p0',g_set{2},...
  'x0',[0;0;0.5677245913;0.5677245912]);
g_cont.add_points('id','EP5','msg','saddle (+/0)','p0',g_set{2},...
  'x0',[0;0;1.021904264;0]);
g_cont.add_points('id','EP6','msg','saddle (+/0)','p0',g_set{2},...
  'x0',[0;0;0;1.021904264]);
g_cont.add_points('id','EP7','msg','unstable spiral (+/0)','p0',g_set{2},...
  'x0',[0;0.5627387296;0.2970885324;0]);
g_cont.add_points('id','EP8','msg','unstable spiral (+/0)','p0',g_set{2},...
  'x0',[0;0.5067755656;0;0.05347593583]);
g_cont.add_points('id','EP9','msg','unstable spiral (+/0)','p0',g_set{2},...
  'x0',[0.1353378225;0;0.05347593583;0]);
g_cont.add_points('id','P1','msg','initial point','p0',g_set{2},...
  'x0',[0.04427262854;0.4841031437;0.04531858969;0.04531858969]);

g_cont.select_point('P1',true); %select the last parameter values in GRIND

%change the matcont/GRIND settings
g_cont.set('par1','alpha','par2','beta');

%Continue Equilibrium to a EP - Equilibrium points curve
g_cont.cont('EP1','EP');



%Continue Hopf to a H - Hopf bifurcation curve
g_cont.cont('H1','H');

%change the matcont/GRIND settings
g_cont.set('stateranges',[0 NaN],...
    'parranges1',[0 2]);

g_cont.set('Backward',false,...
    'StepAmp',10,...
    'MaxNumPoints',50);

%Continue Hopf to a LC - Limit cycle curve (slow)
g_cont.cont('H1','LC');   

g_cont.set('par2','beta')

%Continue Hopf-Hopf to a NS1 - Neimark-Sacker curve (very slow)
g_cont.cont('HH1','NS1');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')

