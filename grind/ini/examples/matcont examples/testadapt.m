%Example MatCont: a modification of the fast subsystem of the Morris-Lecar equations
%the Morris-Lecar equations were introduced  as a model for the electrical activity in
%the barnacle giant muscle fiber. Example from matcont
%
use m_adapt

global g_cont;

%if the previous results have not saved to a mat file
if ~exist('test.mat','file')
    %create object "g_cont" that uses the matcont engine
    g_cont=grind_matcont;
    
    %first add  initial point
    g_cont.add_points([0;0;0]);
    
    
    %change the matcont/GRIND settings
    g_cont.set('par1','alpha');
    
    %Continue Equilibrium to a EP - Equilibrium points curve
    g_cont.cont('EP1','EP');
    
    %change the matcont/GRIND settings
    g_cont.set('Backward',false,...
        'MaxNumPoints',200,...
        'NTST',20,...
        'parranges1',[-2 2]);
    
    %Continue Hopf bifurcation to a LC - Limit cycle curve
    g_cont.cont('H1','LC');
    
    %change the matcont/GRIND settings
    g_cont.set('NTST',40,...
        'MaxNumPoints',250);
    
    %Continue Period Doubling (flip) bifurcation to a LC - Limit cycle curve
    g_cont.cont('PD+1','LC');
    
    %change the matcont/GRIND settings
    g_cont.set('par2','beta');
    
    %Continue Period Doubling (flip) bifurcation to a  curve
    g_cont.cont('PD+1','PD');
    
    %save a lengthy session to a mat file to save the figures and run
    g_cont.saveas('test.mat');
end

%load previously saved session
g_cont=grind_cont.load_session('test.mat');

g_cont.show;

g_cont.plot;

disp('You can use <a href="matlab:conteq">conteq</a> to continue working with this session')
