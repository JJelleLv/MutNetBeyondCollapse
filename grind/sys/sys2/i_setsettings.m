function  i_setsettings(settings)
global g_grind t;
% settings.solver=[g_grind.ndays;t;g_grind.tstep;g_grind.solver.iters;g_grind.solver.backwards;...
%     g_grind.solver.addmode;g_grind.solver.opt.RelTol;g_grind.solver.opt.AbsTol;...
%     g_grind.solver.opt.StepSize;find(strcmpi(solverlist, solvername));g_grind.solver.opt.NonNegative(:);isempty(g_grind.solver.opt.Jacobian)];
% 
g_grind.ndays=settings.solver(1);
t=settings.solver(2);
g_grind.tstep=settings.solver(3);
g_grind.solver.iters=settings.solver(4);
if g_grind.solver.backwards~=settings.solver(5)
    g_grind.solver.backwards=settings.solver(5);
    if ~isempty(g_grind.solver.opt.Jacobian)
        g_grind.solver.opt.Jacobian=i_getodehandle('Jacobian');
    end

end

g_grind.solver.addmode=settings.solver(6);
if settings.solver(7)>0
    g_grind.solver.opt.RelTol=settings.solver(7);
end

if settings.solver(8)>0
   g_grind.solver.opt.AbsTol=settings.solver(8);
end

if settings.solver(9)>0
   g_grind.solver.opt.StepSize=settings.solver(9);
end

solverlist = solver('list');
if settings.solver(10)>0
   solver(solverlist{settings.solver(10)});
end

j = 0;
for i = 1:numel(g_grind.pars)
   p = evalin('base', char(g_grind.pars{i}));
   s=size(p);
   jj=prod(s);
   assignin('base',char(g_grind.pars{i}),reshape(settings.pars(j+1:j+jj),s(1),s(2)));
   j=j+jj;
end

%jj2=numel(NP);
i_keep(settings.initvar,settings.permanent);
%,settings(j+jj+1:j+jj+jj2));




