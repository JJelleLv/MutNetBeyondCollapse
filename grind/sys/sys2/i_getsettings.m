function [settings] = i_getsettings(N0,ndays)
global g_grind t;
if nargin<2
    ndays=g_grind.ndays;
end
if isempty(g_grind.tstep)
    g_grind.tstep=nan;
end


solverlist =  solver('-list'); 
solvername=solver('name');
if isempty(g_grind.solver.opt.StepSize)
    g_grind.solver.opt.StepSize=0.1;
end

settings.solver=[ndays;t;g_grind.tstep;g_grind.solver.iters;g_grind.solver.backwards;...
    g_grind.solver.addmode;g_grind.solver.opt.RelTol;g_grind.solver.opt.AbsTol;...
    g_grind.solver.opt.StepSize;find(strcmpi(solverlist, solvername));isempty(g_grind.solver.opt.Jacobian);g_grind.solver.opt.NonNegative(:)];
s=sprintf('%s(:);',g_grind.pars{:});
settings.pars=transpose(evalin('base',sprintf('[%s];',s)));
if nargin == 0
   N0 = i_initvar;
end

settings.initvar=N0;
s={};
for i=1:length(g_grind.permanent)
    s=[s {g_grind.permanent{i}.name}];
end

s=sprintf('%s(:);',s{:});
s=sprintf('[%s]',s(1:end-1));
settings.permanent=evalin('base',s);
% for i = 1:length(g_grind.pars)
%    p = evalin('base', char(g_grind.pars{i}));
%    jj=numel(p);
%    j2=j+jj;
%    while j2 + 1 > maxn
%       addn=j2 + 1;
%       maxn = maxn + addn;
%       settings = [settings; zeros(addn,1)];
%    end

%    settings(j+1:j2)=p(1:jj);
%    j=j2;
% end


% settings=[settings(1:j);N0];
% j=j+length(N0);
% if isfield(g_grind,'permanent')
%     for i=1:size(g_grind.permanent)
%      p = evalin('base', char(g_grind.permanent{i}.name));
%      jj=numel(p);
%      j2=j+jj;
%      while j2 + 1 > maxn
%         addn=j2 + 1;
%         maxn = maxn + addn;
%         settings = [settings; zeros(addn,1)];
%      end

%      settings(j+1:j2)=p(1:jj);
%      j=j2;
%     end

% end

% g_grind.parvalues=cell(size(g_grind.pars));
% for i=1:length(g_grind.pars)
%     g_grind.parvalues{i}=evalin('base',g_grind.pars{i});
% end


