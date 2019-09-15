%scenario is a struc:
%required fields:
%scenario.pars = list of parameters to change
%scenario.parvalues = for each parameter the value in each step
%scenario.nstabilizing = number of steps for stabilizing
%scenario.nwriting = number of steps for output, or [nsteps and noutput]
%
%optional fields
%scenario.initialstabil(0)= number of steps for initial stabilizing
%scenario.stepcallback ([]) = callback function after each step
%scenario.reset (false) = reset each step to the current initial conditions. can also be a vector of the same
%size as scenario.parvarlues
%scenario.avoidextiction (true) = little input to avoid extinction?
function [data] = i_stepscenario(varargin)
global g_grind g_Y g_t t;
if nargin==1
    scenario=varargin{1};
else
    scenario=struct(varargin{:});
end

if ~isfield(scenario, 'avoidextiction')
   scenario.avoidextiction = true;
end

if ~isfield(scenario, 'reset')
   scenario.reset = false;
end

if ~isfield(scenario, 'stepcallback')
   scenario.stepcallback = [];
end

if ischar(scenario.pars)
    scenario.pars={scenario.pars};
end

data.pars = scenario.pars;
data.parvalues = scenario.parvalues;
oldpars=cell(length(data.pars));
for i1 = 1:length(data.pars)
   oldpars{i1} =  evalin('base', data.pars{i1});
end

data.perm = [];

%try
%   oldi = 1;
%  newi = 0;

t1 = t;
N0 = i_initvar;
OldN0 = N0;
oldtstep = g_grind.tstep;
hasperm = ~isempty(g_grind.permanent);
if ~g_grind.solver.isdiffer
   stabilsteps = 2;
else
   stabilsteps = NaN;
end

writesteps = oldtstep;
if length(scenario.nwriting) == 2
   writesteps = scenario.nwriting(2);
   scenario.nwriting = scenario.nwriting(1);
end

if isnan(writesteps)
    writesteps = round(scenario.nwriting);
end

if scenario.nwriting==0
    writesteps=0;
elseif writesteps<2
    writesteps=2;
end


nstat = size(N0, 1);
nsteps= size(scenario.parvalues,1);
   if writesteps==0
      data.Y = zeros(1, g_grind.statevars.dim, nsteps);
      data.t = zeros(1,1, nsteps);
   else
       data.Y = zeros(writesteps + 1, nstat, nsteps);
       data.t = zeros(writesteps + 1,1, nsteps);
   end

   if hasperm
      time(1, '-s')
      data.perm = zeros(size(data.Y,1), g_grind.permanent{end}.to,nsteps);
   else
      data.perm = [];
   end

if  isequal(scenario.stepcallback,@i_waitbar)
   i_waitbar(0, nsteps, 'grind', 'Calculating',0.5)
end

initstab=isfield(scenario, 'initialstabil')&&(scenario.initialstabil > 0);
for i = 1:nsteps
   %  waitbar(i/nsteps,wb);
   if length(scenario.reset)==1&&scenario.reset||length(scenario.reset)==nsteps&&scenario.reset(nsteps)
       N0=OldN0;
   end

   g_grind.solver.opt.OutputFcn = [];
   for i1 = 1:length(data.pars)
      multassignin('base', data.pars{i1},data.parvalues(i,i1));
   end

   if i==1&&initstab
      g_grind.tstep = stabilsteps;
      i_ru(t1, scenario.initialstabil, N0, 0);
      N0 = transpose(g_Y(size(g_Y, 1), :));
   end

   if isfield(scenario, 'nstabilizing')&&(scenario.nstabilizing > 0)
      g_grind.tstep = stabilsteps;
      i_ru(t1, scenario.nstabilizing, N0, 0);
      N0 = transpose(g_Y(size(g_Y, 1), :));
      drawnow;
   else
       g_t=[];
   end
   g_grind.tstep = writesteps;
   if scenario.nwriting == 0
      g_t = g_t(end);
      t1=g_t;
      g_Y = transpose(N0);
   else
      if ~isempty(g_t)
         t1 = g_t(size(g_t, 1));
      else
         t1 = 0;
      end
      i_ru(t1, scenario.nwriting, N0, 0);
      t1 = g_t(size(g_t, 1));
   end

   N0 = transpose(g_Y(size(g_Y, 1), :));
   if any(isnan(N0))
       warning('grind:paranal','Some values have become NaN(Not-a-Number) in step %d, replaced NaN with original initial values in next step',i);
       N0(isnan(N0))=OldN0(isnan(N0));
   end

   %little invasion to avoid hanging in a trivial equilibrium (tricky
   %if the model can have negative numbers
   if scenario.avoidextiction&&~g_grind.solver.isinteger
      for j = 1:nstat
         epsil = 0.001;
         if ((N0(j) < epsil) && (N0(j) > 0))
            N0(j) = epsil;
         elseif (N0(j) >  - epsil)  && (N0(j) < 0)
            N0(j) =     - epsil;
         elseif N0(j) == 0
            N0(j) = sign(OldN0(j)) * epsil;
            else %very small positive disturbance
            N0(j) = N0(j) + epsil * 0.001 * sign(N0(j));
         end
      end
   end

   if size(g_Y,1)<size(data.Y,1)
       warning('grind:paranal','Failed run at step %d, filled with NaN values',i);
       yy=nan(size(data.Y,1),size(data.Y,2));
       yy(1:size(g_Y,1),:)=g_Y;
       g_Y=yy;
       tt=nan(size(data.Y,1),1);
       tt(1:size(g_t,1),:)=g_t;
       g_t=tt;
   end

   data.Y(:,:,i) = g_Y;
   data.t(:,:,i) = g_t;
   if hasperm
      pperm = defpermanent('-g', []);
      data.perm(:,:,i) =  pperm;
   end

   if isa(scenario.stepcallback,'function_handle')
      if isequal(scenario.stepcallback,@i_waitbar)
        i_waitbar(1);
      elseif scenario.stepcallback(i,scenario)
          break;
      end

   end

   if isfield(scenario,'perturbcallback')&&isa(scenario.perturbcallback,'function_handle')
        N0=scenario.perturbcallback(N0,i,scenario,data.parvalues(i,:));
   end

   drawnow;
end

for i1 = 1:length(data.pars)
   multassignin('base', data.pars{i1},oldpars{i1});
end

g_grind.tstep=oldtstep;
if isequal(scenario.stepcallback,@i_waitbar)
    i_waitbar([]);
end

i_keep(OldN0);
