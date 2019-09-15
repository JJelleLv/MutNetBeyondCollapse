%FOKKERPLANCK   Analyse Fokker-Planck equation based on 1D differential equation
%   By using a Fokker-Planck equation you can analyse the time-evolution of 
%   the probability density function of a stochastic differential equation. 
%   This command only works for one-equation differential equations.
%   Assume the following stochastic differential equation:
%   dX = F(X) dt + sigma dW
%   where X is a state variable, F(X) is the deterministic part and dW is the 
%   Wiener process (Gaussian white noise) times the function sigma.
%   The Fokker-Planck PDE is then defined as:
%   dP/dt=- d(F(X) P)/dX + sigma^2/2*dP^2/d^2X
%   P is the probability density function.
%   We solve this PDE using MATLAB's pdepe function with absorbing or reflecting
%   boundaries. Use simtime to define the simulation time, the initial 
%   condition of the original model is used(g_fokkerplanck.p0 can be used 
%   to define other than default initial conditions)
%      
%   Usage:
%   FOKKERPLANCK - opens a dialog box in which you can enter all necessary 
%   information.
%   FOKKERPLANCK SIGMA - sigma is defined, for the rest use default values.
%   FOKKERPLANCK('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'bl' [absorbing | reflecting] - kind of left boundary (default reflecting)
%     'br' [absorbing | reflecting] - kind of right boundary (default reflecting)
%     'diffsigma' [number] - derivative of sigma to x (leave blanc for numerical solving) (default '')
%     'lim' [number and length(number)==2] - the limits of the statevariable (default defined in ax or [0 10])
%     'npoints' [integer>0] - the number of points for the probability density function (default=300)
%     'sigma' [number] - the standard deviation of the noise (default 1)
%   FOKKERPLANCK('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-p' - plot: redraw the figures of the last run.
%
%  
%   See also potential, plotdiff, <a href="matlab:help pdepe">pdepe</a>  
%
%   Reference page in Help browser:
%      <a href="matlab:commands('fokkerplanck')">commands fokkerplanck</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function fokkerplanck(varargin)
evalin('base','global g_fokkerplanck');
global g_grind g_fokkerplanck;
if isfield(g_grind,'xaxis')
    theax=g_grind.xaxis.lim;
else
    theax=[0 10];
end
fieldnams={'sigma', 'n', 'the standard deviation of the noise (default 1)',1;...
   'diffsigma', 'n', 'derivative of sigma to x (leave blanc for numerical solving) (default '''')','';...
   'npoints', 'i>0', 'the number of points for the probability density function (default=300)',300;...
   'lim', 'n&length(n)==2', 'the limits of the statevariable (default defined in ax or [0 10])',theax;...
   'bl', 'e[a+bsorbing|r+eflecting]', 'kind of left boundary (default reflecting)','reflecting';...
   'br', 'e[a+bsorbing|r+eflecting]', 'kind of right boundary (default reflecting)','reflecting'}';
args=i_parseargs(fieldnams,'sigma','-p',varargin);
if g_grind.statevars.dim > 1||g_grind.solver.isdiffer
   error('grind:fokkerplanck:nodiffer','This command works only for 1D differential equations');
end
if g_grind.solver.nonautonomous||g_grind.solver.isstochastic
   error('grind:fokkerplanck:nonautonomous','This command works only for autonomous deterministic 1D differential equations');
end

if isfield(args,'sigma')
    args.sigma=i_checkstr(args.sigma);
end
if any(strcmp(args.opts,'-p'))
      makefigures;
      return;
end
g_grind.hfun.curr = i_getodehandle(1,'');
avar = i_statevars_names(1);
if isempty(g_fokkerplanck)
   g_grind.fokkerplanck.npoints = 300;
   if strcmp(g_grind.xaxis.var, avar)
      g_grind.fokkerplanck.lim = g_grind.xaxis.lim;
   elseif strcmp(g_grind.yaxis.var, avar)
      g_grind.fokkerplanck.lim = g_grind.yaxis.lim;
   elseif strcmp(g_grind.zaxis.var, avar)
      g_grind.fokkerplanck.lim = g_grind.zaxis.lim;
   else
      g_grind.fokkerplanck.lim = [0 10];
   end
   g_grind.fokkerplanck.bl = 0;
   g_grind.fokkerplanck.br = 0;
   g_grind.fokkerplanck.sigma = 1;
   g_fokkerplanck.diffsigma = [];
end
if nargin > 0
   f = fieldnames(g_grind.fokkerplanck);
   f1 = fieldnames(args);
   for i = 1:length(f)
      if any(strcmpi(f{i}, f1))
         g_grind.fokkerplanck.(f{i}) = args.(f{i});
      end
   end
else
   prompt={'Noise term sigma (can be an equation)','Derivative of sigma to x (leave blanc for numerical solution)','Range for state variable',...
      'What kind left boundary (reflecting/absorbing):','What kind right boundary (reflecting/absorbing):',...
      'Number of points in x'};
   name = 'Input for Fokker-Planck analysis';
   numlines = 1;
   defaultanswer = cell(1, 6);
   if ischar(g_grind.fokkerplanck.sigma)
      defaultanswer{1} = g_grind.fokkerplanck.sigma;
   else
      defaultanswer{1} = num2str(g_grind.fokkerplanck.sigma);
   end
   defaultanswer{2}  = '';
   defaultanswer{3} = mat2str(g_grind.fokkerplanck.lim);
   if g_grind.fokkerplanck.bl == 1
      defaultanswer{4} = 'absorbing';
   else
      defaultanswer{4} = 'reflecting';
   end
   if g_grind.fokkerplanck.br == 1
      defaultanswer{5} = 'absorbing';
   else
      defaultanswer{5} = 'reflecting';
   end
   defaultanswer{6} = num2str(g_grind.fokkerplanck.npoints);
   answer = inputdlg(prompt, name, numlines, defaultanswer);
   if isempty(answer)
      error('grind:fokkerplanck','Cancelled');
   end
   g_grind.fokkerplanck.sigma = answer{1};
   g_grind.fokkerplanck.diffsigma = answer{2};
   g_grind.fokkerplanck.lim = str2num(answer{3});  %#ok<ST2NM>
   g_grind.fokkerplanck.bl = answer{4};
   g_grind.fokkerplanck.br = answer{5};
   g_grind.fokkerplanck.npoints = str2double(answer{6});
end
update_g_grind;
x = linspace(g_grind.fokkerplanck.lim(1), g_grind.fokkerplanck.lim(2), g_grind.fokkerplanck.npoints);
if isnan(g_grind.tstep)
   t = linspace(0, g_grind.ndays, 1000);
else
   t = linspace(0, g_grind.ndays, g_grind.tstep);
end
%options=odeset('RelTol',1e-7,'AbsTol',1e-4);
g_fokkerplanck.t = t;
if ~isfield(g_fokkerplanck, 'p0')||length(g_fokkerplanck.p0)<=1
   g_fokkerplanck.p0  = i_initvar;
else
   g_fokkerplanck.p0 = interp1(g_fokkerplanck.X, g_fokkerplanck.p0, x);
   g_fokkerplanck.p0 = g_fokkerplanck.p0 / trapz(x, g_fokkerplanck.p0); %normalize to 1
end
g_fokkerplanck.X = x;
g_fokkerplanck.Ys = pdepe(0, @pdex1pde, @pdex1ic, @pdex1bc, x, t);
[~, aY] = plotdiff(avar, [g_fokkerplanck.X(1) g_fokkerplanck.X(end)], length(g_fokkerplanck.X));
g_fokkerplanck.isaddles=find(diff(sign(aY)) == 2);
if isempty(g_fokkerplanck.isaddles)||g_fokkerplanck.isaddles(1) > 1
   g_fokkerplanck.isaddles = [1; g_fokkerplanck.isaddles];
end
if g_fokkerplanck.isaddles(end) < length(g_fokkerplanck.X)
   g_fokkerplanck.isaddles = [g_fokkerplanck.isaddles ; length(g_fokkerplanck.X)];
end


% Boundary value problem (not analysing time)
%mean exit probability
if ~isempty(g_fokkerplanck.isaddles)&&(~g_grind.fokkerplanck.bl||~g_grind.fokkerplanck.br)
   % (we don't have an exit if there are no saddles and all boundaries are
   % reflecting)
   isabs = ones(size(g_fokkerplanck.isaddles));
   isabs(1) = g_grind.fokkerplanck.bl;
   isabs(end) = g_grind.fokkerplanck.br;
   g_fokkerplanck.uxl = cell(1, length(g_fokkerplanck.isaddles) - 1);
   g_fokkerplanck.uxr = g_fokkerplanck.uxl;
   g_fokkerplanck.Txl = g_fokkerplanck.uxl;
   g_fokkerplanck.Txb = g_fokkerplanck.uxl;
   g_fokkerplanck.Txr = g_fokkerplanck.uxl;
   g_fokkerplanck.xs = g_fokkerplanck.uxl;
   for j = 1:length(g_fokkerplanck.isaddles) - 1
      xran =  linspace(g_fokkerplanck.X(g_fokkerplanck.isaddles(j)), g_fokkerplanck.X(g_fokkerplanck.isaddles(j + 1)), g_grind.fokkerplanck.npoints);
      leftabs = isabs(j);
      rightabs = isabs(j + 1);
      g_fokkerplanck.xs{j} = xran;
      if leftabs&&rightabs
         g_fokkerplanck.uxl{j} = i_exitprob(xran, true);
         if numel(g_fokkerplanck.uxl{j}) > 1
            g_fokkerplanck.Txl{j} = i_meanexittime(xran, leftabs, rightabs, g_fokkerplanck.uxl{j});
            g_fokkerplanck.uxr{j} = 1 - g_fokkerplanck.uxl{j};
            g_fokkerplanck.Txr{j} = i_meanexittime(xran, leftabs, rightabs, g_fokkerplanck.uxr{j});
            g_fokkerplanck.Txb{j} = i_meanexittime(xran, leftabs, rightabs, 1);
         else
            g_fokkerplanck.Txl{j} = NaN;
            g_fokkerplanck.Txb{j} = NaN;
            g_fokkerplanck.Txr{j} = NaN;
            g_fokkerplanck.uxr{j} = NaN;
         end
      elseif leftabs
         g_fokkerplanck.uxl{j} = ones(size(xran));
         g_fokkerplanck.Txl{j} = i_meanexittime(xran, leftabs, rightabs, 1);
         g_fokkerplanck.Txr{j} = NaN;
         g_fokkerplanck.Txb{j} = NaN;
         g_fokkerplanck.uxr{j} = NaN;
      else
         g_fokkerplanck.uxr{j} = ones(size(xran));
         g_fokkerplanck.Txr{j} = i_meanexittime(xran, leftabs, rightabs, 1);
         g_fokkerplanck.Txl{j} = NaN;
         g_fokkerplanck.Txb{j} = NaN;
         g_fokkerplanck.uxl{j} = NaN;
      end
   end
end

makefigures;


% --------------------------------------------------------------
function u0 = pdex1ic(x)
if length(g_fokkerplanck.p0) == 1
   D_x = g_fokkerplanck.X(2) - g_fokkerplanck.X(1);
   [~, ii] = min(abs(g_fokkerplanck.X - g_fokkerplanck.p0));
   p0 = g_fokkerplanck.X(ii);
   if abs(x - p0) < D_x / 2||(abs(x - p0) == D_x / 2&&x > p0)
      u0 = 1 / D_x; %the integral should always be 1 if p0 is in the range of x
      if (x==g_fokkerplanck.X(1))||(x==g_fokkerplanck.X(end))
         %if p0 is at the borders we need to multiply by 2 because we have
         %a half distribution.
         u0 = u0 * 2;
      end
   else
      u0 = 0;
   end
else
   u0 = interp1(g_fokkerplanck.X, g_fokkerplanck.p0, x);
end
% --------------------------------------------------------------
end
function [pl, ql, pr, qr] = pdex1bc(~, ul, ~, ur, ~)
%global g_grind;
%reflecting boundaries pl=0 ql=1 pr=0 qr=1
%absorbing boundaries pl=ul ql=0 pr=ur qr=0
if g_grind.fokkerplanck.bl == 1 %kind of left boundary?
   %absorbing
   pl = ul;
   ql = 0;
else
   %reflecting boundaries
   pl = 0;
   ql = 1;
end
if g_grind.fokkerplanck.br == 1 %kind of left boundary?
   %absorbing
   pr = ur;
   qr = 0;
else
   %reflecting boundaries
   pr = 0;
   qr = 1;
end
end
end %main function
% --------------------------------------------------------------
function [g_c, g_f, g_s] = pdex1pde(g_x, t, g_u, g_DuDx)
global g_grind
if ischar(g_grind.fokkerplanck.sigma)
   for g_i = 1:length(g_grind.fokkerplanck.sigmavars)
      ix = i_getno(g_grind.fokkerplanck.sigmavars{g_i});
      if ix.isvar
         eval(sprintf('%s=g_x(%d);', g_grind.fokkerplanck.sigmavars{g_i}, ix.no));
      else
         eval(sprintf('global %s;', g_grind.fokkerplanck.sigmavars{g_i}));
      end
   end
   g_sigma = eval(g_grind.fokkerplanck.sigma);
   if ~isempty(g_grind.fokkerplanck.diffsigma)
      g_diffsigma =  eval(g_grind.fokkerplanck.diffsigma);
   else
      g_deltasigma = 1e-8;
      for g_i = 1:length(g_grind.fokkerplanck.sigmavars)
         ix = i_getno(g_grind.fokkerplanck.sigmavars{g_i});
         if ix.isvar
            eval(sprintf('%s=g_x(%d)+g_deltasigma;', g_grind.fokkerplanck.sigmavars{g_i}, ix.no));
         end
      end
      g_diffsigma = (eval(g_grind.fokkerplanck.sigma) - g_sigma) / g_deltasigma;
   end
else
   g_sigma = g_grind.fokkerplanck.sigma;
   g_diffsigma = 0;
end
g_Fx = feval(g_grind.hfun.curr, t, g_x);
%Fx=x*(1-x/K)-C*x^p/(x^p+H^p);
g_c = 1;
%we have to take the derivative of (1/2)*sigma(x)^2*P(x) to x using the
%product and chain rule:
%sigma(x)*P(x)*(diff(sigma(x), x))+(1/2)*sigma(x)^2*(diff(P(x), x))
g_f  = -g_Fx * g_u + g_sigma*g_u*g_diffsigma+ (g_sigma^2) / 2* g_DuDx;
g_s = 0;
end

function addreplaydata(H)
global g_fokkerplanck;
hax=findobj(H,'type','axes');
ud = get(hax, 'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @onstart;
ud.replay.onend = [];
ud.replay.settings=struct('tvar','t','tlim',[g_fokkerplanck.t(1) g_fokkerplanck.t(end)],'numt',length(g_fokkerplanck.t));
set(hax, 'userdata', ud);
end
function onstart(hax)
global  g_fokkerplanck;
if ishandle(hax)
   ud = get(hax, 'userdata');
   if isfield(g_fokkerplanck, 't')&&~isempty(g_fokkerplanck.t)
      ud.replay.settings=struct('tvar','t','tlim',[g_fokkerplanck.t(1) g_fokkerplanck.t(end)],'numt',length(g_fokkerplanck.t));
   end
   set(hax, 'userdata', ud);
   i_figure(get(hax, 'parent'));
end
end
function t = replaycallback(hax,avar, relt)
global g_fokkerplanck;
t = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, 't')
   ud = get(hax, 'userdata');
   t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
   ser = get(hax, 'children');
   if isfield(ud, 'replay')&&isfield(g_fokkerplanck,'t')
      it=find(g_fokkerplanck.t >= t, 1, 'first');
      set(ser(end), 'YData', g_fokkerplanck.Ys(it,:,1));
   end
end
end
function update_g_grind
global g_grind
if ischar(g_grind.fokkerplanck.bl)
   d = str2double(g_grind.fokkerplanck.bl);
   if isnan(d)
      if strncmpi(g_grind.fokkerplanck.bl, 'a', 1)
         g_grind.fokkerplanck.bl = 1;
      else
         g_grind.fokkerplanck.bl = 0;
      end
   else
      g_grind.fokkerplanck.bl = d;
   end
end
if ischar(g_grind.fokkerplanck.br)
   d = str2double(g_grind.fokkerplanck.br);
   if isnan(d)
      if strncmpi(g_grind.fokkerplanck.br, 'a', 1)
         g_grind.fokkerplanck.br = 1;
      else
         g_grind.fokkerplanck.br = 0;
      end
   else
      g_grind.fokkerplanck.br = d;
   end
end
s = str2double(g_grind.fokkerplanck.sigma);
if ~isnan(s)
   g_grind.fokkerplanck.sigma = s;
end
drawnow;
if ischar(g_grind.fokkerplanck.sigma)
   %    g_grind.fokkerplanck.sigmavars={};
   s1 = parsed_equation(g_grind.fokkerplanck.sigma);
   g_grind.fokkerplanck.sigmavars = s1.symvar;
   %     [s,vartype]=parseeq(sigma1);
   %     g_grind.fokkerplanck.sigmavars=unique(s(vartype==50));
else
   g_grind.fokkerplanck.sigmavars = {};
end
end
function makefigures
global g_fokkerplanck g_grind;
avar = i_statevars_names(1);
% Extract the first solution component as u.
u = g_fokkerplanck.Ys(:, :, 1);

hfig = i_makefig('fokkerplanck1');
set(hfig, 'Name','Fokker-Planck pdf')
%figure
h = area(g_fokkerplanck.X, u(end, :));
ud = get(gca, 'userdata');
ud.meta=struct('func','fokkerplanck','xname',avar,'yname','pdf','zname','');
set(gca, 'userdata', ud);
set(h, 'FaceColor', [0.8 0.8 1]);
i_plotdefaults(hfig);
xlabel(g_grind.statevars.names{1})
ylabel('probability density function')
addreplaydata(hfig);
replayall;
hfig = i_makefig('fokkerplanck2');
set(hfig, 'Name','Fokker-Planck time plot')
nt = size(g_fokkerplanck.Ys, 1);
g_fokkerplanck.t = g_fokkerplanck.t(1:nt);
Ps = zeros(nt, length(g_fokkerplanck.isaddles) - 1);
for i = 1:nt
   for j = 1:length(g_fokkerplanck.isaddles) - 1
      ran = g_fokkerplanck.isaddles(j):g_fokkerplanck.isaddles(j + 1);
      if length(ran) > 1
         Ps(i, j) = trapz(g_fokkerplanck.X(ran), g_fokkerplanck.Ys(i, ran));
      end
   end
end

plot(g_fokkerplanck.t, Ps);
ud = get(gca, 'userdata');
ud.meta=struct('func','fokkerplanck','xname','t','yname','probability','zname','');
set(gca, 'userdata', ud);
i_plotdefaults(hfig);
xlabel('Time (t)')
ylabel('Probability')
leg = cell(length(g_fokkerplanck.isaddles) - 1, 1);
for j = 1:length(g_fokkerplanck.isaddles) - 1
   leg{j} = sprintf('%s in [%g,%g]',avar,g_fokkerplanck.X(g_fokkerplanck.isaddles(j)),g_fokkerplanck.X(g_fokkerplanck.isaddles(j + 1)));
end
legend(leg);
ylim([0 1]);
xlim([g_fokkerplanck.t(1) g_fokkerplanck.t(end)]);
replayall('-plot','t');


isabs = ones(size(g_fokkerplanck.isaddles));
isabs(1) = g_grind.fokkerplanck.bl;
isabs(end) = g_grind.fokkerplanck.br;
bothabs = isabs(1:end - 1) & isabs(2:end);
if any(bothabs)
   hfig = i_makefig('fokkerplanck3');
   set(hfig, 'Name','Fokker planck -Exit prob.')
   hold on;
   cleft = 'r';
   cright = 'b';
   for i = 1:length(isabs) - 1
      if length(g_fokkerplanck.uxl{i}) > 1
         plot(g_fokkerplanck.xs{i}, g_fokkerplanck.uxl{i}, ['-' cleft])
      else
         plot(g_fokkerplanck.xs{i}, zeros(size(g_fokkerplanck.xs{i})), ['-' cleft])
      end
      if length(g_fokkerplanck.uxr{i}) > 1
         plot(g_fokkerplanck.xs{i}, g_fokkerplanck.uxr{i}, ['-' cright])
      else
         plot(g_fokkerplanck.xs{i}, zeros(size(g_fokkerplanck.xs{i})), ['-' cright])
      end
      hfig1 = cleft;
      cleft = cright;
      cright = hfig1;
   end
   ud = get(gca, 'userdata');
   ud.meta=struct('func','fokkerplanck','xname',avar,'yname','P exit','zname','');
   set(gca, 'userdata', ud);
   i_plotdefaults(hfig);
   plotsaddles(isabs);
   xlabel(avar);
   ylabel('Exit probability');
   ylim([ - 0.005 1.005]);
end
if any(isabs)
   hfig = i_makefig('fokkerplanck4');
   set(hfig, 'Name','Fokker planck - Mean first passage time')
   hold on;
   cleft = 'r';
   cright = 'b';
   for i = 1:length(isabs) - 1
      h1 = 0;
      if isabs(i)&&isabs(i+1)
         h1 = 1;
      end
      if length(g_fokkerplanck.Txl{i}) > 1
         plot(g_fokkerplanck.xs{i}(1 + h1:end-h1), g_fokkerplanck.Txl{i}(1 + h1:end-h1) ./ g_fokkerplanck.uxl{i}(1 + h1:end-h1), ['-' cleft])
      end
      if length(g_fokkerplanck.Txr{i}) > 1
         plot(g_fokkerplanck.xs{i}(1 + h1:end-h1), g_fokkerplanck.Txr{i}(1 + h1:end-h1) ./ g_fokkerplanck.uxr{i}(1 + h1:end-h1), ['-' cright])
      end
      if length(g_fokkerplanck.Txb{i}) > 1
         plot(g_fokkerplanck.xs{i}, g_fokkerplanck.Txb{i}, ':k')
      end
       ud = get(gca, 'userdata');
       ud.meta=struct('func','fokkerplanck','xname',avar,'yname','P exit','zname','');
       set(gca, 'userdata', ud);
       hfig1 = cleft;
      cleft = cright;
      cright = hfig1;
   end
   i_plotdefaults(hfig);
   plotsaddles(isabs);
   xlabel(avar);
   ylabel('Mean first passage time');
end
%Psend=zeros(1,length(g_fokkerplanck.isaddles) - 1);
%Txmean=Psend;
for j = 1:length(g_fokkerplanck.isaddles) - 1
   ran = g_fokkerplanck.isaddles(j):g_fokkerplanck.isaddles(j + 1);
   if length(ran) > 1
      fprintf('\nBasin %d from %s=%g to %g\n',j,avar, g_fokkerplanck.X(g_fokkerplanck.isaddles(j)),g_fokkerplanck.X(g_fokkerplanck.isaddles(j + 1)))
      Psend = trapz(g_fokkerplanck.X(ran), g_fokkerplanck.Ys(end, ran));
      fprintf('Probability to be in basin %d: %g\n', j, Psend);
      if length(g_fokkerplanck.Txb{j}) > 1
         Tx = interp1(g_fokkerplanck.xs{j}, g_fokkerplanck.Txb{j}, g_fokkerplanck.X(ran));
      elseif length(g_fokkerplanck.Txl{j}) > 1
         Tx = interp1(g_fokkerplanck.xs{i}, g_fokkerplanck.Txl{j}, g_fokkerplanck.X(ran));
      elseif length(g_fokkerplanck.Txr{j}) > 1
         Tx = interp1(g_fokkerplanck.xs{j}, g_fokkerplanck.Txr{j}, g_fokkerplanck.X(ran));
      else
         Tx = nan;
      end
      Txmean =  trapz(g_fokkerplanck.X(ran),g_fokkerplanck.Ys(end, ran) .* Tx) / Psend;
      fprintf('Mean exit time from basin %d: %g\n', j, Txmean);
   end
end
end
function plotsaddles(isabs)
global g_fokkerplanck;
ys = get(gca, 'ylim');
cleft = 'r';
cright = 'b';
dx = (g_fokkerplanck.X(end) - g_fokkerplanck.X(1)) / 100;
xlim([g_fokkerplanck.X(1) - dx g_fokkerplanck.X(end) + dx]);
for i = 1:length(isabs)
   xx = g_fokkerplanck.X(g_fokkerplanck.isaddles(i));
   if isabs(i)
      plot([xx, xx], ys, [':', cleft]);
      %             if i==1||i==length(isabs)
      %                 text(xx,ys(2),'absorbing','color',cleft);
      %             else
      %                 text(xx,ys(2),'saddle','color',cleft);
      %             end
   else
      plot([xx, xx], ys, ['-', cleft]);
      %             text(xx,ys(2),'reflecting','color',cleft);
      
   end
   h = cleft;
   cleft = cright;
   cright = h;
end
end
