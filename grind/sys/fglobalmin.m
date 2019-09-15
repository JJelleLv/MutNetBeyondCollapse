%FGLOBALMIN Scalar bounded global function minimization.
%   Using the very robust Shuffled Complex Evolution (SCE-UA) global optimization algorithm.
%   written by Q. Duan, 9/2004 interface adapted by E. van Nes
%
%References:
%  Duan, Q., S. Sorooshian, and V.K. Gupta, 1992. Effective and efficient global optimization 
%    for conceptual rainfall-runoff models. Water Resources Research 28: 1015-1031.
%  Duan, Q., S. Sorooshian, and V.K. Gupta, 1994. Optimal use of the SCE-UA global 
%    optimization method for calibrating watershed models. Journal of Hydrology 158:265-284.
%
%Usage:
%  [bestx, bestf, exitflag, output] = fglobalmin(funfcn, ax, bx, options, P1,P2..)
%  bestx = fglobalmin(funfcn, ax, bx)
%  bestx = fglobalmin(funfcn, ax, bx, [], P1,P2..)
%
%Arguments:
%  funfcn = the function to be evaluated (can be inline function or m file)
%  ax = the lower bound of the parameters (can be a vector for each parameter)
%  bx = the upper bound of the parameters (can be a vector for each parameter)
%  options = options, see optimset or below
%  P1,P2.. = additional parameters for funfcn
%Output:
%  bestx    = the optimized parameter array at the end;
%  bestf    = the objective function value corresponding to the optimized parameters
%  exitflag = convergence? 1 = yes, 0 = no
%  output   = structure with details about convergence
%
%Options:
%  options.Display ('off'|'final'|'iter') = display nothing, final result or each evolution ('final')
%  options.MaxFunEvals = maximum number of function evaluations (numberOfVariables*1000)
%  options.MaxIter = maximum number of evolution loops before convergence(5)
%  options.MaxPCGIter (or options.NCmplx)= number of complexes (sub-populations) (4)
%  options.TolFun = the percentage change allowed in MaxIter loops before convergence (0.1)
%  options.TolX = convergence if normalized parameter space less than this value (0.001)
%  options.TypicalX = optional vector with first trial ([])
%
% 
%  optimset('fglobalmin') returns the default options
%
%See also: <a href="matlab:help fminbnd">fminbnd</a>, <a href="matlab:help fminsearch">fminsearch</a>, <a href="matlab:help optimset">optimset</a>, <a href="matlab:help optimget">optimget</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('fglobalmin')">commands fglobalmin</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [bestx, bestf, exitflag, output] = fglobalmin(funfcn, ax, bx, options, varargin)
%FGLOBALMIN Scalar bounded global function minimization.
% Using the very robust Shuffled Complex Evolution (SCE-UA) global optimization algorithm.
% written by Q. Duan, 9/2004 interface adapted by E.van Nes
%
% References:
%   Duan, Q., S. Sorooshian, and V.K. Gupta, 1992. Effective and efficient global optimization 
%     for conceptual rainfall-runoff models. Water Resources Research 28: 1015-1031.
%   Duan, Q., S. Sorooshian, and V.K. Gupta, 1994. Optimal use of the SCE-UA global 
%     optimization method for calibrating watershed models. Journal of Hydrology 158:265-284.
%
% Usage:
%   [bestx, bestf, exitflag, output] = fglobalmin(funfcn, ax, bx, options, P1,P2..)
%   bestx = fglobalmin(funfcn, ax, bx)
%   bestx = fglobalmin(funfcn, ax, bx, [], P1,P2..)
%
% Arguments:
%   funfcn = the function to be evaluated (can be inline function or m file)
%   ax = the lower bound of the parameters (can be a vector for each parameter)
%   bx = the upper bound of the parameters (can be a vector for each parameter)
%   options = options, see optimset or below
%   P1,P2.. = additional parameters for funfcn
% Output:
%   bestx    = the optimized parameter array at the end;
%   bestf    = the objective function value corresponding to the optimized parameters
%   exitflag = convergence? 1 = yes, 0 = no
%   output   = structure with details about convergence
%
% Options:
%   options.MaxFunEvals = maximum number of function evaluations (numberOfVariables*1000)
%   options.TolX = convergence if normalized parameter space less than this value (0.001)
%   options.MaxIter = maximum number of evolution loops before convergency (5)
%   options.TolFun = the percentage change allowed in MaxIter loops before convergency (0.1)
%   options.MaxPCGIter (or options.NCmplx)= number of complexes (sub-populations) (4)
%   options.TypicalX = optional vector with first trial ([])
%   options.Display ('off'|'final'|'iter') = display nothing, final result or each evolution ('final')

%  
%   optimset('fglobalmin') returns the default options
%
% See also: FMINBND,FMINSEARCH,OPTIMSET,OPTIMGET
defaultopt = optimset('display','final','maxiter',5,...
   'MaxFunEvals',1000,'TolX',0.001,'TolFun',0.1,'MaxPCGIter',4);
if nargin < 4
   options = [];
end

if nargin == 1 && nargout <= 1 && isequal(funfcn, 'defaults')
   bestx = defaultopt;
   return
end

if nargin < 3
   error('GRIND:fglobalmin:ArgError','Please specify upper and lower boundaries (ax and bx)');
end
if length(ax) ~= length(bx)
   error('GRIND:fglobalmin:rangerror','Number of upper bounds is not equal to the number of lower bounds');
end

%numberOfVariables = length(ax);

if isfield(options,'NCmplx') %more logical name for ngs, but is not in optimset
   options = optimset(options, 'MaxPCGIter', options.NCmplx);
end

options = optimset(defaultopt, options);

printtype = optimget(options, 'Display');
if ~ischar(printtype)
   displ = printtype;
else
   switch printtype
    case {'none','off','0'}
      displ = 0;
    case {'iter','2'}
      displ = 2;
    case {'final','1'}
      displ = 1;
    case 'simplex'
      displ = 3;
    otherwise
      displ = 1;
   end
end
kstop = optimget(options, 'maxiter');
maxn = optimget(options, 'maxfunevals');
peps  = optimget(options, 'TolX'); %peps
pcento =  optimget(options, 'TolFun'); %Pcento
ngs = optimget(options, 'MaxPCGIter'); %Ngs
x0 = optimget(options, 'typicalx'); % vector with first trial
outputfcn = optimget(options,'OutputFcn');
hasoutputfcn= ~isempty(outputfcn);
    
if ischar(pcento)
   pcento = eval(pcento);
end
if ischar(ngs)
   ngs = eval(ngs);
end

if ischar(maxn)
   maxn = eval(maxn);
end
if ischar(kstop)
   kstop = eval(kstop);
end

if length(x0) ~= length(ax)
   x0 = [];
end

% Convert to inline function as needed.
if ischar(funfcn)
    funfcn=str2func(funfcn);
end
    %funfcn = fcnchk(funfcn, length(varargin));

% Initialize SCE parameters:
exitflag = 1; %assume we converge
nopt = length(ax);
npg = 2 * nopt + 1;
nps = nopt + 1;
nspl = npg;
%mings = ngs;
npt = npg * ngs;

bound = bx - ax;

% Create an initial population to fill array x(npt,nopt):
x = zeros(npt, nopt);
for i = 1:npt
   x(i, :) = ax + rand(1, nopt) .* bound;
end

if ~isempty(x0)
   x(1, :) = x0;
end


icall = 0;
xf=zeros(1,npt);
for i = 1:npt
   xf(i) = feval(funfcn, x(i, :), varargin{:});
   %  xf(i) = functn(nopt,x(i,:));
   icall = icall + 1;
end


% Sort the population in order of increasing function values;
[xf, idx] = sort(xf);
x = x(idx, :);

% Record the best and worst points;
bestx = x(1, :); bestf = xf(1);
worstx = x(npt, :); worstf = xf(npt);
%BESTF=bestf; BESTX=bestx;ICALL=icall;

% Compute the standard deviation for each parameter
%xnstd = std(x);

% Computes the normalized geometric range of the parameters
gnrng = exp(mean(log((max(x) - min(x)) ./ bound)));
if displ > 1
   disp('The Initial Loop: 0');
   disp(['Best F  : ' num2str(bestf)]);
   disp(['Best X  : ' num2str(bestx)]);
   disp(['Worst F : ' num2str(worstf)]);
   disp(['Worst X : ' num2str(worstx)]);
   disp(' ');
   
   % Check for convergency;
   if icall >= maxn
      exitflag = 0;
      if displ > 0
         disp('*** Optimization search terminated because the limit');
         disp('on the maximum number of trials ');
         disp(maxn);
         disp('has been exceeded.  Search was stopped at trial number:');
         disp(icall);
         disp('of the initial loop!');
      end
   end
   
   if gnrng < peps
      if disp > 0
         disp('The population has converged to a prespecified small parameter space');
      end
   end
end
% Begin evolution loops:
nloop = 0;
criter = [];
criter_change = 1e5;

while icall < maxn && gnrng > peps && criter_change > pcento
   nloop = nloop + 1;
   stop=false;
   % Loop on complexes (sub-populations);
   for igs = 1: ngs
      
      % Partition the population into complexes (sub-populations);
      k1 = 1:npg;
      k2 = (k1 - 1) * ngs + igs;
      cx(k1, :) = x(k2, :);
      cf(k1) = xf(k2);
      
      % Evolve sub-population igs for nspl steps:
      for loop = 1:nspl
         
         % Select simplex by sampling the complex according to a linear
         % probability distribution
         lcs=zeros(1,nps);
         lcs(1) = 1;
         for k3 = 2:nps
            for iter = 1:1000
               lpos = 1 + floor(npg+0.5-sqrt((npg+0.5)^2 - npg * (npg+1) * rand));
               %idx=find(lcs(1:k3 - 1) == lpos);  if isempty(idx) 
               if any(lcs(1:k3 - 1) == lpos)
                   break;
               end
            end
            lcs(k3) = lpos;
         end
         lcs = sort(lcs);
         
         % Construct the simplex:
        % s = zeros(nps, nopt);
         s=cx(lcs, :); sf = cf(lcs);
         
         [snew, fnew, icall] = cceua(funfcn, s, sf, ax, bx, icall, varargin{:});
         if hasoutputfcn
            stop=outputfcn(x, bestx, bestf);
            if stop
               exitflag=0;
               icall=maxn+1;
               break;
            end
        end
         
         % Replace the worst point in Simplex with the new point:
         s(nps, :) = snew; sf(nps) = fnew;
         
         % Replace the simplex into the complex;
         cx(lcs, :) = s;
         cf(lcs) = sf;
         
         % Sort the complex;
         [cf, idx] = sort(cf); cx=cx(idx, :);
         
         % End of Inner Loop for Competitive Evolution of Simplexes
      end

      % Replace the complex back into the population;
      x(k2, :) = cx(k1, :);
      xf(k2) = cf(k1);
      if stop
          break;
      end      
      % End of Loop on Complex Evolution;
   end
   
   % Shuffled the complexes;
   [xf, idx] = sort(xf); x=x(idx, :);
   
   % Record the best and worst points;
   bestx = x(1, :); bestf = xf(1);
   worstx = x(npt, :); worstf = xf(npt);
   
   % Compute the standard deviation for each parameter
   %xnstd = std(x);
%    if hasoutputfcn
%        stop=outputfcn(x, bestx, bestf);
%        if stop
%            exitflag=0;
%            icall=maxn+1;
%        end
%    end
   % Computes the normalized geometric range of the parameters
   gnrng = exp(mean(log((max(x) - min(x)) ./ bound)));
   if displ > 1
      disp(['Evolution Loop: ' num2str(nloop) '  - #calls - ' num2str(icall)]);
      disp(['Best F  : ' num2str(bestf)]);
      disp(['Best X  : ' num2str(bestx)]);
      disp(['Worst F : ' num2str(worstf)]);
      disp(['Worst X : ' num2str(worstx)]);
      disp(' ');
      
      % Check for convergency;
      if icall >= maxn
         if displ > 0
            disp('*** Optimization search terminated because the limit');
            disp(['on the maximum number of trials ' num2str(maxn) ' has been exceeded!']);
            
         end
         exitflag = 0;
      end
      
      if gnrng < peps
         if displ > 0
            disp('The population has converged to a prespecified small parameter space');
         end
      end
   end
   criter = [criter; bestf];
   if (nloop >= kstop)
      criter_change = abs(criter(nloop) - criter(nloop - kstop + 1)) * 100;
      criter_change = criter_change / mean(abs(criter(nloop - kstop + 1:nloop)));
      if (displ > 0) && (criter_change < pcento)
         disp(['The best point has improved in last ' num2str(kstop) ' loops by ', ...
            'less than the threshold ' num2str(pcento) '%']);
         disp('Convergency has achieved based on objective function criteria!!!')
      end
   end
   
   % End of the Outer Loops
end
if displ > 0
   disp(['Search was stopped at trial number: ' num2str(icall)]);
   disp(['Normalized geometric range = ' num2str(gnrng)]);
   disp(['The best point has improved in last ' num2str(kstop) ' loops by ', ...
      num2str(criter_change) '%']);
end
output.iterations = nloop;
output.funcCount = icall;
output.bestf = bestf;
output.bestx = bestx;
output.worstf = worstf;
output.worstx = worstx;
output.criterion_change_perc = criter_change;
output.size_parameterspace = gnrng;
output.algorithm = 'SCE-UA (Q. Duan)';


% END of Subroutine sceua
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [snew, fnew, icall] = cceua(funfcn, s, sf, bl, bu, icall,  varargin)
%  This is the subroutine for generating a new point in a simplex
%
%   s(.,.) = the sorted simplex in order of increasing function values
%   s(.) = function values in increasing order
%

[nps, nopt] = size(s);
n = nps;
%m = nopt;
alpha = 1.0;
beta = 0.5;

% Assign the best and worst points:
%sb = s(1, :);
%fb = sf(1);
sw = s(n, :); 
fw = sf(n);

% Compute the centroid of the simplex excluding the worst point:
ce = mean(s(1:n - 1, :));

% Attempt a reflection point
snew = ce + alpha * (ce-sw);

% Check if is outside the bounds:
ibound = 0;
s1 = snew - bl; 
%idx = find(s1 < 0);if ~isempty(idx)
if any(s1(:)<0)
   ibound = 1; 
end
s1 = bu - snew; 
%idx = find(s1 < 0); if ~isempty(idx)
if any(s1(:)<0)
   ibound = 2; 
end

if ibound  >= 1
   snew = bl + rand(1, nopt) .* (bu - bl);
end
fnew = feval(funfcn, snew, varargin{:});
%fnew = functn(nopt,snew);

icall = icall + 1;

% Reflection failed; now attempt a contraction point:
if fnew > fw
   snew = sw + beta * (ce-sw);
   % fnew = functn(nopt,snew);
   fnew = feval(funfcn, snew, varargin{:});
   icall = icall + 1;
   
   % Both reflection and contraction have failed, attempt a random point;
   if fnew > fw
      snew = bl + rand(1, nopt) .* (bu - bl);
      %fnew = functn(nopt,snew);
      fnew = feval(funfcn, snew, varargin{:});
      icall = icall + 1;
   end
end

% END OF CCE
return;


