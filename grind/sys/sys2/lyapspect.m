%LYAPSPECT   Calculate Lyapunov exponent spectrum
%   Calculate the spectrum of Lyapunov exponents (lambda) and the Lyapunov dimension 
%   (=Kaplan-Yorke dimension). These parameters expresses the sensitivity to initial conditions. 
%   If there is one (clearly) positive Lyapunov exponent the model is chaotic. More positive 
%   Lyapunov exponents means "hyperchaos". The algorithm is based on:
%   A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano, "Determining Lyapunov Exponents 
%   from a Time Series," Physica D, Vol. 16, pp. 285-317, 1985.
%   The Jacobian is approximated numerically if the user has not entered the equations 
%   (see <a href="matlab:help enterjac">enterjac</a>). The use of symbolic Jacobians is recommended (enterjac -sym).
%
%
%   Usage:
%   LYAPSPECT - calculate lyapunov spectrum with g_grind.ndays days. Default step = 0.1, StepJac=0.3
%   LYAPSPECT NDAYS STEP - calculate NDAYS days. The step for QR algorithm=0.1.
%   RES=LYAPSPECT(N,STEP) = returns a structure with the last exponents (RES.LAMBDA) and 
%   the last Lyapunov dimension (RES.DIMENSION). Furthermore it includes the data of convergence: the times (RES.TIMES), lyapunov exponents (RES.LAMBDAS) for each step.
%       
%   LYAPSPECT('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'ndays' [number>0] - number of days to analyse (set default with <a href="matlab:help simtime">simtime</a>).
%     'nout' [number>=0] - step for output to the command window (default=100)
%     'step' [number>0] - step for the QR algorithm (default=0.5)
%     'symbolic' [logical] - if true, a symbolic Jacobian is generated if possible and necessary (default=true)
%
%   See also lyapunov, enterjac, lorenzmap, takens, poincaremap, poincaresect
%
%   Reference page in Help browser:
%      <a href="matlab:commands('lyapspect')">commands lyapspect</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [res] = lyapspect(varargin)
%(ndays,step,stepjac)
global g_grind t;
%tic;
if isfield(g_grind,'ndays')
    ndays=g_grind.ndays;
else
    ndays=1000;
end
fieldnams={'ndays', 'n>0', 'number of days to analyse (set default with <a href="simtime.htm">simtime</a>).',ndays;...
   'step', 'n>0', 'step for the QR algorithm (default=0.5)',0.5;...
   'nout', 'n>=0', 'step for output to the command window (default=100)',100;...
   'symbolic', 'l', 'if true, a symbolic Jacobian is generated if possible and necessary (default=true)',true}';
args=i_parseargs(fieldnams,'ndays,step','',varargin);
if ~isfield(args,'ndays')
    args.ndays = g_grind.ndays;
end
if ~isfield(args,'symbolic')
    args.symbolic= ~g_grind.statevars.vector;
end
if ~isfield(args,'step')
    args.step = 0.5;
end
if ~isfield(args,'nout')
    args.nout=100;
end
if args.symbolic&&isempty(g_grind.syms.Jacobian)
    enterjac('-sym')
end
N0 = i_initvar;
if g_grind.solver.nonautonomous
    if strcmp(questdlg('Might not find correct lyapunov exponents in nonautonomous equation, continue?','Error','Yes','No','No'),'No')
        error('GRIND:lypaspect:NonAutonomous','Make equation autonomous by adding tau''=1 and replacing all t by tau');
    end
end
if ~args.symbolic
    rhs_ext_fcn=i_getodehandle('ExtendedSyst','numonly');
else
    rhs_ext_fcn=i_getodehandle('ExtendedSyst');
end
fcn_integrator=str2func(g_grind.solver.name);
[tim_1,lambda_1]=i_wolf_lyapunov(numel(N0),rhs_ext_fcn,fcn_integrator,t,args.step,args.ndays,N0,args.nout);
hfig = i_makefig('lyapspect');
plot(tim_1, lambda_1);
i_plotdefaults(hfig);
xlabel('time');
ylabel('\lambda''s');
Lambda = sort(lambda_1(end,:), 'descend');
%To calculate the Lyapunov dimension (or Kaplan-Yorke dimension)
LESum = Lambda(1);
LD = 0;
d=numel(N0);
if (d > 1 && Lambda(1) > 0)
    for N = 1:d - 1
        if Lambda(N + 1) ~= 0
            LD = N + LESum / abs(Lambda(N + 1));
            LESum = LESum + Lambda(N + 1);
            if LESum < 0
                break;
            end
        end
    end
end
if nargout>0
    res.time=tim_1;
    res.lamdas=lambda_1;
    res.lambda=lambda_1(end,:);
    res.dimension=LD;
else
    for i = 1:d
        fprintf('lambda(%g) = %g\n', i, lambda_1(end,i));
    end
    fprintf('Lyapunov dimension (or Kaplan-Yorke dimension) = %g\n', LD);
end

function [times,lambdas]=wolf_lyapunov(d,rhs_ext_fcn,fcn_integrator,t,step,ndays,N0,nout)

L = d + d^2;
A = eye(d);
IC = [N0(:); A(:)];
T = t - step;
%g_grind.lyapspect.L = d^2;
% g_grind.lyapspect.dtJac = args.stepjac * args.step;
% g_grind.lyapspect.d = d;
% g_grind.lyapspect.num = isempty(g_grind.syms.Jacobian);
% g_grind.lyapspect.J = i_calcjac(g_grind.lyapspect.num, 1, N0);
% g_grind.lyapspect.tJac = t + g_grind.lyapspect.dtJac;
% g_grind.lyapspect.odehandle=i_getodehandle('normal');
sumR = zeros(d, 1);
lyaps = zeros(ceil(ndays / step), d);
times = zeros(1, ceil(ndays / step) - 1);
lambdas = zeros(ceil(ndays / step) - 1, d);
for i = 1:ceil(ndays / step)
   T = T + step;
   [~,X]=feval(fcn_integrator, rhs_ext_fcn,[T,T+0.5*step,T+step],IC');
   IC = X(size(X, 1), :);
   P = reshape(IC(d + 1:L), d, d);
   %QR decomposition
   [A, R] = qr(P);
   IC(d + 1:L) = A(:);
   R = abs(diag(R));
   for j = 1:size(A, 1)
      if R(j) < 1E-40
         R(j) = 0;
      else
         R(j) = log(R(j));
      end
   end
   lyaps(i, :) = R';
   sumR = sumR + R;
   if T > t
      lambdas(i - 1, :) =  sumR' / T;
      times(i - 1) = T;
   end
end
% hfig = figure;
% plot(times, lambdas);
% i_plotdefaults(hfig);
% xlabel('time');
% ylabel('\lambda''s');
% lyaps = lyaps ./ T;
% lambda = sum(lyaps);
% Lambda = sort(lambda, 'descend');
% %To calculate the Lyapunov dimension (or Kaplan-Yorke dimension)
% LESum = Lambda(1);
% LD = 0;
% if (d > 1 && Lambda(1) > 0)
%    for N = 1:d - 1
%       if Lambda(N + 1) ~= 0
%          LD = N + LESum / abs(Lambda(N + 1));
%          LESum = LESum + Lambda(N + 1);
%          if LESum < 0
%             break;
%          end
%       end
%    end
% end;
% for i = 1:d
%    fprintf('lambda(%g) = %g\n', i, lambda(i));
% end;
% fprintf('Lyapunov dimension (or Kaplan-Yorke dimension) = %g\n', LD);
% if nargout > 0
%    lambda1 = lambda;
% end;
% g_grind.lyapspect = [];
% toc;
