%LYAPMAX Simple Ellner algorithm for Lyapunov exponent
%   Use simtime to get a time series with a small fixed time step
%
%   Usage:
%   lyapsimple - gives the lyapunov

function [LE,Jac]=lyapmax(varargin)
global g_grind t g_t g_Y;
fieldnams={'step', 'n>0', 'step for the lyapunov exponent ',1}';
args=i_parseargs(fieldnams,'step','',varargin);
if nargin == 0
   if g_grind.solver.isdiffer || isnan(g_grind.tstep)
      args.step = 1;
   else
      args.step =  g_grind.ndays / g_grind.tstep;
   end
end
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(t, g_grind.ndays, N0, 1);
end
tt = (g_t(1):args.step:g_t(end))';
if length(tt)~=length(g_t)
    YY = interp1(g_t, g_Y, tt);
else
    YY=g_Y;
end
Jac = zeros(length(tt), length(N0)^2);
U = ones(size(N0));
epsilon=1E-6;
U=U/norm(U)*epsilon;
LE = 0;
donum = isempty(g_grind.syms.Jacobian);
for i = 1:length(tt)
   N0 = YY(i, :)';
   J = i_calcjac(donum, 1,N0);
   Jac(i, :) = J(:)';
   if g_grind.solver.isdiffer
      U = J * U;
   else
       %the time step should be small enough for this
      U = expm(args.step * J) * U;
   end
   maxU = norm(U);%max(abs(U));
   LE = LE + log(maxU/epsilon);
   U = U / maxU*epsilon;
end
LE = LE / (tt(end) - tt(1));
disp('LE = ');
disp(LE);
