%RETURNTIME   Estimates the time necessary to reach a stable node
%   From the current initial conditions, the time is estimated till the change in the 
%   model is negligible. This command analyses the results of the last run. (if there 
%   is no last run or if parameters have changed, it calls RU).
%  
%   Usage:
%   RETURNTIME - estimates the return time based on the current SIMTIME settings and a default
%   value for the tolerance (1E-8).
%   RET=RETURNTIME - store the return time in the variable RET.
%   RETURNTIME ABSTOL - if the distance to the equlibrium < ABSTOL, the model is returned to equilibrium.
%   RETURNTIME ABSTOL NDAYS - use ABSTOL and simulate NDAYS time units to find an equilibrium.
%   RETURNTIME('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'abstol' [number>0] - if the distance to the equlibrium < ABSTOL, the model is returned to equilibrium.
%     'method' [change | equil | both] - method to determine if equilibrium is reached 'both', 'change' or 'equil'
%     'ndays' [number>0] - number of days for run.
%  
%   
%   See also returntime2d, ru, simtime, findeq
%
%   Reference page in Help browser:
%      <a href="matlab:commands('returntime')">commands returntime</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [rettime,YY] = returntime(varargin)
%(err1, maxt, methodopt), methods: both,change,equil
global g_Y t g_t g_grind;
fieldnams={'abstol', 'n>0', 'if the distance to the equlibrium &lt; ABSTOL, the model is returned to equilibrium.',1E-8;...
   'ndays', 'n>0', 'number of days for run.',g_grind.ndays;...
   'method', 'e[change|equil|both]', 'method to determine if equilibrium is reached ''both'', ''change'' or ''equil''','both'}';
args=i_parseargs(fieldnams,'abstol,ndays,method','',varargin);
i_parcheck;
if ~isfield(args,'abstol')
  args.abstol = 1E-8;
end
opts=solver('-properties');
if ~opts.hasfixedstep&&opts.AbsTol>args.abstol
    solver('abstol',args.abstol/10,'reltol',args.abstol/10);
end
if ~isfield(args,'method')
    args.method='both';
end
if isfield(args,'ndays')
    args.ndays= g_grind.ndays;
end
switch lower(args.method)
    case 'change'
       method=0;
    case 'equil'
       method=1;
    otherwise
       method=2;
end 
if method==2
   err2=args.abstol*2;
elseif method==1
   err2=args.abstol;
   args.abstol=1E40; %very large                       number: will never been reached
elseif method==0
   err2=1E40;
end
N0 = i_initvar;
if ~isnan(args.ndays)||i_settingschanged(N0)
   i_ru(t, args.ndays, N0, 1);
end
if length(g_t)-g_grind.tstep<2
   tt=g_t;
   YY=g_Y;
else
   tt=g_t(1):1:g_t(size(g_t,1));
   YY=interp1(g_t,g_Y,tt);
   if size(g_Y,2)==1
      YY=transpose(YY);
   end   
end
iret = 1;
conv=0;
for i = 1:size(YY, 1) - 1
   differ = mean(abs(YY(i, :) - YY(i + 1, :)));
   diff2 = mean(abs(YY(i, :) - YY(size(YY,1), :)));
   conv= (differ<args.abstol) & (diff2<err2);
   if conv
 %     disp(sprintf('iret %0.5g diff %0.5g diff2 %0.5g',[iret,diff,diff2]));
      break;
   end
   iret = i;
end
if ~conv
   rettime=tt(length(tt))-tt(1);
   if (nargout == 0)
      disp('no equilibrium reached (cannot detect cycles)');
   end
else
   rettime =tt(iret) - tt(1);
end

