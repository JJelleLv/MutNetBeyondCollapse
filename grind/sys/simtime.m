%SIMTIME   Set the simulation duration
%   Enter the following information:
%   - Start time - sets the start time (can also be assigned as <a href="matlab:commands t">t</a> =..)
%   - Number of units to simulate - number of time units for default runs
%   - Number of values for output (leave empty for maximal) - you can reduce
%     the number of steps in the output, for example as Poincare sections for
%     yearly cycles. Leave this value empty for maximal output.
%
%    Usage:
%    SIMTIME - enter the data interactively.
%    SIMTIME T NDAYS  - sets the starting time to T and the number of time units to NDAYS
%    SIMTIME T NDAYS TSTEP - sets the number of steps for output to TSTEP also.
%    SIMTIME('t',10,'ndays',100,'tstep',1000) - arguments as parameter/value pairs.
%    SIMTIME('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'ndays' [number>0] - number of days for a run
%     't' [number] - start time
%     'tstep' [number>2|isnan(number)] - number of steps for output (NaN = maximum)
%
%
%   See also ru, time
%
%   Reference page in Help browser:
%      <a href="matlab:commands('simtime')">commands simtime</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function simtime(varargin)
global t g_grind;
fieldnams={'t', 'n', 'start time',0;...
   'ndays', 'n>0', 'number of days for a run',1000;...
   'tstep', 'n>2|isnan(n)', 'number of steps for output (NaN = maximum)',NaN}';
args=i_parseargs(fieldnams,'t,ndays,tstep',{},varargin);
if isfield(args,'t')
   t = args.t;
end
if isfield(args,'ndays')&&~isempty(args.ndays)
   g_grind.ndays = args.ndays;
end
if nargin == 3
   g_grind.tstep = args.tstep;
end
if nargin == 0
   i_parcheck;
   prompt = {'Start time:', ...
      'Number of time units to simulate', ...
      'Number of values for output (leave empty for maximal)'};
   answer = cell(3);
   answer{1} = num2str(t);
   answer{2} = num2str(g_grind.ndays);
   if ~isnan(g_grind.tstep)
      answer{3} = num2str(g_grind.tstep);
   else
      answer{3} = '';
   end
   answer = inputdlg(prompt, 'Simulator time', 1, answer);
   if ~isempty(answer)
      t = i_checkstr(answer{1});
      g_grind.ndays = i_checkstr(answer{2});
      if ~isempty(strtrim(answer{3}))
         g_grind.tstep = i_checkstr(answer{3});
      else
         g_grind.tstep = NaN;
      end
   end
end
