%KE   Keep the final state after a run
%   Fix the final state of the last simulation as new starting
%   point.
%  
%   Usage:
%   KE - saves the final state
%   KE -start (-s) - saves the start state
%   KE T - saves the state closest to time T (no interpolation is used)
%   KE('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     't' [number] - saves the state closest to time t
%   KE('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-s' - saves the start state
%      
%  
%   See also val, ru
%
%   Reference page in Help browser:
%      <a href="matlab:commands('ke')">commands ke</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function ke(varargin)
global g_Y g_t g_grind;
if (nargin>0)
   fieldnams={'t', 'n', 'saves the state closest to time t',[]}';
   args=i_parseargs(fieldnams,'t','-s',varargin);
   if strcmp(args.opts,'-s')
      ndx=1;
   elseif isfield(args,'t')&&~isempty(args.t)
      aa=abs(g_t-args.t);
      ndx=find(aa==min(aa));
   end
else
   ndx=size(g_Y,1);
end
if ~isempty(g_Y)
   i_keep(g_Y(ndx, :));
   if ~isempty(g_grind.permanent)
      defpermanent('-s',defpermanent('-p',g_t(ndx))); 
   end
else
   warning('GRIND:ke:norun','No simulation run to keep, use <a href="matlab:time">time</a> or <a href="matlab:ru">ru</a> to simulate first');
end
