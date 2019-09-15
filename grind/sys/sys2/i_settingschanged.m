function [changed,changedpars] = i_settingschanged(N0, ndays)
global g_grind g_Y g_t;
changedfields={};
if nargin == 0
   settings = i_getsettings;
elseif nargin==1
   settings = i_getsettings(N0);
elseif nargin == 2
   settings = i_getsettings(N0,ndays);
end

if isempty(g_grind.checks.lastsettings)
   changed = 1;
elseif g_grind.solver.addmode
   % in addmode always run, update the run if the initial conditions change
   nperm = 0;
   if isfield(g_grind, 'permanent')
      for i = 1:size(g_grind.permanent)
         p = evalin('base', char(g_grind.permanent{i}.name));
         nperm = nperm + numel(p);
      end

   end

   if (nargin > 0) && isdifferent(N0,g_grind.checks.lastsettings.initvar) || (size(g_Y, 1) < 3) || (size(g_Y, 1)~=length(g_t))
      g_Y = [];
      g_t = [];
   end

   changed = 1;
elseif isempty(g_t) || ((nargin == 2) && (abs(g_t(length(g_t)) - settings.solver(1) + settings.solver(2)) > 1))
   changed = 1;
   [~,changedfields] = struccmp(settings,g_grind.checks.lastsettings);
else
   [thesame,changedfields] = struccmp(settings,g_grind.checks.lastsettings);
   if ~thesame
       g_grind.checks.no=g_grind.checks.no+uint16(ismember(fieldnames(settings),changedfields));
   end

   changed=~thesame|| (size(g_Y, 1) < 3);
end

if nargout > 1
   if ~changed
      changedpars = 0;
   else
      changedpars= any(strcmp(changedfields,'pars'));
   end

end


function res = isdifferent(A, B)
res=~min(size(A) == size(B));
if isempty(B)&&isempty(A)
   res = 0;
elseif ~res
   compar = A == B;
   res = ~min(compar);
   if res
      res = ~min(compar + isnan(A) .* isnan(B));
   end

end

