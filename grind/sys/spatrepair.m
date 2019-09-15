%SPATREPAIR = time to recover from a disturbance in a spatial model
%
%[res,repairing]=repairtime(err,dist,ref)
%     err= tolerated difference between the two points
%     dist= coordinates of disturbed point
%     ref = coordinates of reference point
%
% result: res = time to recovery (NaN = no recovery)
%         repairing: 1 = repairing, 0 expanding/expanded
%         ntrial = maximum number of appended runs
%
function [Reptime,Repairing] = spatrepair(err, dist, ref, ntrial)
global g_grind g_Y g_t t;
if ~g_grind.statevars.vector
   error('GRIND:sparrepair:NoVector','designed for vector/matrix variables only')
end
if nargin < 1
   err = 0.01;
else
   err =  i_checkstr(err);
end
if (nargin < 2)||isempty(dist)
   dist = round([g_grind.statevars.dims{1}.dim1, g_grind.statevars.dims{1}.dim2] ./ 2);
else
   dist =  i_checkstr(dist);
end
if nargin  < 4
   ntrial = 1;
else
   ntrial =  i_checkstr(ntrial);
end
if (nargin < 2)||isempty(ref)
   ref = [1, 1];
else
   ref =  i_checkstr(ref);
end
oldmode = g_grind.solver.addmode;
isok = 0;
trial = 1;
while ~isok
   N0 = i_initvar;
   if i_settingschanged(N0)
      i_ru(t, g_grind.ndays, N0, 1);
   end
   refs = zeros(size(g_Y,1),length(g_grind.statevars.vectnames));
   dists = zeros(size(g_Y,1),length(g_grind.statevars.vectnames));
   for i = 1:length(g_grind.statevars.vectnames)
      inddist = g_grind.statevars.dims{i}.from + sub2ind([g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2], dist(1), dist(2))-1;
      indref = g_grind.statevars.dims{i}.from + sub2ind([g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2], ref(1), ref(2))-1;
      refs(:, i) = g_Y(:, indref);
      dists(:, i) = g_Y(:, inddist);
   end
   differ = ((refs - dists).^2).^0.5;
 %  diffs = max(diff,[],2);
   f = find(differ < err);
   refchange = max(((refs(1, :) - refs(size(refs, 1), :)).^2).^0.5);
   distchange = max(((dists(1, :) - dists(size(dists, 1), :)).^2).^0.5);
   totchange = max(((g_Y(round(size(dists, 1) * 0.95), :) - g_Y(size(dists, 1), :)).^2).^0.5);
   repairing = 1;
   if ~isempty(f)
      reptime = g_t(f(1)) - t;
      if refchange > err
         repairing = 0;
         fprintf('Global regime shift in %d time units\n', reptime);
      else
         fprintf('Repaired in %d time units\n', reptime);
      end
   else
      if totchange < err
         reptime = Inf;
         repairing = 0;
         disp('Stable but not repaired')
      else
         reptime = NaN;
         if (distchange < refchange) || (distchange < totchange)
            repairing = 0; 
            fprintf('Still expanding after %d time units\n', g_t(end))
         else
            fprintf('Still repairing/expanding after %d time units\n', g_t(end))
         end
      end
   end
   isok = ~isnan(reptime) | (trial >= ntrial);
   if ~isok
      g_t = [g_t(1, 1); g_t(end - 1:end, 1)];
      g_Y = [g_Y(1, :); g_Y(end - 1:end, :)];
      addmode('on','-silent');
      trial=trial+1;
   end
end
g_grind.solver.addmode = oldmode;
if nargout > 0
   Reptime = reptime;
   Repairing = repairing;
end
