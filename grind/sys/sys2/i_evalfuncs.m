%before using this function run ones as i_evalfuncs to update g_grind
function [g_l_F] = i_evalfuncs(g_l_t)
global g_Y g_grind g_t t;
t0=t;
g_l_F = {};
if nargin == 0
   %check the dimensions of the functions
   if isfield(g_grind, 'funcnames')&&isfield(g_grind.funcnames,'dims')
      %Note that this shortcut implies that the dimensions of functions may not change, see setdimension
      return;
   end

   if ~isempty(g_grind.pars)
      eval(i_globalstr(g_grind.pars));
   end

   if g_grind.statevars.vector
      eval(i_globalstr(g_grind.statevars.vectnames));
   else
      eval(i_globalstr(g_grind.statevars.names))
   end

%   global t;
%   g_l_t = 1;
else
   if ~isempty(g_t) && (g_l_t <= length(g_t)) && (g_l_t > 0) &&~isempty(g_grind.funcnames.names)
      t = g_t(g_l_t); 
   else
      return;
   end

end

if ~isempty(g_grind.pars)
   eval(i_globalstr(g_grind.pars));
end

if ~isempty(g_Y) && (nargin ~= 0)
   if g_grind.statevars.vector
      for g_l_j = 1:length(g_grind.statevars.dims)
         eval(sprintf('%s = reshape(g_Y(g_l_t,g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to),g_grind.statevars.dims{%d}.dim1,g_grind.statevars.dims{%d}.dim2);', ...
            g_grind.statevars.vectnames{g_l_j}, g_l_j, g_l_j, g_l_j, g_l_j));
      end

   else
      for g_l_j = 1:g_grind.statevars.dim
         eval(sprintf('%s = g_Y(g_l_t,%d);',i_statevars_names(g_l_j),g_l_j));
      end

   end

end

NO = i_initvar;
eval(g_grind.funcs);
i_keep(NO);
for g_l_j = 1:length(g_grind.funcnames.names)
   eval(sprintf('g_l_F{g_l_j}=%s;', g_grind.funcnames.names{g_l_j}));
end

if nargin == 0
   g_l_i = 1;
   for g_l_j = 1:length(g_l_F)
      g_grind.funcnames.dims{g_l_j}.dim1 = size(g_l_F{g_l_j}, 1);
      g_grind.funcnames.dims{g_l_j}.dim2 = size(g_l_F{g_l_j}, 2);
      g_grind.funcnames.dims{g_l_j}.from = g_l_i;
      g_l_i = g_l_i + numel(g_l_F{g_l_j});
      g_grind.funcnames.dims{g_l_j}.to = g_l_i - 1;
   end

end

t=t0;

