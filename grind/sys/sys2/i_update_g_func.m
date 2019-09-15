function i_update_g_func
global g_func g_grind g_t;
if isempty(g_func)
   if ~isempty(g_grind) && isfield(g_grind, 'funcnames')&&~isempty(g_grind.funcnames.names)
      i_evalfuncs; %update g_grind.funcnames
      g_func = zeros(length(g_t), g_grind.funcnames.dims{end}.to);
      for g_l_t = 1:length(g_t)
         g_l_F = i_evalfuncs(g_l_t);
         for g_l_i =  1:length(g_grind.funcnames.names)
            g_l_f = g_grind.funcnames.dims{g_l_i};
            g_func(g_l_t, g_l_f.from:g_l_f.to) = transpose(g_l_F{g_l_i}(:));
         end

      end

   end

end

