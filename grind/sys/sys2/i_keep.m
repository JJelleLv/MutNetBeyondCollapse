% KE (KEEP)
% Fix the final state of the last simulation as new starting point
function i_keep(N0,NP)
global g_grind;
if ~isempty(N0)
   if g_grind.statevars.vector
      if size(N0, 2) > size(N0, 1)
         N0 = transpose(N0);
      end

      d = 1;
      for j = 1:length(g_grind.statevars.vectnames)
            avar = g_grind.statevars.vectnames{j};
            dim1 = g_grind.statevars.dims{j}.dim1;
            dim2 = g_grind.statevars.dims{j}.dim2;
            dim=dim1*dim2;
            if dim2>1
               assignin('base', avar, reshape(N0(d:d + dim-1),dim1,dim2))
            else
               assignin('base', avar, N0(d:d + dim-1));
            end

            d = d + dim;
      end

    else
      for i = 1:g_grind.statevars.dim
         multassignin('base', g_grind.statevars.names{i}, N0(i));
      end

   end
   if (nargin==2)&&~isempty(g_grind.permanent)
      defpermanent('-s',NP);
   end

   if g_grind.solver.haslags&&isfield(g_grind.solver,'newhist')&&~isempty(g_grind.solver.newhist)
       if isstruct(g_grind.solver.history)&&isstruct(g_grind.solver.newhist)
           g_grind.solver.history.x=[g_grind.solver.history.x g_grind.solver.newhist.x(2:end)];
           g_grind.solver.history.y=[g_grind.solver.history.y g_grind.solver.newhist.y(:,2:end)];
           g_grind.solver.history.yp=[g_grind.solver.history.yp g_grind.solver.newhist.yp(:,2:end)];
       else
          g_grind.solver.history=g_grind.solver.newhist;
       end

        g_grind.solver.newhist=[];
   end

end

    
