function [dim1, dim2] = i_getvardims(avar)
global g_grind;
if ~g_grind.statevars.vector
   dim1=1;
   dim2=1;
else
   for k = 1:size(g_grind.statevars.dims, 2)
      if strcmp(avar,g_grind.statevars.vectnames{k})
         dim1=g_grind.statevars.dims{k}.dim1;
         dim2=g_grind.statevars.dims{k}.dim2;
         return;
      end

   end

   dim1=[];
   dim2=[];
end

