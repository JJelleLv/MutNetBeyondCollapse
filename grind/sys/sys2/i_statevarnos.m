function [sfrom, sto] = i_statevarnos(avar)
global g_grind;
sfrom=[];
fcolon = strfind(avar, ':');
fbrack = strfind(avar, '(');
if g_grind.statevars.vector && (isempty(fbrack) || ~isempty(fcolon))
   if ~isempty(fbrack)
      avar1 = avar(1:fbrack(1)-1);
   else
      avar1 = avar;
   end
   sto = 0;
   for i = 1:length(g_grind.statevars.dims)
      sfrom = sto + 1;
      sto = sto + g_grind.statevars.dims{i}.dim1 * g_grind.statevars.dims{i}.dim2;
      if strcmp(g_grind.statevars.vectnames{i}, avar1)
         if ~isempty(fcolon)
          %  l = sub2ind2d([g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2], avar);
             l = i_findindices(avar,[g_grind.statevars.dims{i}.dim1, g_grind.statevars.dims{i}.dim2]);
             sfrom = l - 1 + sfrom;
            sto=[];
         end

         break
      end

   end

else
   sto = 0;
   for i = 1:g_grind.statevars.dim
      sto = sto + 1;
      if strcmp(i_statevars_names(i), avar)
         sfrom = sto;
         break
      end

   end

end

