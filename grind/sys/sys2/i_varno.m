%i_varno returns the number of a statevar and [] if no statevar
function [varno, vecno] = i_varno(avar)
global g_grind;
varno = [];
vecno  = [];
start = 0;
if g_grind.statevars.vector
   %check for efficiency, especially necessary for grids
   found = 0;
   f1 = strfind(avar, '(');
   f2 = strfind(avar, ')');
 %  f3= strfind(avar,':');
   if (length(f1)==1)&&(length(f2)==1)&&(f2==length(avar))%&&isempty(f3)
      vecvar = avar(1:f1(1) - 1);
   else
      vecvar = avar;
   end
   for i = 1:size(g_grind.statevars.vectnames, 2)
      if strcmp(vecvar, g_grind.statevars.vectnames{i})
         found = 1;
         vecno = i;
         d = g_grind.statevars.dims{i}.dim1;
         break;
      end

      start = start + g_grind.statevars.dims{i}.dim1 * g_grind.statevars.dims{i}.dim2;
   end

   if found
      f1 = strfind(avar, '(');
      if isempty(f1)
         varno = start + 1;
      else
         f = strfind(avar, ',');
         f2 = strfind(avar, ')');
         f3 = strfind(avar, ':');
         if isempty(f2)
            error('GRIND:NoStatevar','"%s" is not a valid state variable',avar);
         end

         if isempty(f3)
            if ~isempty(f)
               f3 = [f(1), f2(1)];
            else
               f3 = f2(1);
            end

         end

         if isempty(f)
            varno = start + str2double(avar(f1(1) + 1:f3(1) - 1));
         else
            varno = start + str2double(avar(f1(1) + 1:f3(1) - 1)) + (str2double(avar(f(1) + 1:f3(2) - 1)) - 1) * d;
         end

      end

   end

else
   varno=find(strcmp(avar, g_grind.statevars.names));
   if isempty(varno)
       f1 = strfind(avar, '(1)');
       if ~isempty(f1)
           avar=avar(1:f1(1)-1);
       end

       varno=find(strcmp(avar, g_grind.statevars.names));
   end

%    for k = 1:length(g_grind.statevars.names)
%       if strcmp(avar, char(g_grind.statevars.names{k}))
%          varno = k;
%          return;
%       end
%    end
end

