function s = i_othervars(N0, iX, iY, iZ)
global g_grind t;
if nargin == 2
   iY = -1;
end

if nargin < 4
   iZ = -1;
end

if g_grind.solver.nonautonomous
    s = sprintf('t=%g ',t);
else
    s='';
end

n= g_grind.statevars.dim;
if n>20
   n=20;
   shortened=1;
else
   shortened=0;
end

for i = 1:n
   if (i ~= iX) && (i ~= iY) && (i ~= iZ)
      if isempty(s)
         s = [i_statevars_names(i) ' = ' num2str(N0(i))];
      else
         s = sprintf('%s and %s = %0.4g',s, i_statevars_names(i),N0(i));
      end

   end

end

if shortened
   s=[s ' and more...'];
end


