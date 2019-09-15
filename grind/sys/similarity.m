function dist = similarity(V1, V2, method, par1)
if nargin < 3
   method = 'EUCLEDEAN';
end
method = upper(method);
switch method
 case {'EUCLEDEAN','EUCL','E'}
   dist = sqrt(sum((V1 - V2).^2));
 case {'CITYBLOCK','CITY','C'}
   dist = sum(abs(V1 - V2));
 case {'CHEBYCHEV','CHEBY','CH'}
   dist = max(abs(V1 - V2));
 case {'ABSSINE','SINE','S'}
   if nargin == 4
      dist = sqrt(1 - (sum(V1 .* V2) ./ sqrt(sum(V1.^2)*par1)).^2); 
   else
      dist = sqrt(1 - (sum(V1 .* V2) ./ sqrt(sum(V1.^2)*sum(V2.^2))).^2);
   end
   %Dist := Sqrt(1 - Sqr(SumXY / Sqrt(SumX2 * SumY2)));
 otherwise
   error('GRIND:similarity:UnknownSim','Unknown similarity measure');
end


