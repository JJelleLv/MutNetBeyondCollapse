%PARLOOKUP   find a value in the first column of a table
%  This function finds a key value in the first column of a table. It returns the remainder.
% 
% 
%  Usage:
%  PARLOOKUP(Table,Key,extrapolate), lookup a key in the table (first column) with or without
%  extrapolation outside the data range
%
%  Example:
%  A=[1 3;5 6;10 5;15 3];
%  parlookup(A,4) returns 5.25 (interpolated)
%  parlookup(A,20) returns 0 (outside range)
%  parlookup(A,20,1) returns 3 (extrapolated)
%
%  See also insim
%
%   Reference page in Help browser:
%      <a href="matlab:commands('parlookup')">commands parlookup</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function apar = parlookup(tabl, key1, extrapolate, NaNisZero)
if nargin < 4
   NaNisZero = 1;
end
[N, C] = size(tabl);
if ~isempty(tabl)
if key1 < tabl(1, 1) || key1 > tabl(N, 1)
   if (nargin < 3) || ~extrapolate
      apar = NaN .* zeros(1, size(tabl, 2) - 1);
   elseif key1 < tabl(1, 1)
      apar = tabl(1, 2:C);
   else
      apar = tabl(N, 2:C);
   end
else
     i = 1;
     while (i<N) && tabl(i, 1) <= key1
        i = i + 1;
     end
%      i=find(tabl(:, 1) > key1,1);
%     if isempty(i)
%         i=N;
%     end
   apar = tabl(i - 1, 2:C)+(tabl(i, 2:C)-tabl(i - 1, 2:C))./(tabl(i, 1) - tabl(i - 1, 1))*(key1-tabl(i - 1, 1));
%    a = (tabl(i, 2:C) - tabl(i - 1, 2:C)) ./ (tabl(i, 1) - tabl(i - 1, 1));
%    b = tabl(i - 1, 2:C) - a .* tabl(i - 1, 1);
%    par = a .* key1 + b;
%    par-par1
end
else
    apar=NaN;
end
if NaNisZero
   apar(isnan(apar)) = 0;
end

