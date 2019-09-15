% BYYEARCLASS = order results by yearclass
%   Special function to order the results of age structured models, such that each year class
%   can be followed.
%
%   Example:
%   out -1 byyearclass(lengths);
%
%   See also:
%   OUT, AX
%
function res = byyearclass(A, lag)
global g_Y;
if nargin < 2
   lag = 1;
end
s = size(A);
maxA = max(A);
s1 = size(g_Y);
res = ones(s) * NaN;
k = 0;
i1 = 1;
%find the shifts
% if most of the differences between the current and the next are larger than shifted, then there has been
% a shift of the yearclasses (boxcartrain)
Ynext = g_Y(2:end, :);
shftY = abs(g_Y(1:end-1,:) - Ynext(:, [1 + lag:s1(2) + 1 - lag 1:lag]));
diffY = abs(g_Y(1:end-1, :) - Ynext)-1E-4;
shifts = find(sum(shftY > diffY, 2) < s1(2) * 0.9);
for j = 1:length(shifts)
   i = shifts(j);
   res(i1:i, :) = A(i1:i, 1 + mod((1:s(2)) - 1 + k, s(2)));
   k = k + lag;
   i1 = i + 1;
end
if i1 < s(1)
   res(i1:s(1)-1, :) = A(i1:s(1)-1, 1 + mod((1:s(2)) - 1 + k, s(2)));
end
for j = 1:length(shifts)
   i = shifts(j);
   % one value is NaN to disconnect different classes
   % (a bit ticky)
   cond =  (abs(res(i + 1, :) - res(i, :)) > maxA * 0.5);
   res(i, cond) = NaN;
end

