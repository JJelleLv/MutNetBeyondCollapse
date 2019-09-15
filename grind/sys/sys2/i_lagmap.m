function Res = i_lagmap(M)
s = size(M);
Res = zeros(s(1), s(2));
if (s(1)>0)&&(s(2)>0)
   Res(s(1), :) = NaN+zeros(1, s(2));
   Res(1:s(1) - 1, :) = M(2:s(1), :);
end

