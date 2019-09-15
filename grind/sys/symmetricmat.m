function M = symmetricmat(A, down2up)
M = A;
siz = size(M, 1);
if siz ~= size(M, 2)
   error('GRIND:symmetricmat:square','Matrix must be square');
end
if (nargin > 1) && down2up
   for i = 1:siz
      for j = i + 1:siz
         M(i, j) = M(j, i);
      end
   end
else
   for i = 1:siz
      for j = i + 1:siz
         M(j, i) = M(i, j);
      end
   end
end
