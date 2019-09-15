%sub2ind2d = 2D version of sub2ind
% index=sub2ind2d(size(N),r1,r2,c1,c2) gives the indices of
% N(r1:r2,c1:c2)
%
%Alternative usage:
%   index=sub2ind2d(size(N),'N(10:20,10:20)')
function index = sub2ind2d(sizeN, r1, r2, c1, c2)
if nargin == 2
   avar = r1;
   f = strfind(avar, '(');
   f1 = strfind(avar, ':');
   f2 = strfind(avar, ',');
   f3 = strfind(avar, ')');
   if isempty(f)
      f = 0;
   end
   if isempty(f3)
      f3 = length(avar) + 1;
   end
   if isempty(f2)
      r1 = str2double(avar(f(1) + 1:f1(1) - 1));
      r2 = str2double(avar(f1(1) + 1:f3(1) - 1));
      c1 = -1;
      c2 = -1;
   else
      r1 = str2double(avar(f(1) + 1:f1(1) - 1));
      r2 = str2double(avar(f1(1) + 1:f2(1) - 1));
      c1 = str2double(avar(f2(1) + 1:f1(2) - 1));
      c2 = str2double(avar(f1(2) + 1:f3(1) - 1));
   end
   if isempty(r1)
      r1 = 1;
   end
   if isempty(r2)
      r2 = sizeN(1);
   end
   if isempty(c1)
      c1 = 1;
   end
   if isempty(c2)
      c2 = sizeN(2);
   end
end
rows = r1:r2;
cols = c1:c2;
b1 = repmat(rows, length(cols), 1);
b2 = repmat(transpose(cols), 1, length(rows));
if c1 == -1
   index = transpose(sub2ind(sizeN, b1(:)));
else
   index = transpose(sub2ind(sizeN, b1(:), b2(:)));
end

