%NEIGHBORCELLS   Get all 4 or 8 neighbors (optionally random order)
%Get all neighbors as a three dimensional array (neighbors are always in the 3rd dimension). 
%If you get the neighbors of a vector only two values are given.
%
%  Usage:
%  NEIGHBORCELLS(N,4,BORDER,VALUE) - Get 4 neighbors. If BORDER is 1 the first/last column/row 
%  borders to itself (Neumann boundaries, you can also set a value to the derivative), if BORDER is 
%  0 the first column/row borders the last (periodic) BORDER is 2 for Dirichlet boundary conditions 
%  (set to a fixed value) (see also: <a href="matlab:help leftcells">leftcells</a>).
%  NEIGHBORCELLS(N,6) - Get 6 neighbors of hexagonial grid: all cells [r,c] connect to the 4 neighbors
%   even column cells connect to [r+1,c-1] and [r+1,c+1] and  odd column cells connect to [r-1,c-1] [r-1,c+1]
%  NEIGHBORCELLS(N,8) - Get 8 neighbors, BORDER =0 by default.
%  NEIGHBORCELLS(N,-4) - Get 4 neighbors in random order.
%  NEIGHBORCELLS(N,-8) - Get 8 neighbors in random order.
%  NEIGHBORCELLS(N,-2,BORDER) - Get 2 neighbors (of a vector) in random order.
%
%  See also model, rightcells, leftcells, upcells, downcells
%
%   Reference page in Help browser:
%      <a href="matlab:commands('neighborcells')">commands neighborcells</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function Result = neighborcells(N, nneighbors, bordered, value)
siz = size(N);
if nargin < 2
   if any(siz) == 1
      nneighbors = 2;
   else
      nneighbors = 4;
   end
end
if nargin < 3
   if any(siz)==1&&nneighbors>=0&&nneighbors<=2&&(nargin>1)
      bordered = nneighbors;
      nneighbors = 2;
   else
      bordered = 1;
   end
end
if nargin < 4
   value = 0;
end
switch nneighbors
 case 2
   if siz(1) ~= 1
      Result = cat(3, upcells(N, bordered, value), downcells(N, bordered, value));
   else
      Result = cat(3, leftcells(N, bordered, value), rightcells(N, bordered, value));
   end
 case 4
   Result = zeros(siz(1), siz(2), 4);
   Result(:,:,1) = leftcells(N, bordered, value);
   Result(:,:,2) = rightcells(N, bordered, value);
   Result(:,:,3) = upcells(N, bordered, value);
   Result(:,:,4) = downcells(N, bordered, value);
 case 6 %hexagonal grid
   %all cells [r,c] connect to [r,c-1] [r-1,c] [r,c+1], [r+1,c]
   %even column cells (r,c) connect to [r+1,c-1] [r+1,c+1]
   %odd column cells (r,c) connect to [r-1,c-1] [r-1,c+1]
   Result = zeros(siz(1), siz(2), 6);
   L = leftcells(N, bordered, value);
   R = rightcells(N, bordered, value);
   Result(:, :, 1) = L;
   Result(:, :, 2) = R;
   Result(:,:,3) = upcells(N, bordered, value);
   Result(:,:,4) = downcells(N, bordered, value);
   LU = upcells(L, bordered, value);
   LD = downcells(L, bordered, value);
   LU(:, 2:2:siz(2)) = LD(:, 2:2:siz(2)); %odd cells go left - down
   Result(:, :, 5) = LU;
   RU = upcells(R, bordered, value);
   RD = downcells(R, bordered, value);
   RU(:, 2:2:siz(2)) = RD(:, 2:2:siz(2)); %odd cells go right - down
   Result(:, :, 6) = RU;
   %   Result = cat(3, leftcells(N, bordered, value), rightcells(N, bordered, value), upcells(N, bordered, value), downcells(N, bordered, value));
 case 8
   L = leftcells(N, bordered, value);
   R = rightcells(N, bordered, value);
   Result = cat(3, L, R, upcells(N, bordered, value), downcells(N, bordered, value),upcells(L, bordered, value), ...
      downcells(L, bordered, value),    upcells(R, bordered, value),  downcells(R, bordered, value));
 case  {  - 2,  - 4,  - 6,  - 8}
   Result = neighborcells(N,   - nneighbors, bordered, value) ;
   Result = getrand3(Result);
 otherwise
   error('neigborcell:grind:nneighbors', 'Number of neighborcells must be 2, 4, 6 or 8 (use random order for something in between)');
end
function res = getrand3(A)
siz = size(A);
[~, ndx] = sort(rand(siz), 3);
res = A(sub2ind(siz, repmat(transpose(1:siz(1)), [1, siz(2), siz(3)]), repmat((1:siz(2)), [siz(1), 1, siz(3)]), ndx));

