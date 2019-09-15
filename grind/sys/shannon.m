function H = shannon(x, dim, threshold)
%SHANNON   Shannon-Weaver diversity index.
%   For matrices, SHANNON(X) is a column vector containing the shannon index of
%   each row. For vectors, SHANNON(X) is the shannon index of the whole vector.
%   (also called the Shannon Wiener index)
%
%   SHANNON(X,DIM, threshold) takes the index along the dimension DIM of X.
%
%   H = -sum(p(i) * ln(p(i)))
%
%   in which p(i)=X(i)/sum(X)
%
if nargin < 3
   threshold = 0.01;
end
if nargin < 2
   % Determine which dimension SHANNON will use
   if min(size(x)) == 1
      dim = find(size(x) == 1); dim=dim(1);
   else
      dim  = 1;
   end
   if isempty(dim)
      dim = 1;
   end
end
if dim == 2
   sums = repmat(sum(x, 1), size(x, 1), 1);
   dim2 = 1;
else
   sums = repmat(sum(x, 2), 1, size(x, 2));
   dim2 = 2;
end
p = x ./ sums;
lnp = zeros(size(p));
lnp(x > threshold) = log(p(x > threshold));
plnp = p .* lnp;
H = -sum(plnp, dim2);
