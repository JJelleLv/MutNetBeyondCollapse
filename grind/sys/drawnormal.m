%A=drawnormal(dim1,dim2,mu,sd)
%
%draw a matrix (A) from normal distributions
%
function A=drawnormal(dim1,dim2,mu,sd)
if nargin==3
   sd=mu;
   mu=dim2;
   dim2=dim1;
end
A=randn(dim1,dim2)*sd+mu;
