%SETDIAGON - set the diagonal of a matrix to a value
%
% Usage:
% A = setdiagon(A,1) - sets the diagonal of matrix A to 1
%
% Example:
% A = setdiagon(rand(3),0)
%
% A =
%
%         0    0.6093    0.1789
%    0.8542         0    0.3465
%    0.1177    0.4292         0
function A=setdiagon(A,aval)
if nargin<2
    aval=1;
end
%efficient way using indexing
%works also for non-square matrices and aval can be a vector
siz=size(A);
A(1:siz(1)+1:min(siz(1)*(siz(1)+1),numel(A)))=aval;
%A=A+diag((aval-diag(A))); %a bit faster but not flexible


