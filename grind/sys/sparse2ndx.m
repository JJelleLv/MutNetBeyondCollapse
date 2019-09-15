% SPARSE2NDX - Create index with non-zeros values of sparse matrix
%
% Usage:
%   NDX=SPARSE2NDX(A)
%
% See also:
% NDX2SPARSE, CREATE_NETWORK
function ndxs=sparse2ndx(A)
A=triu(A); %use the upper triangle of A (as A should be symmetric)
[ndx1, ndx2]=ind2sub(size(A), find(A>0)); % get the subscript values (ndx1, ndx2) where A == 1
ndxs = [ndx1, ndx2]; %save the subscripts
 