%NDX2SPARSE - create sparse matrix with ones on the indexed values
%
% Usage:
%   A=NDX2SPARSE(NDX)
%
% See also:
% SPARSE2NDX, CREATE_NETWORK
function A=ndx2sparse(ndxs)
A=sparse([ndxs(:,1);ndxs(:,2)],[ndxs(:,2);ndxs(:,1)],1);
A(A>1)=1;