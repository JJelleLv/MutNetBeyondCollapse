%MAT2CONNECT - matrix with connections between neighbors
%  This efficient function returns a sparse matrix defining the connections between a matrix 
%  the matrix gives the connections between the elements of the space matrix
%  if A= the space matrix and M the connections then
%  d*(M*A(:)) defines the diffusion between 4 neighbors.
%
% Usage:
%   M=MAT2CONNECT([SIZX,SIZY],FUNHANDLE,FUNARG1,FUNARG2 ...) SIZX is the 
%   size in X direction, SIZY is the size in the Y direction, FUNHANDLE is
%   handle to a function that finds the connections of each element (e.g. leftcells, rightcells, 
%   neighborcells), FUNARG1, FUNARG2 .. list of arguments of the FUNHANDLE (the first argument is assumed
%   to be the matrix
%   M=MAT2CONNECT(A,FUNHANDLE,FUNARG1,FUNARG2 ...) if numel(A)>2 A is
%   assumed to be the matrix.
%
%
%  See also neighborcells, leftcells, rightcells, upcells, downcells, connectmat    
%
function sp=mat2connect(A,matfunct,varargin)
if numel(A)==2
    siz=A;
else
    siz=size(A);
end

psiz=prod(siz);
%make a matrix with indices
A1=reshape(transpose(1:psiz),siz);
%apply the matfunc to the matrix with indices
%to get the connecting indices
if iscell(matfunct)
   %if matfunct is a cell, make 3 dim matrix
   A2=zeros(siz(1),siz(2),numel(matfunct));
   for i=1:numel(matfunct)
       f=matfunct{i};
       A2(:,:,i)=f(A1,varargin{:});
   end
else
   %just a simple function
   A2=matfunct(A1,varargin{:});
end

%make sparse result with true'th for connections (excluding self)
sp=sparse(psiz,psiz);
for i=1:size(A2,3)
    AA=A2(:,:,i);
    sp=max(sp,sparse(A1(:),AA(:),1,psiz,psiz));
end

