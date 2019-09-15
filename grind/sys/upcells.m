%UPCELLS   Shift up the elements of a matrix/vector
%  Efficient function to shift the elements of a matrix one row up. This
%  function is used to access neighbors of a matrix or vector state variable. 
%  Optionally, the first row is neighboring to the last row.
%
%  Usage:
%  Usage:
%  UPCELLS(N,BORDER) - shift the matrix N 1 position to the right. If BORDER is 'neumann' the first row 
%  borders to itself (Neumann boundaries, you can also set a value to the derivative), if BORDER is 
%  'periodic' the first row borders the last.
%  UPCELLS(N,'dirichlet',VALUE) - if BORDER is 'dirichlet' the boundary condition VALUE is used at the border.
%
%
%  See also model, rightcells, leftcells, downcells, neighborcells
%
%   Reference page in Help browser:
%      <a href="matlab:commands('upcells')">commands upcells</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function result=upcells(N,bordered,aValue)
if nargin<2
    bordered='neumann';
end
if isnumeric(bordered) %backward compatibility
    switch bordered
        case 0
            bordered='periodic';
        case 1
            bordered='neumann';
        case 2
            bordered='dirichlet';
    end
end
if nargin<3
    aValue=0;
end
s1=size(N,1);

switch bordered
    case 'neumann'
        result=[N(1,:)-aValue;N(1:s1-1,:)];
    case 'periodic'
        result=[N(s1,:);N(1:s1-1,:)];
    case 'dirichlet'
        result=[zeros(1,size(N,2))+aValue;N(1:s1-1,:)];
    otherwise
        error('grind:upcells','Unknown boundary condition')
end

