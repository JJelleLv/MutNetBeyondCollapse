%RIGHTCELLS   Shift right the elements of  a matrix
%  Efficient function to shift the elements of a matrix one column to the right. This
%  function is used to access neighbors of a matrix state variable. 
%  Optionally, the first column is neighboring the last column.
%
%  Usage:
%  RIGHTCELLS(N,BORDER) - shift the matrix N 1 position to the right. If BORDER is 'neumann' the last column 
%  borders to itself (Neumann boundaries, you can also set a value to the derivative), if BORDER is 
%  'periodic' the last column borders the first (periodic).
%  RIGHTCELLS(N,'dirichlet',VALUE) - if BORDER is 'dirichlet' the boundary condition VALUE is used at the border.
%
%  See also model, leftcells, upcells, downcells, neighborcells
%
%   Reference page in Help browser:
%      <a href="matlab:commands('rightcells')">commands rightcells</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function result=rightcells(N,bordered,aValue)

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
s2=size(N,2);

switch bordered
    case 'neumann'
        result=[N(:,2:s2),N(:,s2)+aValue];
    case 'periodic'
        result=[N(:,2:s2),N(:,1)];
    case 'dirichlet'
        result=[N(:,2:s2),zeros(size(N,1),1)+aValue];
    otherwise
        error('grind:rightcells','Unknown boundary condition')
end

