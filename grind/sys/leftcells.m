%LEFTCELLS   shift left the elements of  a matrix
%  Efficient function to shift the elements of a matrix one column to the left. This
%  function is used to access neighbors of a matrix state variable. 
%  Optionally, the first column is neighboring the last column.
%
%  Usage:
%  LEFTCELLS(N,BORDER) - shift the matrix N 1 position to the right. If BORDER is 'neumann' the first column 
%  borders to itself (Neumann boundaries, you can also set a value to the derivative), if BORDER is 
%  'periodic' the first column borders the last.
%  LEFTCELLS(N,'dirichlet',VALUE) - if BORDER is 'dirichlet' the boundary condition VALUE is used at the border.
%
%
%  See also model, rightcells, upcells, downcells, neighborcells
%
%   Reference page in Help browser:
%      <a href="matlab:commands('leftcells')">commands leftcells</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function result=leftcells(N,bordered,aValue)
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
        result=[N(:,1)-aValue,N(:,1:s2-1)];
    case 'periodic'
        result=[N(:,s2),N(:,1:s2-1)];
    case 'dirichlet'
        result=[zeros(size(N,1),1)+aValue,N(:,1:s2-1)];
    otherwise
        error('grind:leftcells','Unknown boundary condition: %s',bordered)
end

