%DOWNCELLS   shift down the elements of a matrix/vector
%  Efficient function to shift the elements of a matrix one row down. This
%  function is used to access neighbors of a matrix or vector state variable. 
%  Optionally, the first row is neighboring to the last row.
%
%  Usage:
%  DOWNCELLS(N,BORDER) - shift the matrix N 1 position to the right. If BORDER is 'neumann' the last row 
%  borders to itself (Neumann boundaries, you can also set a value to the derivative), if BORDER is 
%  'periodic' the last row borders the first.
%  DOWNCELLS(N,'dirichlet',VALUE) - if BORDER is 'dirichlet' the boundary condition VALUE is used at the border.
% 
%
%  See also model, rightcells, leftcells, upcells, neighborcells
%
%   Reference page in Help browser:
%      <a href="matlab:commands('downcells')">commands downcells</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function result=downcells(N,bordered,aValue)
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
        result=[N(2:s1,:);N(s1,:)+aValue];
    case 'periodic'
        result=[N(2:s1,:);N(1,:)];
    case 'dirichlet'
        result=[N(2:s1,:);zeros(1,size(N,2))+aValue];
    otherwise
        error('grind:leftcells','Unknown boundary condition')
end

