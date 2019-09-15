%Dyy   Double partial derivative to spatial y (=column) coordinate
%   Defines the double gradient in spatial coordinates [d2N/dy2] to be used for example for diffusion. By default
%   it uses the simple diffusion discretisation:
%   Dyy = Nj-1 - 2*Nj + Nj+1
%   Different boundary conditions can be selected using the command <a href="matlab:help definespace">definespace</a>.
%
%  
%   Usage:
%   Dyy(N) defines in a model a double partial derivative of the variable N with respect to 
%     the spatial y (=row) coordinates (except for single row variables, then y = column)..
%
%   See also model, Dxx, Dx, Dy, neighborcells, definespace
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('Dyy')">commands Dyy</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function differ=Dyy(N)
global g_grind;
if ~isfield(g_grind,'space')||isempty(g_grind.space)
    warning('grind:Dyy','No space defined, set to defaults');
    definespace('-d');
end

siz=size(N);
if siz(2)==1&&siz(1)~=1 %if vector x is always the longest dimension
     differ=rightcells(N,g_grind.space.bc{2},g_grind.space.bcvalue(2))-2*N+leftcells(N,g_grind.space.bc{1},g_grind.space.bcvalue(1));
else
     differ=upcells(N,g_grind.space.bc{3},g_grind.space.bcvalue(3))-2*N+downcells(N,g_grind.space.bc{4},g_grind.space.bcvalue(4));
end

