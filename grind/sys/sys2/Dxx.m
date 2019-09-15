%Dxx   Double partial derivative to spatial x (=column) coordinate
%   Defines the double gradient in spatial coordinates [d2N/dx2] to be used for example for diffusion. By default
%   it uses the simple diffusion discretisation:
%   Dxx = Nj-1 - 2*Nj + Nj+1
%   Different boundary conditions can be selected using the command <a href="matlab:help definespace">definespace</a>.
%
%  
%   Usage:
%   Dxx(N) defines in a model a double partial derivative of the variable N with respect to 
%     the spatial x (=column) coordinates (except for single row variables, then x = row).
%
%   See also model, Dyy, Dx, Dy, neighborcells, definespace
%
%   Reference page in Help browser:
%      <a href="matlab:commands('Dxx')">commands Dxx</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function differ=Dxx(N)
global g_grind;
if ~isfield(g_grind,'space')||isempty(g_grind.space)
    warning('grind:Dxx','No space defined, set to defaults');
    definespace('-d');
end

siz=size(N);
if siz(2)==1&&siz(1)~=1 %if vector x is always the longest dimension
     differ=upcells(N,g_grind.space.bc{3},g_grind.space.bcvalue(3))-2*N+downcells(N,g_grind.space.bc{4},g_grind.space.bcvalue(4));
else
     differ=rightcells(N,g_grind.space.bc{2},g_grind.space.bcvalue(2))-2*N+leftcells(N,g_grind.space.bc{1},g_grind.space.bcvalue(1));
 end

