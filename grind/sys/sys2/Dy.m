%Dy   Single partial derivative to spatial y (=row) direction
%   Defines the gradient in spatial coordinates [dN/dy] to be used for example for advection. Three methods are implemented:
%   Center: centered derivative = 1/2*(Nj+1 - Nj-1) (default)
%   Forward: forward derivative =  Nj+1 - Nj
%   Backward: backward derivative = N - Nj-1
%   These methods and different boundary conditions can be selected using the command <a href="matlab:help definespace">definespace</a>.
% 
%   Usage:
%   Dy(N) defines in a model a partial derivative of the variable N with respect to 
%     the spatial y (=row) coordinates (except for single row variables, then y = column).
%  
%
%   See also model, Dx, Dxx, Dyy, neighborcells, definespace 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('Dy')">commands Dy</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function differ=Dy(N)
global g_grind;
if ~isfield(g_grind,'space')||isempty(g_grind.space)
    warning('grind:Dy','No space defined, set to defaults');
    definespace('-d');
end

scheme=find(strcmp(g_grind.space.diffscheme,{'center','forward','backward'}));
siz=size(N);
if siz(2)==1&&siz(1)~=1 %if vector x is always the longest dimension
       differ=leftrightdiff(N,scheme,g_grind.space.bc,g_grind.space.bcvalue);
else
       differ=updowndiff(N,scheme,g_grind.space.bc,g_grind.space.bcvalue);    
end

function differ=leftrightdiff(N,scheme,bc,bcval)
    switch scheme
        case 1
            differ=(rightcells(N,bc{2},bcval(2))-leftcells(N,bc{1},bcval(1)))/2;
        case 2
            differ=rightcells(N,bc{2},bcval(2))-N;
        case 3
            differ=N-leftcells(N,bc{1},bcval(1));
    end

 function differ=updowndiff(N,scheme,bc,bcval)
   switch scheme
        case 1
            differ=(downcells(N,bc{3},bcval(3))-upcells(N,bc{4},bcval(4)))/2;
        case 2
            differ=downcells(N,bc{3},bcval(3))-N;
        case 3
            differ=N-upcells(N,bc{4},bcval(4));
   end
