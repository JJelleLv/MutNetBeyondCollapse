% 
%function [data, y] = coco_test_PD(prob, data, u) %#ok<INUSL>
%test function for Period doubling (flip) bifurcation (co-dim 1, discrete)
%
function [data,y] = coco_test_PD(~, data, u)
x = u(data.ep_eqn.x_idx);
p = u(data.ep_eqn.p_idx);
J=data.dfdxhan(x,p); %jacobian of f(x)-x
y=det(J+2*eye(size(x,1)));
