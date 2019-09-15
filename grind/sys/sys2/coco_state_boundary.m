% 
%function [data, y] = coco_state_boundary(prob, data, u) %#ok<INUSL>
%Monitor if state variables become negative
%
%
%
%[data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
% prob = coco_add_func(prob, 'coco_state_boundary', @coco_state)boundary, data.ep_eqn, ...
%   'singular/regular?', 'EP', 'uidx', uidx);
% prob = coco_add_event(prob, 'EPS', 'boundary','EP.minX', 0);
function [data,y] = coco_state_boundary(~, data, u, N0ranges)
x = u(data.x_idx);
ndx=~isnan(N0ranges(:,1));
if ~isempty(ndx)&&any(ndx)
    xlow= min(x(ndx)-N0ranges(ndx,1));
else
    xlow=inf;
end

ndx=~isnan(N0ranges(:,2));
if ~isempty(ndx)&&any(ndx)
   xhigh= min(N0ranges(ndx,2)-x(ndx));
else
    xhigh=inf;
end

y=min(xlow,xhigh);
end
