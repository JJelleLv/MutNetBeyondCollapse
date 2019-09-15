function h=i_figure(h,varargin)
if nargin==0||isempty(h)
    h=figure;
elseif ~ishandle(h)||strcmp(get(h,'visible'),'on')
    h=figure(h);
else
    set(0,'currentfigure',h);
end

if ~isempty(varargin)
    set(h,varargin{:});
end

