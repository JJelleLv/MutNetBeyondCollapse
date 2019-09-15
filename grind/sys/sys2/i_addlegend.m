function s = i_addlegend(adds)
[~, ch] = legend;
s = {};
t=findobj(ch,'type','text');
if ~isempty(t)
   s1 = get(t, 'string');
   s = cellstr(char(s1));
end

if ischar(adds)
   s = [s {adds}];
else
   s = [s; adds];
end

