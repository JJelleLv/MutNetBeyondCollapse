function [v,errmsg]=i_iscolor(x)
errmsg='';
if nargin==0
    v='color code';
    return;
end
v=x;
if ischar(x)
    colorcodes=struct('code',{'b',   'g',   'r',     'c',   'm',    'y',    'k', 'w'},...
    'color',{[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[1 1 0],[0 0 0],[1 1 1]});
    f=find(strcmp(x,{colorcodes(:).code}));
    if ~isempty(f)
        x=colorcodes(f(1)).color;
    else
        try
            x=eval(x);
        catch
            errmsg='argument needs to be a valid color code';
        end
    end
end
if isnumeric(x)
    if length(x)~=3
        errmsg='length of the color code should be 3';
    elseif ~all(x<=1&x>=0)
        errmsg='color code should be 3 values between 0 and 1';
    else
        v=x;
    end
end

