function [v,errmsg]=i_isid(x)
errmsg='';
if nargin==0
    v='identifier';
    return;
end
if ischar(x)&&any(x==' ')
    x=regexp(x,' ','split');
end
if ischar(x)&&any(x=='{')
    x=eval(x);
end
v=x;
if iscell(x)
    for i=1:length(x)
       [~,errmsg]=i_isid(x{i});
       if ~isempty(errmsg)
           return;
       end
    end
elseif ~ischar(x)
    errmsg='argument needs to be a valid identifier';
else
    v=x;
    xx=regexp(x,'\<[A-Za-z][A-Za-z0-9_]*','match','once');
    if ~strcmp(xx,x)
        errmsg='argument needs to be a valid identifier';
    end
end
