function [validated_x,errmsg]=i_is_out_equation(x)
%just the the default check for an equation, but if the equation contains
%short cuts (as we can use in out) they are handled
errmsg='';
if nargin==0
    validated_x='equation or shortcut';
    return;
end
if iscell(x)
    for i=1:length(x)
        [~,errmsg]=i_is_out_equation(x{i});
        if ~isempty(errmsg)
            return
        end
    end
    validated_x=x;
    return;
end
x1=i_getoutlist(x); %-mean -std etc is translated (see command "out")
if length(x1)>1||~(~isempty(x1)&&strcmp(x1{1},x))
    validated_x=x;
    return;
end
x=outf('changeshortcut', x);
[validated_x,errmsg]=i_validat(x,'q',{});