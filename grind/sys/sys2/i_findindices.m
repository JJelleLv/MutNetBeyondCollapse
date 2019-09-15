function ndx = i_findindices(var, dims)
if length(dims) == 1
   dims = [dims, 1];
end

xxx = zeros(dims);
var=strrep(strrep(var,'}',')'),'{','(');
f = strfind(var, '(');
if isempty(f)
   ndx=find(xx == 0);
else
   eval(sprintf('xxx%s =1;', var(f(1):end)));
   ndx=find(xxx == 1);
end


