function result = tempsum(fun, critvalue, resets)
if nargin == 0
   fun = 'T';
end
if nargin < 2
   critvalue = 0;
end
if nargin < 3
   resets = [];
end
global g_t;
if ischar(fun)&&~strncmp(fun,'??',2)
   funval = outfun(fun);
   x = g_t;
else
   funval = fun;
   x = 1:length(funval);
end
sumval = max(funval, critvalue) - critvalue;
if length(resets)==1
   resets=1:resets(1):x(end);
end
resets = sort(resets);
ndx=zeros(1,length(resets)+1);
for i = 1:length(resets)
   xx = abs(x - resets(i));
   f=find(xx == min(xx));
   ndx(i) = f(1);
end
ndx(end) = length(funval)+1;
result = nan(length(funval),1);
k=0;
for i = 2:length(ndx)
   if ndx(i - 1) < ndx(i)
      csum= cumsum(sumval(ndx(i - 1):ndx(i) - 1));
      result(k+1:k+length(csum),1) = csum;
      k=k+length(csum);
   end
end

