function g_X2=i_period_map(t,g_X1,varargin)
global g_grind;
per=g_grind.solver.map.tperiod;
if ischar(per)
  per=evalin('base',per);
end

tsample=g_grind.solver.map.tsample;
if ischar(tsample)
   tsample=evalin('base',tsample);
end

ts=tsample+[t;t+per/2;t+per];
g_X2=zeros(size(g_X1));
for i=1:size(g_X1,2) %in case you have vectorized your model this is a loop (not better than non-vectorized)
   [~, Y] = feval(str2func(g_grind.solver.map.solver),i_getodehandle(7,'') , ts(:,i), g_X1(:,i),g_grind.solver.map.opt);
    g_X2(:,i)=transpose(Y(end,:));
end

