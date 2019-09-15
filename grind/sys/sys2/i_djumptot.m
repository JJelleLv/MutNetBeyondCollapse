function res=i_djumptot(t,varargin)
global g_grind
res=cell(size(g_grind.solver.djump.args));
for j=1:numel(res);
    res{j}=i_djump(j,t,varargin{j}{:});
end

