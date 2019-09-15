function out=matcontmodel
global  g_grind;
activepars=sprintf(',%s',g_grind.pars{:});
out{9}=[];
out{1} = @my_init;
out{2} = i_getodehandle('singlestep',activepars);
if ~isempty(g_grind.syms.Jacobian)
    out{3} = i_getodehandle('Jacobian',activepars);
end
if ~isempty(g_grind.syms.Jacobianp)
    out{4} = i_getodehandle('Jacobianp',activepars);
end
if ~isempty(g_grind.syms.Hessian)
    out{5} = i_getodehandle('Hessian',activepars);
end
if ~isempty(g_grind.syms.Hessianp)
    out{6} = i_getodehandle('Hessianp',activepars);
end
if ~isempty(g_grind.syms.der3)
    out{7} = i_getodehandle('der3',activepars);
end
if ~isempty(g_grind.syms.der4)
    out{8} = i_getodehandle('der4',activepars);
end
if ~isempty(g_grind.syms.der5)
    out{9} = i_getodehandle('der5',activepars);
end
end

function [tspan,y0,options] = my_init
global g_grind
han = matcontmodel;
y0=zeros(1,g_grind.statevars.dim);
options = odeset('Jacobian',han(3),'JacobianP',han(4),'Hessians',han(5));
tspan = [0 10];

end


