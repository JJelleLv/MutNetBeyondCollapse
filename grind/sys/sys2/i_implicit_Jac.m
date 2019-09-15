function [dFdy,dFdyp] = i_implicit_Jac(t,y,yp)
global g_grind;
if ~isfield(g_grind.solver,'Jhandles')||isempty(g_grind.solver.Jhandles)
    s = sprintf('%s,',g_grind.syms.Jacobian_y{:});
%     s(s=='''')='#';
%     parsed = parsed_equation(s);
%     for i = 1:g_grind.statevars.dim
%         parsed.changevar(i_statevars_names(i), sprintf('g_X1(%d)', i));
%     end
% 
%     for i = 1:g_grind.statevars.dim
%         parsed.changevar([i_statevars_names(i) '#'], sprintf('g_X2(%d)', i));
%     end

    g_grind.solver.Jhandles{1}= evalin('base',sprintf('@(t,g_X1,g_X2)reshape([%s],[%d,%d]);',s,g_grind.statevars.dim,g_grind.statevars.dim));
     s = sprintf('%s,',g_grind.syms.Jacobian_yp{:});
%     s(s=='''')='#';
%     parsed = parsed_equation(s);
%     for i = 1:g_grind.statevars.dim
%         parsed.changevar(i_statevars_names(i), sprintf('g_X1(%d)', i));
%     end
% 
%     for i = 1:g_grind.statevars.dim
%         parsed.changevar([i_statevars_names(i) '#'], sprintf('g_X2(%d)', i));
%     end

    g_grind.solver.Jhandles{2}= evalin('base',sprintf('@(t,g_X1,g_X2)reshape([%s],[%d,%d]);',s,g_grind.statevars.dim,g_grind.statevars.dim));
end

f1= g_grind.solver.Jhandles{1};
dFdy=f1(t,y,yp);
f2= g_grind.solver.Jhandles{2};
dFdyp=f2(t,y,yp);
