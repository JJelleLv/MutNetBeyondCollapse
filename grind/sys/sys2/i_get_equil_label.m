function [lab,eigenval]=i_get_equil_label(N0)
global g_grind
if nargin<1
    N0=i_initvar;
end
han=i_getodehandle(5,'');
isequil=han(N0)<0.001;
[~, eigenval] = i_eigen(isempty(g_grind.syms.Jacobian),g_grind.solver.iters,N0);
[stable, issaddle, isspiral] =  i_stability(eigenval, g_grind.solver.isdiffer);
if ~isequil
    s1='no equilibrium';
elseif issaddle && ~isspiral
    s1 = 'saddle';
else
    if stable
        s1 = 'stable ';
    else
        s1 = 'unstable ';
    end
    if isspiral
        s1 = sprintf('%sspiral', s1);
    else
        s1 = sprintf('%snode', s1);
    end
end
if max(abs(N0)) < 1E-3
    s = '0';
elseif abs(min(N0)) < 1E-3
    s = '+/0';
elseif min(N0) < 0
    s = '--';
else
    s = '++';
end
if any(~isreal(N0))
    eq=round(N0*1000)/1000;
    s2 = sprintf('%g%+gi ', real(eq),imag(eq));
else
    s2 = sprintf('%g ', round(N0 * 10000) / 10000);
end
lab = sprintf('%s (%s): %s', s1,s, s2);
