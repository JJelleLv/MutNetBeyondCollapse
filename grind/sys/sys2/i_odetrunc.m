function status = i_odetrunc(~, y, flag)
global g_grind;
status = 0;                             % Assume stop button wasn't pushed.
if nargin < 3 || isempty(flag)           % odephas2(t, y)
    if g_grind.truncate
        iX=i_getno(g_grind.xaxis.var);
        if iX.isvar
            status=any(y(iX.no,:)<g_grind.xaxis.lim(1))||any(y(iX.no,:)>g_grind.xaxis.lim(2));
        end
        if ~status
            iY=i_getno(g_grind.yaxis.var);
            if iY.isvar
                status=any(y(iY.no,:)<g_grind.yaxis.lim(1))||any(y(iY.no,:)>g_grind.yaxis.lim(2));
            end
        end
        if ~status
            iZ=i_getno(g_grind.zaxis.var);
            if iZ.isvar
                status=any(y(iZ.no,:)<g_grind.zaxis.lim(1))||any(y(iZ.no,:)>g_grind.zaxis.lim(2));
            end
        end
    end
end  