function tilman
global g_grind;
if ~isfield(g_grind, 'tilman')
    % error('GRIND:tilman:NoCritVal','Tilman critical values not defined');
end
[H, new] = i_makefig('phase2');
if new
    set(H, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
    set(H, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
end
set(H, 'Name', 'Phase plane');
oldhold = ishold;
hold('on');
ud = get(gca, 'userdata');
ud.meta=struct('func','tilman','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname',g_grind.yaxis.var,'ylim',g_grind.yaxis.lim,'zname','');
set(gca, 'userdata', ud);
i_plotdefaults;
set(gca, 'XLim', g_grind.xaxis.lim);
set(gca, 'YLim', g_grind.yaxis.lim);
xlabel(i_disptext(g_grind.xaxis.var));
ylabel(i_disptext(g_grind.yaxis.var));
if isfield(g_grind, 'tilman')
    critvals = evalin('base', g_grind.tilman.critvals);
    rcs = evalin('base', g_grind.tilman.rcs);
    maxax = [g_grind.xaxis.lim(2), g_grind.yaxis.lim(2)];
    plot([critvals(1, 1), critvals(1, 1), maxax(1)], [maxax(2), critvals(2, 1), critvals(2, 1)], 'b');
    if size(critvals,2)>1
        plot([critvals(1, 2), critvals(1, 2), maxax(1)], [maxax(2), critvals(2, 2), critvals(2, 2)], 'r');
  %      legend('R^*(1)','R^*(2)','C(1)','C(2)')
    end
    if size(critvals,2)>2
        plot([critvals(1, 3), critvals(1, 3), maxax(1)], [maxax(2), critvals(2, 3), critvals(2, 3)], 'g');
 %       legend('R^*(1)','R^*(2)','R^*(2)','C(1)','C(2)','C(3)')
    end
    eq1 = transpose(critvals(:, 1));
    if (size(critvals,2)==2) && (size(critvals,1)==2)
        if all(critvals(:, 1) < critvals(:, 2)) || all(critvals(:, 1) > critvals(:, 2))
            eq1 = transpose(critvals(:, 1));
            eq2 = transpose(critvals(:, 2));
        else
            eq1 = max(critvals,[],2);
            eq2 = eq1;
        end
    end
end
ends = (maxax(1) - eq1(1)) .* rcs(1) + eq1(2);
plot([eq1(1), maxax(1)], [eq1(2), ends(1)], 'b:');

%annotation('arrow','units','normalized',[eq1(1),eq1(1)*0.9],[eq1(2),eq1(2)*0.9]);
if size(critvals,2)>1
    ends = (maxax(1) - eq2(1)) .* rcs(2) + eq2(2);
    plot([eq2(1), maxax(1)], [eq2(2), ends(1)], 'r:');
end
if size(critvals,2)>1
    legend('R^*(1)','R^*(2)','C(1)','C(2)')
end
if size(critvals,2)>2
    legend('R^*(1)','R^*(2)','R^*(2)','C(1)','C(2)','C(3)')
end
legend off
if oldhold
    hold('on');
else
    hold('off');
end

% function [critval, rc] = getcritval(iW, iR, S, alim, npoints)
% global g_grind;
% oudS = evalin('base', S);
% N0 = i_initvar;
% N0(g_grind.tilman.res) = 9999;
% N0(g_grind.tilman.spec) = 0;
% R = alim(1):(alim(2) - alim(1))  / npoints:alim(2);
% dR = zeros(1, size(R, 2));
% dW = dR;
% N0(iW) = 1;
% for i = 1:size(X, 2)
%    N0(iR) = R(i);
%    assignin('base', S, R(i));
%    dR(i) = yy(iR);
%    dW(i) = yy(iW);
% end
% assignin('base', S, oudS);
