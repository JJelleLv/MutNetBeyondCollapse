function i_plotdefaults(hfig, apen)
global g_grind;
if nargin < 2
    if ~isempty(g_grind) && isfield(g_grind, 'pen')
        apen = g_grind.pen;
    else
        apen = i_nextpen([]);
    end
    
end

if nargin < 1
    hfig = gcf;
end

try
    hax=get(hfig,'CurrentAxes');
catch
    hax=findobj(hfig,'type','axes');
end
if isempty(hax)
    hax = gca;
end

Cbar=findobj(hfig,'tag','Colorbar');
if ~isempty(Cbar)
    set(Cbar, 'FontSize', apen.fontsize);
end

set(hax, 'FontSize', apen.fontsize);
set(hax, 'LineWidth', apen.linewidth);
set(hax, 'TickDir', 'out');
set(hax,'TickLength',[0.015 0.025])
%set(H, 'TickDir', apen.tickdir);
set(hax, 'Box', apen.box);
t=findobj(hax,'type','text');
if ~isempty(t)
    if ~isfield(apen, 'fontsize')
        apen.fontsize = 14;
    end
    
    set(t, 'fontsize', apen.fontsize);
end

t1= [findobj(hfig,'tag','conteq-text');findobj(hfig,'tag','contbif-text');findobj(hfig,'tag','label-text')];
set(t1, 'fontsize',10);
if isempty(Cbar) || ~any(hax == Cbar)
    View = get(hax, 'view');
    if (apen.fontsize >= 14)  && min(View == [0 90])
        set(hax,'units','normalized');
        if ~isempty(Cbar)
            P = [0.1500    0.1900    0.6455    0.7350];
        else
            P = [0.1300    0.1100    0.7750    0.8150];
        end
        
        pos1 = get(hax,'position');
        if (pos1(3)>0.6)&&(pos1(4)>0.6)
            set(hax, 'position', transform(P,apen.fontsize));
            if ~isempty(Cbar)
                P = [0.831375 0.13 0.058125 0.795];
                set(Cbar, 'position', transform(P,apen.fontsize))
            end
            
        end
        
    end
    
end

function P = transform(P, fontsize)
%P = [left, bottom, width, height]
%where left and bottom define the distance from the lower-left corner of the screen
if fontsize <= 18
    P(1) = P(1) + 0.01;
    P(2) = P(2) + 0.03;
    P(3) = P(3) - 0.02;
    P(4) = P(4) - 0.04;
elseif fontsize <= 22
    P(2) = P(2) + 0.04;
    P(4) = P(4) - 0.04;
elseif fontsize <= 28
    P(1) = P(1) + 0.02;
    P(3) = P(3) - 0.02;
    P(2) = P(2) + 0.08;
    P(4) = P(4) - 0.08;
else
    P(1) = P(1) + 0.08;
    P(3) = P(3) - 0.08;
    P(2) = P(2) + 0.16;
    P(4) = P(4) - 0.16;
end

return

