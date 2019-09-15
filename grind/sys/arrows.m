%ARROWS   Draw arrows in a figure to show directions of lines
%   This command works on all simple line plots. The arrows
%   indicate the direction in which a line is drawn.
%   NOTE: in 2D phase planes, arrows can be added manually by pointing with the 
%   mouse and pressing Ctrl-Y. Similarly, arrows can be deleted by pressing Shift-Y in 
%   a phase plane.
%
%   Usage:
%   ARROWS - Adds arrows in the current figure (last series that was added, or 
%   a series that is selected)
%   ARROWS delete - delete the arrows.
%   ARROWS nearest [X, Y] find the closest point from the point [X,Y] in a line 
%   and add an arrow.
%   ARROWS delnearest [X, Y] delete the nearest arrow from the point [X,Y]
%
%
%   See also null, phasclick
%
%   Reference page in Help browser:
%      <a href="matlab:commands('arrows')">commands arrows</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function arrows(hax,narrow, arg2)
nearest = [];
if nargin==1&&ischar(hax)&&strncmp(hax,'??',2)
    return;
end
if nargin < 1||nargin==1&&ishandle(hax)
    if nargin==0
        hax=[];
    end
    narrow = 5;
else
    if ~ishandle(hax)
        if nargin>1
            arg2=narrow;
        end
        narrow=hax;
        hax=[];
    end
end
if isempty(hax)
    hfig = get(0, 'CurrentFigure');
    if isempty(hfig)
        error('GRIND:arrows:NoFigure','No figure available');
    else
        hax = get(hfig, 'CurrentAxes');
    end
else
    hfig=get(hax,'parent');
end
switch narrow
    case {'nearest','nearest_normalized'}
        nearest = i_checkstr(arg2);
        if strcmp(narrow,'nearest_normalized')
            nearest=adjust_to_axis(hax,nearest);
        end
        narrow = 1;
    case 'delete'
        ars=findall(hfig,'tag','Arrow');
        delete(ars);
        return;
    case {'delnearest','delnearest_normalized'}
        ars=findall(hfig,'tag','Arrow');
        if ~isempty(ars)
            nearest = i_checkstr(arg2);
            if strcmp(narrow,'delnearest_normalized')
               nearest=adjust_to_axis(hax,nearest);
            end
            mindist = 1E30;
            iarr = 1;
            for i = 1:length(ars)
                ud = get(ars(i), 'userdata');
                dist = min(sqrt((ud.x(1) - nearest(1))^2 + (ud.y(1) - nearest(2))^2),...
                    sqrt((ud.x(2) - nearest(1))^2 + (ud.y(2) - nearest(2))^2));
                if dist < mindist
                    mindist = dist;
                    iarr = i;
                end
            end
            delete(ars(iarr));
        end
        return;
    otherwise
        narrow = i_checkstr(narrow);
end


iser = 1;
series = findobj(hax,'type','line') ;
if ~isempty(series)
    for i = 1:length(series)
        if strcmp(get(series(i),'selected'),'on')
            iser = i;
        end
    end
    X = get(series(iser), 'xdata');
    Y = get(series(iser), 'ydata');
    Z = get(series(iser), 'zdata');
    n = length(X);
    if n < narrow
        narrow = n;
    end
    if isempty(Z)
        if ~isempty(nearest)
            mindist = 1E30;
            iarr = 1;
            for jser = 1:length(series)
                X = get(series(jser), 'xdata');
                Y = get(series(jser), 'ydata');
                for i = 1:length(X)
                    dist =  sqrt((X(i) - nearest(1))^2 + (Y(i) - nearest(2))^2);
                    if dist < mindist
                        mindist = dist;
                        iarr = i;
                        iser = jser;
                    end
                end
            end
            X = get(series(iser), 'xdata');
            Y = get(series(iser), 'ydata');
            col = get(series(iser), 'color');
            if iarr == size(X)
                iarr = iarr - 1;
            end
            make_arrow(hax,[X(iarr), Y(iarr)], [X(iarr+1), Y(iarr + 1)],col);
        else
            col = get(series(iser), 'color');
            for i = 1:narrow
                iarr = floor((i-1) / narrow * (n - 1)) + 1;
                make_arrow(hax,[X(iarr), Y(iarr)], [X(iarr+1), Y(iarr + 1)],col);
            end
        end
    else
        axis(axis);
        col = get(series(iser), 'color');
        for i = 1:narrow
            iarr = floor((i-1) / narrow * (n - 1)) + 1;
            make_arrow(hax,[X(iarr), Y(iarr), Z(iarr)], [X(iarr+1), Y(iarr + 1), Z(iarr+1)],col);
        end
    end
end
function   point = adjust_to_axis(hax,point)
xlim=get(hax,'xlim');
ylim=get(hax,'ylim');
if numel(point)==2
    point(1)=xlim(1)+(xlim(2)-xlim(1))*point(1);
    point(2)=ylim(1)+(ylim(2)-ylim(1))*point(2);
else
    error('grind:arrow','Z not supported')
end

function make_arrow(hax,X, Y, col)

%h = arrow(X,Y,15, 'BaseAngle', 30);
%set(h, 'facecolor', col);
%set(h, 'edgecolor', col);
%
set(get(hax,'parent'),'CurrentAxes',hax);
axannotation('arrow',[X(1),Y(1)],[X(2),Y(2)],'Color', col, 'Tag','Arrow', 'HeadLength', 12, 'HeadWidth', 9);
%
%A=arrowline([X(1),Y(1)],[X(2),Y(2)],'arrowsize', 500);
%A=struct(A);
%set(A.arrowhead,'FaceColor', [0 0 0], 'EdgeColor',[0 0 0])
%delete(A.line);
%delete(A.fullline);
