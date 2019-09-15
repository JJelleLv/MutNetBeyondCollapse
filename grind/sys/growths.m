%GROWTHS   Growth of each state variable
%   Make filled contour plots with the growth of each state
%   variable in the 2D phase plane.
%
%   Usage:
%   GROWTHS - determine the growth in a grid of 40x40 points.
%   GROWTHS NPOINTS - determine the growth in a grid of NPOINTSxNPOINTS points.
%   GROWTHS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'npoints' [integer>3 and length(integer)<=2] - determines the size of the grid
%
%   See also null
%
%   Reference page in Help browser:
%      <a href="matlab:commands('growths')">commands growths</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function growths(varargin)
%(igrid)
global g_grind;
fieldnams={'npoints', 'i>3&length(i)<=2', 'determines the size of the grid',40}';
args=i_parseargs(fieldnams,'npoints','',varargin);
i_parcheck;
if ~isfield(args,'npoints')
    args.npoints=40;
end
iX = i_getno(g_grind.xaxis.var,true);
iY = i_getno(g_grind.yaxis.var,true);
if g_grind.statevars.dim>20
    disp('Too many state variables only the first 20 shown');
end
if ~(iX.isvar||iX.ispar||iX.isext) || ~(iY.isvar||iY.ispar||iX.isext)
    ax('-?');
    error('GRIND:growths:NoStatevars','Cannot create growths plot if there are no state variables on the axes');
end
N0 = i_initvar;
if ~(isempty(iX.no) || isempty(iY.no))
    [X, Y, Vect] = i_vector(args.npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim);
    %  H = figure(i_figno('growths') + 1);
    for i = 1:min(20,g_grind.statevars.dim)
        [hfig,isnew] = i_makefig('growths',i);
        if isnew
            set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
            set(hfig, 'KeypressFcn',@(hobj,ev)i_callb('keypressed',hobj));
            set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
        end
        set(hfig, 'Name', ['Growth of ', char(i_statevars_names(i))]);
        set(hfig, 'doublebuffer','on');
        if ~isoctave&&verLessThan('matlab','8.6')
            set(hfig,'renderer','ZBuffer'); %in fact this is solving a window bug
        end
        V=Vect(:, :, i);
        if ~any(isreal(V(:)))
            V=real(V);
            warning('grind:growth:complex','Imaginary part of solution ignored');
        end
        if all(isnan(V(:)))
            warning('grind:growths:nonumber','No valid numbers of %s to plot',i_statevars_names(i));
            break;
        end
        surf(X, Y, V)
        ud = get(gca, 'userdata');
        ud.meta=struct('func','growths','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname',g_grind.yaxis.var,'ylim',g_grind.yaxis.lim,'zname',sprintf('%s''',i_statevars_names(i)));  
        set(gca, 'userdata', ud);
        i_plotdefaults(hfig);
        oldhold=ishold;
        hold('on');
        contour3(X, Y, V, [0, 0], 'k');
        if ~oldhold
            hold off;
        end
        xlabel(i_disptext(g_grind.xaxis.var));
        ylabel(i_disptext(g_grind.yaxis.var));
        zlabel('growth');
        set(gca, 'XLim', g_grind.xaxis.lim);
        set(gca, 'YLim', g_grind.yaxis.lim);
        %shading flat;
        shading('interp');   % plot a surface
        %       if getrelease<14 % for some reason matlab 7.5 cannot do this lightning
        %           light;
        %           lighting('gouraud');
        %           material('dull');%/shiny/metal
        %       end
        colorbar;
        if g_grind.statevars.dim > 2
            htitle=title(i_disptext(['Valid for ' i_othervars(N0, iX.no, iY.no)]));
            set(htitle,'fontweight','normal');
        end
    end
else
    i_errordlg('No state variable on the y axis');
    error('GRIND:growths:NoStateVar','growths: No state variable on the y axis');
end
