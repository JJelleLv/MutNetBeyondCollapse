%NULL3   3D phase space with nullclines planes
%   Create nullclines in the 3D state space (it is also allowed to have a
%   parameter on one of the axes).
%   Note that there should be state variables or parameters on all 
%   three axes
%
%   Usage:
%   NULL3  - creates nullcline planes using default settings
%   NULL3 NPOINTS - creates nullclines based on a 3D grid of NPOINTS*NPOINTS*NPOINTS points (default=25).
%   NULL3('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'hax' [handle] - handle of the axis (used by replayall)
%     'npoints' [all(integer>=10) and length(integer)<=3] - size of the 3D matrix for calculations
%     'vect' [number] - use extenal data
%  
%
%   See also ax, null, phas
%
%   Reference page in Help browser:
%      <a href="matlab:commands('null3')">commands null3</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function null3(varargin)
%(npoints, opt)
global g_grind;
fieldnams={'npoints', 'all(i>=10)&length(i)<=3', 'size of the 3D matrix for calculations',25;...
   'hax', 'h', 'handle of the axis (used by replayall)',get(get(0,'CurrentFigure'),'CurrentAxes');...
   'vect', 'n', 'use extenal data',[]}';
args=i_parseargs(fieldnams,'npoints', '',varargin);
i_parcheck;
if ~isfield(args,'npoints')||isempty(args.npoints)
    args.npoints = 25;
end
if ~isfield(args,'hax')
    args.hax=get(get(0,'currentfigure'),'currentaxes');
end
hasdata = 0;
if  isfield(args, 'vect')
    %plug in data from a user-defined function (used for R* evaluations)
    Vs = args.vect;
    hasdata = 1;
    args.npoints = size(Vect);
end
if length(args.npoints)==1
    args.npoints=repmat(args.npoints,3,1);
elseif length(args.npoints)==2
    args.npoints=[args.npoints(:);args.npoints(1)];
end
iX = i_getno(g_grind.xaxis.var,true);
iY = i_getno(g_grind.yaxis.var,true);
iZ = i_getno(g_grind.zaxis.var,true);
if ~(iX.isvar || iX.ispar) || ~(iY.isvar||iY.ispar) ||~(iZ.isvar || iZ.ispar)
    ax('-?');
    i_errordlg('Cannot create 3D null-isoclines if there are no state variables on the axes.')
    error('GRIND:null3:NoStatevariables','Cannot create 3D null-isoclines if there are no state variables on the axes');
end
N0 = i_initvar;
Zaxis = g_grind.zaxis.lim;
if abs(Zaxis(1)) < 0.0001
    Zaxis(1) = 0.0001;
end
oldZ = i_getparvar(iZ);
Xs = zeros(args.npoints(1), args.npoints(2), args.npoints(3));
Ys = Xs;
Zs = Xs;
%Vs = Xs;
try
    if ~hasdata
        Vs=zeros(args.npoints(1), args.npoints(2), args.npoints(3),g_grind.statevars.dim);
        i_waitbar(0, args.npoints(3), 'null3', 'Calculating',0.5)
        N0 = i_initvar;
        for z1 = 1:args.npoints(3)
            i_waitbar(1);
            pz = (z1 - 1) * (Zaxis(2) - Zaxis(1)) / args.npoints(3) + Zaxis(1);
            N0 =  i_setparvar(iZ, N0, pz);
            [Xs(:, :, z1), Ys(:, :, z1),  Vs(:,:,  z1,:)] = i_vector(args.npoints(1:2), iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, N0);
            Zs(:, :, z1) = pz + zeros(args.npoints(1),args.npoints(2));
        end
        i_waitbar([]);
    else
        [Xs, Ys, Zs] = meshgrid(linspace(g_grind.xaxis.lim(1), g_grind.xaxis.lim(2), args.npoints(1)), ...
            linspace(g_grind.yaxis.lim(1), g_grind.yaxis.lim(2), args.npoints(2)), ...
            linspace(g_grind.zaxis.lim(1), g_grind.zaxis.lim(2), args.npoints(3)));
    end
    [hfig, new] = i_makefig('phase3');
    if new
        set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
        set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
        hold('on');
    end
    i_plotdefaults;
    set(hfig,'name', '3D phase space');
    
    oldhold = ishold;
    hold('on');
    leg = cell(min(100,size(Vs,4)), 1);
    col = [0 0 1; 1 0 0;0.5 1 0; 0 1 0; 0 0 0.5; 0.5 0 0; 0.5 1 0; 0 0.5 0];
    for i = 1:min(100,size(Vs,4))
        h = patch(isosurface(Xs, Ys, Zs, Vs(:,:,:,i), 0));
        %   camlight('left'); lighting('phong');
        %    alpha(0.85); %opaquenesss
        set(h,'FaceColor', col(mod(i-1, size(col, 1))+1, :), 'EdgeColor','none');
        leg{i} = [i_statevars_names(i) '''' '=0'];
    end
    camlight('left'); lighting('phong');
    material('metal');
    ud = get(gca, 'userdata');
    ud.meta=struct('func','null3','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname',g_grind.yaxis.var,'ylim',g_grind.yaxis.lim,...
        'zname',g_grind.zaxis.var,'zlim',g_grind.zaxis.lim);
    set(gca, 'userdata', ud);
    if new
        legend(leg{:});
    end
    set(gca, 'View', [322.5, 30]);
    box('on');
    xlabel(i_disptext(g_grind.xaxis.var));
    ylabel(i_disptext(g_grind.yaxis.var));
    zlabel(i_disptext(g_grind.zaxis.var));
    set(gca,'XLim', g_grind.xaxis.lim);
    set(gca,'YLim', g_grind.yaxis.lim);
    set(gca,'ZLim', g_grind.zaxis.lim);
    if g_grind.statevars.dim > 3
        htitle = title(i_disptext(['Valid for ' i_othervars(N0, iX.no, iY.no, iZ.no)]));
        set(htitle,'fontweight','normal');
    end
    i_setparvar(iZ, N0, oldZ);
    if oldhold
        hold on;
    end
catch err
    %   err=lasterror;
    N0 = i_setparvar(iZ, N0, oldZ);
    i_keep(N0);
    rethrow(err);
end



