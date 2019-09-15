%PHAS   Phase space 
%   Create or open a phase space plot, showing the last run. If
%   there is no run, or if parameters/state variables have been
%   changed, the model is run first.
%
%   Usage:
%   PHAS - if there is a variable for the z-axis selected, create
%   a 3D phase space,
%   else if there is a variable for the y-axis, create a 2D phase plane 
%   else create a 1D phase space.
%   PHAS 2 - create a 2D phase plane if there is something on the
%   y-axis, else create
%   a 1D phase space.
%   PHAS 3 - create a 3D phase space
%   PHAS NDAYS - runs for NDAYS days (cannot be 2 or 3)
%   PHAS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'dim' [integer==2|integer==3] - the dimension of the phase plane (2 or 3)
%     'ndays' [number>0] - number of days to run
%   PHAS('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' -add (-a) continue with the previous run.
%     '-n' -nocheck (-n)  do not check for changed option, just show last results.
%     '-p' - makes it possible to replay a paranal
%     '-p1' or '-pa' - show the results of the last run of PARANAL.
%     '-p2' or '-paranal2' - show the results of the last two runs of PARANAL
%     '-r' - rerun the model always.
%     '-s' - silent (-s) do not plot the results.
%     
%
%   See also ax, null, null3, ru
%
%   Reference page in Help browser:
%      <a href="matlab:commands('phas')">commands phas</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function phas(varargin)
global g_grind t g_Y g_t;
fieldnams={'ndays', 'n>0', 'number of days to run',g_grind.ndays;...
   'dim', 'i==2|i==3', 'the dimension of the phase plane (2 or 3)',2}';
res = i_time_options(i_parseargs(fieldnams,'if nargs==1,deffields=''ndays'';else,deffields=''dim,ndays'';end;',...
    {'-r','-a','-n','-s','-p2|-paranal2','-p1|-pa','-p'},varargin));
i_parcheck;
%res = i_time_options(varargin{:});
if (res.ndays==2)||(res.ndays==3)
    nDim = res.ndays;
    res.ndays = g_grind.ndays;
elseif isfield(res,'dim')
    nDim=res.dim;
else
    if isempty(g_grind.zaxis.var)
        nDim = 2;
    else
        nDim = 3;
    end
end

if res.rerun
    %    disp('running');
    i_ru(t, res.ndays, res.N0, 1);
end

if ~res.silent
    if nDim == 3
        if isempty(g_grind.zaxis.var)
            ax('-?');
            i_errordlg('Cannot create 3D plot if there is nothing on the Z axis (see AX)');
            error('GRIND:phas:NoZAxis','Cannot create 3D plot if there is nothing on the Z axis (see AX)');
        end
        [hfig, fignew] = i_makefig('phase3');
        if fignew
            set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
            set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
            hax=gca;
            set(hax, 'View', [322.5, 30]);
            ud = get(hax, 'userdata');
            ud.meta=struct('func','phas','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname',g_grind.yaxis.var,'ylim',g_grind.yaxis.lim,...
                'zname',g_grind.zaxis.var,'zlim',g_grind.zaxis.lim);
            set(hax, 'userdata', ud);
        end
    elseif isempty(g_grind.yaxis.var)
        [hfig, fignew] = i_makefig('phase1');
        if fignew
            set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
            set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
        end
        hax=gca;
        ud = get(hax, 'userdata');
        ud.meta=struct('func','phas','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname','','zname','');
        set(hax, 'userdata', ud);
    else
        [hfig, fignew] = i_makefig('phase2');
        if fignew
            set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
            set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
        end
        hax=gca;
        ud = get(hax, 'userdata');
        ud.meta=struct('func','phas','xname',g_grind.xaxis.var,'xlim',g_grind.xaxis.lim,'yname',g_grind.yaxis.var,'ylim',g_grind.yaxis.lim,...
            'zname','');
        set(hax, 'userdata', ud);
    end
    i_phas(hfig, 0);
end
if res.adding
    g_grind.solver.addmode = false;
end
if ~isempty(res.OldY)
    g_Y = res.OldY;
    g_t = res.Oldt;
end
