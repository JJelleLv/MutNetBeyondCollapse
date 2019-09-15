%RUNGRID   Create grid of trajectories
%   Generates a grid of trajectories in a 2D phase plane. These trajectories start at 
%   regular intervals in the state space.
%
%   Usage:
%   RUNGRID - creates a 5x5 grid and runs with the default number of days
%   RUNGRID NPOINTS - creates a NPOINTSxNPOINTS grid
%   RUNGRID [NX NY] - creates a NX x NY grid. If one of these arguments are smaller than
%   1, it is not varied. Instead the current initial conditions are used.
%   RUNGRID [NX NY] NDAYS - runs for NDAYS.
%   RUNGRID [NX NY] NDAYS TRUE - runs for NDAYS in two directions.
%   RUNGRID('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'bidirectional' [logical] - if true the model runs backwards and forwards in each point
%     'ndays' [number>0] - number of days to run
%     'npoints' [integer>0 and length(integer)<=2] - size of the grid of initial conditions.
%
%   See also ru, phas, null 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('rungrid')">commands rungrid</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function rungrid(varargin)
%(nx, ny, ndays, twodir)
global g_grind;
fieldnams={'npoints', 'i>0&length(i)<=2', 'size of the grid of initial conditions.',[5 5];...
   'ndays', 'n>0', 'number of days to run',g_grind.ndays;...
   'bidirectional', 'l', 'if true the model runs backwards and forwards in each point',false}';
args=i_parseargs(fieldnams,'npoints,ndays,bidirectional','',varargin);
i_parcheck;
if ~isfield(args,'npoints')||isempty(args.npoints)
    args.npoints = [5,5];
end
if length(args.npoints)==1
    args.npoints =[args.npoints,args.npoints];
end
if ~isfield(args,'ndays')||isempty(args.ndays)
    args.ndays = g_grind.ndays;
end
if ~isfield(args,'bidirectional')||isempty(args.bidirectional)
    args.bidirectional=false;
end
oldndays = g_grind.ndays;
g_grind.ndays = args.ndays;
olddraw = g_grind.drawnow;
g_grind.drawnow = 0;
OldN0 = i_initvar;
N0 = OldN0;
hwait=i_waitbar(0,prod(args.npoints),'rungrid','Simulating',0.5,true,'position',[345,72,270,64]);
ix = i_getno(g_grind.xaxis.var);
iy = i_getno(g_grind.yaxis.var);
if isempty(iy.no) && isempty(ix.no)
    i_errordlg('Can not set initial variables, because there are no state variables on the axes');
else
    if (args.npoints(1) > 1) && (args.npoints(2) > 1)
        for X = g_grind.xaxis.lim(1) + 0.01:(g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) / (args.npoints(1) - 1):g_grind.xaxis.lim(2)
            for Y = g_grind.yaxis.lim(1) + 0.01:(g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) / (args.npoints(2) - 1):g_grind.yaxis.lim(2)
                if ~DoRun(N0, ix, X, iy, Y, args.bidirectional)
                    %               return;
                end
                if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
                    disp('Cancelled by user');
                    break;
                end
            end
            if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
                break;
            end
        end
    elseif (args.npoints(1) <= 1) &&(args.npoints(2)>1)
        X = NaN;
        for Y = g_grind.yaxis.lim(1) + 0.01:(g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) / (args.npoints(2) - 1):g_grind.yaxis.lim(2)
            if ~DoRun(N0, ix, X, iy, Y, args.bidirectional)
                %         return;
            end
            if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
                disp('Cancelled by user');
                break;
            end
        end
    elseif (args.npoints(2) <= 1) &&(args.npoints(1)>1)
        Y = NaN;
        for X = g_grind.xaxis.lim(1) + 0.01:(g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) / (args.npoints(1) - 1):g_grind.xaxis.lim(2)
            DoRun(N0, ix, X, iy, Y, args.bidirectional);
            if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
                disp('Cancelled by user');
                break;
            end
        end
    else
        Y=NaN;
        X=NaN;
        DoRun(N0, ix, X, iy, Y, args.bidirectional);
    end
end
i_keep(OldN0);
i_waitbar([]);
g_grind.drawnow =  olddraw;
g_grind.ndays = oldndays;

function Ok=DoRun(N0, ix, X, iy, Y, twodir)
N0 = setaxisval(N0, ix, X);
N0 = setaxisval(N0, iy, Y);
i_keep(N0);
Ok=1;
ru;
if twodir
    backw;
end
i_waitbar(1,'Running');
ud=get(gcf,'userdata');
if ~isempty(ud)&&isfield(ud,'stop')
    Ok=~ud.stop;
end
drawnow;

function N0 = setaxisval(N0, ix, aval)
global g_grind;
if ~isempty(ix.no) && ~isnan(aval)
    if ix.ispar
        evalin('base', g_grind.pars{ix.no});
        assignin('base', g_grind.pars{ix.no}, aval);
    else
        N0(ix.no) = aval;
    end
end

