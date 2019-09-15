%TORUS   Polar coordinates plot 
%  Create a time plot using polar coordinates. The first two axis of the phase plane.
%
%   Usage:
%   TORUS - create a torus plot with a period of 365 time units. The x axis starts at -1.
%   TORUS PERIOD XSTART - create torus plot with a period of PERIOD steps. The x axis 
%   starts at -XSTART.
%   TORUS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'period' [number>0] - the period of the torus.
%     'tstart' [number] - the x axis starts at -XSTART
%      
%          
%   See also ru, ax    
%
%   Reference page in Help browser:
%      <a href="matlab:commands('torus')">commands torus</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function torus(varargin)
%(period,increment)
global g_t g_Y t g_grind;
fieldnams={'period', 'n>0', 'the period of the torus.',365;...
   'tstart', 'n', 'the x axis starts at -XSTART',1}';
args=i_parseargs(fieldnams,'period,tstart','',varargin);
i_parcheck;
if ~isfield(args,'period')
    args.period=365;
end
if ~isfield(args,'tstart')
    args.tstart=1;
end
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(t, g_grind.ndays, N0, 1);
end
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if ~iX.isvar|| ~iY.isvar
   error('GRIND:torus:NoStatevars','Error: there are no state variables on the axes, use ax to set the first 2 axes');
end
hfig=i_makefig('torus');
hld=ishold;
plot3(sin(g_t/args.period*2*pi).*(args.tstart+g_Y(:,iX.no)),cos(g_t/args.period*2*pi).*(args.tstart+ ...
   g_Y(:,iX.no)),g_Y(:,iY.no));
i_plotdefaults(hfig);
if ~isoctave&&verLessThan('matlab','8.4.0')
    set(gca, 'drawmode','fast');
else
    set(gca,'SortMethod','depth');
end

hold on;
box off;
lims=max(abs([get(gca,'xlim'),get(gca,'ylim')]));
zlim=get(gca,'zlim');
plot3([0;0],[0;0],get(gca,'zlim'),'k')
plot3([0;0],[-lims,lims],[zlim(1),zlim(1)],'k');
plot3([-lims,lims],[0;0],[zlim(1),zlim(1)],'k');
xlabel(['sin(t)*' g_grind.xaxis.var]);
ylabel(['cos(t)*' g_grind.xaxis.var]);
zlabel(g_grind.yaxis.var);
t1=0:0.05:6.5;
plot3(sin(t1)*lims,cos(t1)*lims,zlim(1)*ones(1,length(t1)),'k');
plot3(sin(t1)*args.tstart,cos(t1)*args.tstart,zlim(1)*ones(1,length(t1)),'k');
if ~hld
   hold off;
end
us1.period=args.period;
us1.increment=args.tstart;
set(hfig, 'userdata', us1);
