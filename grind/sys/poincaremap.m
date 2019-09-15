%POINCAREMAP   Construct a Poincare map
%   All trajectories that cross a certain surface in the state variable space
%   are mapped to itself, i.e. the subsequent crossings are mapped.
%   POINCAREMAP analyses the results of the last run. (if there is
%   no last run, it calls <a href="matlab:help time">time</a>). Use <a href="matlab:help ru">ru</a> if you want to
%   update the last run.
%
%   Usage:
%   POINCAREMAP - opens a dialog to enter the variables.
%   POINCAREMAP FUN PLANE- makes a map of variable or function FUN when it crosses the plane PLANE. 
%   PLANE is a function that defines the plane (for instance X=9). Only 
%   increasing trajectories are analysed.
%   POINCAREMAP FUN PLANE 0 - Analyses only decreasing trajectories.
%   POINCAREMAP('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - variable or function for poincaremap.
%     'hax' [handle] - handle to axis, used by replayall
%     'increasing' [logical] - increasing or decreasing trajectories
%     'plane' [string] - function that defines the plane (for instance X=9)
%   POINCAREMAP('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-addreplaydata' - add replaydata (used internally)
%
%  
%   See also poincaresect, lorenzmap, takens, lyapunov, makemap
%
%   Reference page in Help browser:
%      <a href="matlab:commands('poincaremap')">commands poincaremap</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function poincaremap(varargin)
%(avar1, aplane, increasing)
global g_grind g_t;
fieldnams={'fun', 's', 'variable or function for poincaremap.',i_statevars_names(1);...
   'plane', 's', 'function that defines the plane (for instance X=9)',sprintf('%s==0',i_statevars_names(1));...
   'increasing', 'l', 'increasing or decreasing trajectories',true;...
   'hax', 'h', 'handle to axis, used by replayall',get(get(0,'CurrentFigure'),'CurrentAxes')}';
args=i_parseargs(fieldnams,'if(hasoption(''-addreplaydata'')),deffields=''hax'';else,deffields=''fun,plane,increasing'';end;','-addreplaydata',varargin);
i_parcheck;
if any(~isfield(args,{'fun','plane'}))&&~strcmp(args.opt,'-addreplaydata')
   prompt={'Variable for mapping','Formula of the plane (e.g. x==4)','Increasing trajectories?(1/0)'};
   if isfield(g_grind,'poincare')
       answer=g_grind.poincare;
   else
       answer={i_statevars_names(1),sprintf('%s==0',i_statevars_names(1)),'1'};
   end
   answer=inputdlg(prompt,'Poincare Map',1,answer);
   args.fun=answer{1};
   args.plane=answer{2};
   args.increasing=i_checkstr(answer{3});
%   errordlg('Not enough arguments, Usage: poincaremap avar1 avar avalue');
%   error('GRIND:poincaremap:ArgError','Not enough arguments, Usage: poincaremap avar1 avar avalue');
end
if ~isfield(args,'increasing')
    args.increasing=true;
end
if any(strcmp(args.opts,'-addreplaydata'))
    poincaresect(varargin);
    return;
end
ud.avar1 = args.fun;
ud.aplane = args.plane;
ud.increasing = args.increasing;
t_poincare = i_poincare(args.plane, args.increasing);
if isempty(t_poincare)
     error('GRIND:poincaremap:intersection','Not any points intersect the plane');
end
poincar=outfun(args.fun);
poincar = interp1(g_t, poincar, t_poincare, 'pchip');
g_grind.poincare={args.fun,args.plane,int2str(args.increasing)}; % also used by poincaremap
hfig = i_makefig('poinmap');
set(hfig,'name',['Poincare map through plane ' args.plane])
set(hfig, 'userdata', ud);
m = i_lagmap(poincar);
h1=plot(poincar, m, 'k.');
ud = get(gca, 'userdata');
ud.meta=struct('func','poincaremap','xname',args.fun,'yname',args.fun,'zname','');  
set(gca, 'userdata', ud);
ud=get(h1,'userdata');
ud.tdata = t_poincare(2:end);
ud.plane = args.plane;
ud.xdata = poincar;
ud.ydata = m;
ud.zdata = [];
set(h1,'userdata',ud);
i_plotdefaults(hfig);
xlabel(i_disptext([args.fun '(t)']));
ylabel(i_disptext([args.fun '(t+1)']));
poincaresect('-addreplaydata',get(hfig,'CurrentAxes'));


