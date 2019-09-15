% RU   Run the model
%    Run the model and update the currently selected plot (the other plots 
%    that need to be updated are closed). If there is no plot selected, a 
%    1D, 2D, 3D phase plane is opened for 1D, 2D and >2D models respectively.
%    This command can also be used to force a new run even if no parameters 
%    were changed.
%
%    Usage:
%    RU - run with the number of days that is in the global variable g_grind.ndays.
%    Use the command <a href="matlab:help simtime">simtime</a> to set the default number of days.
%    RU N - run the model for N days.
%    RU('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'ndays' [number>0] - number of days to run
%    RU('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' -add (-a) continue with the previous run.
%     '-n' -nocheck (-n)  do not check for changed option, just show last results.
%     '-p1' or '-pa' - show the results of the last run of PARANAL.
%     '-p2' or '-paranal2' - show the results of the last two runs of PARANAL
%     '-r' - rerun the model always.
%     '-s' - silent (-s) do not plot the results.
%
%
%    See also simtime, time, phas, null, null3, addmode
%
%   Reference page in Help browser:
%      <a href="matlab:commands('ru')">commands ru</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function ru(varargin)
global t g_grind g_Y g_t;
fieldnams={'ndays', 'n>0', 'number of days to run',1000}';
res = i_time_options(i_parseargs(fieldnams,'ndays',...
    {'-r','-a','-n','-s','-p2|-paranal2','-p1|-pa'},varargin));
i_parcheck;
% if nargin < 2
%    N0 = i_initvar;
% end
i = 1;
while ~ishandle(i) && (i < i_figno('maxno'))
   i = i + 1;
end
if ishandle(i)
   curplot = gcf;
elseif isempty(g_grind.yaxis.var)
   curplot = i_figno('phase1');
   null;
elseif isempty(g_grind.zaxis.var)
   curplot = i_figno('phase2');
   [hfig,fignew]=i_makefig('phase2');
   if fignew
      set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
      set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
      gca;
      i_plotdefaults(hfig);
   end
else
   curplot = i_figno('phase3');
   [hfig,fignew]=i_makefig('phase3');
   if fignew
      set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
      set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
      set(gca, 'View', [322.5, 30]);
      i_plotdefaults(hfig);
   end
end
if (g_grind.drawnow)
   if (curplot == i_figno('phase2'))
      if ~(isempty(i_varno(g_grind.xaxis.var)) || isempty(i_varno(g_grind.yaxis.var)))
         g_grind.solver.opt.OutputFcn = str2func('i_odephas');
      end
   elseif curplot == i_figno('phase3')
      if ~(isempty(i_varno(g_grind.xaxis.var)) || isempty(i_varno(g_grind.yaxis.var)) ...
            ||    isempty(i_varno(g_grind.zaxis.var)))
         g_grind.solver.opt.OutputFcn = str2func('i_odephas');
      end
   end
else
   if ishandle(curplot)
      g_grind.solver.opt.OutputFcn = str2func('i_odespeed');
   else
      g_grind.solver.opt.OutputFcn = [];
   end
end
oldpointer = get(curplot, 'pointer');
set(curplot, 'pointer', 'watch');
if ~res.nocheck
%    disp('running')
   i_ru(t, res.ndays, res.N0, 1);
end
if ishandle(curplot)
   set(curplot, 'pointer', oldpointer);
end
g_grind.solver.opt.OutputFcn = [];
if ~res.silent
    i_phas(curplot, 1);
end
if res.adding
   g_grind.solver.addmode = false;
end


%if ~isempty(out1)
%   g_grind.timevars = out1;
%end
if ~isempty(res.OldY)
   g_Y = res.OldY;
   g_t = res.Oldt;
end
