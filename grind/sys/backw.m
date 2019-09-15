%BACKW   Simulate the model backwards
%   Run the model with negated right-hand sides (=backwards). 
%   For difference equations this is more complicated due to the 
%   discrete steps. Therefore a local optimizer is used to find 
%   the next step back in time.
%   This way you can find unstable equilibria (see also: <a href="matlab:help findeq">findeq</a>) or a separatrix.
%   (see: <a href="matlab:help perturb">perturb</a>).
%
%   Usage:
%   BACKW - run backwards with <a href="matlab:commands g_grind.ndays">g_grind.ndays</a> time units (see also <a href="matlab:help simtime">simtime</a> for setting the simulation duration).
%   BACKW NDAYS - run backwards for NDAYS steps
%   BACKW('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'ndays' [number>0] - number of steps to run backwards
%     'truncate' [logical] - stop running after leaving the phase plane.
%   BACKW('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - add a simulation to the current run.
%     '-n' - normal run.
%     '-p1' or '-pa' - show the results of the last run of PARANAL.
%     '-p2' or '-paranal2' - show the results of the last two runs of PARANAL
%     '-r' - force to rerun the model.
%     '-s' - suppress all output (silent mode)
%
%
%   See also ru, simtime, solver, perturb, addmode
%
%   Reference page in Help browser:
%      <a href="matlab:commands('backw')">commands backw</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function backw(varargin)
%(ndays,dotrunc)
global g_grind;
try
    ndays=g_grind.ndays;
catch %#ok<CTCH>
    ndays=1000;
end
fieldnams={'ndays','n>0','number of steps to run backwards',ndays;...
   'truncate','l','stop running after leaving the phase plane.',true}';
args = i_time_options(i_parseargs(fieldnams,'ndays,truncate',...
    {'-r','-a','-n','-s','-p2|-paranal2','-p1|-pa'},varargin));
i_parcheck;
if ~isfield(args,'ndays')
   args.ndays = g_grind.ndays;
end
if ~isfield(args,'truncate')
   args.truncate=1;
end
g_grind.solver.backwards =  ~g_grind.solver.backwards;
if ~isempty(g_grind.solver.opt.Jacobian)
    g_grind.solver.opt.Jacobian=i_getodehandle('Jacobian');
end
oldtrunc=g_grind.truncate;
g_grind.truncate=args.truncate;
try
   ru(args.ndays,args.opts{:})
   g_grind.solver.backwards = ~g_grind.solver.backwards;
   if ~isempty(g_grind.solver.opt.Jacobian)
      g_grind.solver.opt.Jacobian=i_getodehandle('Jacobian');
   end
   g_grind.truncate=oldtrunc;
catch err
%   err=lasterror;
   g_grind.solver.backwards = ~g_grind.solver.backwards;
   if ~isempty(g_grind.solver.opt.Jacobian)
      g_grind.solver.opt.Jacobian=i_getodehandle('Jacobian');
   end
   g_grind.truncate=oldtrunc;
   rethrow(err);
end
