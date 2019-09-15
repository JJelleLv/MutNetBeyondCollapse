%STABIL   Stabilize
%   Run the current model without showing results and keep the final state as
%   initial value
%
%   Usage:
%   STABIL - runs for 1000 days.
%   STABIL NDAYS - runs for NDAYS days.
%   STABIL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'backwards' [logical] - run backwards? (default=false)
%     'keep' [logical] - keep last value? (default=true)   
%     'ndays' [number>0] - number of days to stabilize
%     'silent' [logical] - suppress messages
%
%   STABIL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-s' - suppress messages
%
%   See also ru, ke
%
%   Reference page in Help browser:
%      <a href="matlab:commands('stabil')">commands stabil</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function stabil(varargin)
global t g_grind g_Y;
fieldnams={'ndays', 'n>0', 'number of days to stabilize',1000;...
   'silent', 'l', 'suppress messages',false;...
   'backwards', 'l', 'run backwards? (default=false)',false;...
   'keep', 'l', 'keep last value? (default=true)',true}';
args=i_parseargs(fieldnams,'ndays,silent',{'-s'},varargin);
if ~isfield(args,'ndays')
   args.ndays = 1000;
end
if any(strcmp('-s',args.opts))
   args.silent=true;
elseif ~isfield(args,'silent')
   args.silent=false;
end
if ~isfield(args,'keep')
    args.keep=true;
end
i_parcheck;
N0 = i_initvar;
oldstep = g_grind.tstep;
try
   if isfield(args,'backwards')&&args.backwards
       g_grind.solver.backwards=true;
       g_grind.tstep=oldstep;
       g_grind.solver.opt.OutputFcn = @i_odetrunc;
       g_grind.truncate=true;
   else
       if ~(g_grind.solver.haslags)
           g_grind.tstep = 2;
       else
           g_grind.tstep=NaN;
       end
       g_grind.solver.opt.OutputFcn = [];
   end
   i_ru(t, args.ndays, N0, 1);
   g_grind.tstep = oldstep;
   if isfield(args,'backwards')&&args.backwards
       g_grind.solver.backwards=false;
       g_grind.truncate=false;
       if args.keep
          ndx=find(~isnan(g_Y(:,1)),1);
          i_keep(g_Y(ndx,:).');
       end
   elseif args.keep
       ke;
   end
   if ~args.silent
      fprintf('Simulated %d days.\n', args.ndays);
   end   
catch err
%   err=lasterror;
   g_grind.tstep = oldstep;
   rethrow(err);
end
