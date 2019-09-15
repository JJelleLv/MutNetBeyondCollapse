%ADDMODE   Sets the simulation mode of the model
%   If the add mode is ON, a new run is appended to the last run, unless the initial
%   conditions change. If add mode is OFF, nothing happens if for instance 
%   TIME is pressed the second time. By default the add mode is OFF.
%
%   Usage:
%   ADDMODE - Toggles add mode ON or OFF. 
%   ADDMODE ON - Set add mode ON.
%   ADDMODE OFF - Set add mode OFF.
%   ADDMODE ON -SILENT - doesn't display anything.
%   ADDMODE('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'silent' [logical] - suppresses messages
%     'state' [logical] - add mode on or off
%   ADDMODE('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - show state of addmode
%     '-r' - reset, resets the current run.
%     '-s' - silent, no output
%
%
%   See also ru, backw, time 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('addmode')">commands addmode</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function addmode(varargin)
%(state, opt)
global g_grind g_Y;
fieldnams={'state','l','add mode on or off',~g_grind.solver.addmode;...
   'silent','l','suppresses messages',false}';
[args,defaults]=i_parseargs(fieldnams,'state,silent','-s,-r,-?',varargin);
i_parcheck;
args=mergestructs(defaults,args);
if any(strcmp(args.opts,'-s'))
    args.silent=true;
end
if ~any(strcmp(args.opts,'-?'))
   g_grind.solver.addmode=args.state;
end
if any(strcmp(args.opts,'-r'))
    g_Y=[];
    if ~args.silent
        disp('GRIND Reset');
    end
    if ~isfield(args,'atate')
        return;
    end
end
if ~args.silent
    if g_grind.solver.addmode
        disp('GRIND add mode is on');
    else
        disp('GRIND add mode is off');
    end
end
