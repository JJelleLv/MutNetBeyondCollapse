%SETEVENT   Insert a discrete event in the queue
%    This command is to create a discrete event. A continuous model can be 
%    interrupted by an event (=function) which is scheduled in an event 
%    queue. It is not necessary that this command is done in the model definition, 
%    but it can only be placed in the lower parameter panel. The settings are restored.
%    Example of an event:
%   
%    %Example of an event function
%    function nextt=runevent(t);
%    global A;
%    A=A+4;
%    nextt=t+365;
%
%    Set the result of the function (nextt) to NaN if the event should be deleted.
%
%    Usage:
%    SETEVENT - Dialog box where you can add or edit events.
%    SETEVENT EVENT FIRSTT - EVENT is the name of the event function, FIRSTT is the 
%    first time that the event should take place. FISTT can be a parameter, to get
%    it updated before each run enter a string with the name of the parameter.
%    SETEVENT('simpleevent',FIRSTT,COMM,NEXTT) - create a simple event (e.g. change of 
%    state variable or parameter), COMM = string with commands (e.g. 'A=A-1; if A<0.1,
%    A=0.1;end;') NEXTT = time lag between the events.
%    SETEVENT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'event' [string] - name of event function
%     'firstt' [number or string] - first time it will run
%     'pars' [general] - parameters for the event function
%   SETEVENT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - Clear all events.
%     '-d' - disable all events.
%     '-e' - ensable all events.
%     '-l' - List all events.
%
%  
%    See also model, djump
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setevent')">commands setevent</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function setevent(varargin)
%(event, dt,varargin)
if nargin == 0
    i_seteventdlg;
    return;
end
fieldnams={'event', 's', 'name of event function','simpleevent';...
   'firstt', 'n#s', 'first time it will run',0;...
   'pars', '', 'parameters for the event function',{}}';
args=i_parseargs(fieldnams,'event,firstt,pars(+)','-c,-l,-e,-d',varargin);
if any(strcmp(args.opts,'-c'))
    i_setevent('clear');
    return;
end
if any(strcmp(args.opts,'-l'))
    i_setevent('list');
    return;
end
if any(strcmp(args.opts,'-e'))
    i_setevent('enable',true);
    return;
end
if any(strcmp(args.opts,'-d'))
    i_setevent('enable',false);
    return;
end
if ~isfield(args,'firstt')
    args.firstt=0;
end
if ~isfield(args,'event')
    args.event='simpleevent';
end
if ~isfield(args,'pars')
    if strcmp(args.event,'simpleevent')
        error('grind:setevent:simpleevent','Error in simpleevent, no parameters supplied to define the event');
    end
    args.pars={};
end
i_setevent('addnew',args.event,args.firstt,args.pars);

