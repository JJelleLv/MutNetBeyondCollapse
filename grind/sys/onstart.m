%ONSTART   Function that should be run once before each simulation
%   If some initiation should be done before each simulation, use onstart. 
%   This function is run before each new simulation. Similar things can be
%   done with <a href="matlab:help setevent">setevent</a>.
%  
%   Usage:
%   ONSTART FUN - run function FUN before each run (FUN may also be a single line of 
%   commands). (if there was already another function the new function is added)
%   ONSTART('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [string] - function/single line of code to be run before starting each run
%   ONSTART('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - clear all commands.
%     '-l' - list the commands.
%  
%   See also setevent, modelpanel
%
%   Reference page in Help browser:
%      <a href="matlab:commands('onstart')">commands onstart</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function onstart(varargin)
global g_grind;
fieldnams={'fun', 's', 'function/single line of code to be run before starting each run',''}';
args=i_parseargs(fieldnams,'fun(+)','-l,-c',varargin);
if nargin == 0
    for j = 1:length(g_grind.onstart.funs)
        evalin('base', g_grind.onstart.funs{j});
    end
    return;
end
if any(strcmp(args.opts, '-c'))
    g_grind.onstart.funs = {};
    return;
end
if any(strcmp(args.opts, '-l'))
    if isempty(g_grind.onstart.funs)
        disp('no onstart commands');
    else
        fprintf('onstart %s\n',g_grind.onstart.funs{:})
    end
    return;
end
if isfield(args,'fun')
    fun=strtrim(sprintf('%s ',args.fun{:}));
    if ~strcontains(fun,';')
        fun=sprintf('%s;',fun);
    end
    i = length(g_grind.onstart.funs) + 1;
    ndx=strcmp(g_grind.onstart.funs, fun);
    if any(ndx)
        i=find(ndx,1);
        warning('MATLAB:onstart:alreadyin','Function "%s" was already in onstart',fun)
    end
    
    g_grind.onstart.funs{i} = fun;
    evalin('base', g_grind.onstart.funs{i});
end
