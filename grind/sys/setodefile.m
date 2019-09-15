%SETODEFILE   Use specific ODEFILE
%   use setodefile to use a specific ODEFILE (is stored in <a href="matlab:commands g_grind.odefile">g_grind.odefile</a>). The names of the
%   parameters should be given as extra arguments of the odefile. These names are extracted from this file into 
%   <a href="matlab:commands g_grind.pars">g_grind.pars</a> and are declared as global variables;
%   This command is only used for special models. Parameters need to be initialized separately (or entered in the "commands" field).
%
%   Usage: SETODEFILE FILE = File is set as the new ODEFILE (see also MATLAB's ODEFILE).
%   SETODEFILE FILE {'statevar1','statevar2'} = File loaded and the list of state variables is used.
%   SETODEFILE FILE VAR command -difference = Treat the file as a difference equation.
%   SETODEFILE('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'commands' [cell] - cell array of strings with parameter initializations and default commands
%     'file' [string] - name of the odefile
%     'var' [identifier] - list of state variables
%   SETODEFILE('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-d' - is difference equation
%
%   See also model, use, g_grind.odefile, g_grind.pars 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setodefile')">commands setodefile</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function ggrind=setodefile(varargin)
%(afile, varlist,opt)
global g_grind
addpath([grindpath filesep 'sys2']);
fieldnams={'file', 's', 'name of the odefile','';...
    'var', 'U1', 'list of state variables',''
    'commands','c','cell array of strings with parameter initializations and default commands',''}';
args=i_parseargs(fieldnams,'file,var,commands','-d',varargin,false,{@i_isid});

gg_grind=i_init_g_grind;
if nargin == 0
    help setodefile;
    error('GRIND:setodefile:NoFile','No file indicated');
end
if any(strcmp(args.opts,'-d'))
    gg_grind.solver.isdiffer = 1;
    gg_grind.solver.name = 'i_differ';
end
if isfield(args,'var')
    if ischar(args.var)
        args.var = eval(args.var);
    end
    gg_grind.statevars.names = args.var;
    gg_grind.statevars.dim = length(args.var);
end
f=strfind(args.file,'''');
if ~isempty(f)
    args.file(f)=[];
end
[~, name] = fileparts(args.file);
gg_grind.odefile = name;
gg_grind.pars = cell(1, 50);
name = [name '.m'];
if ~exist(name, 'file')
    cd(findgrindfile(name));
end
p = 1;
ID = fopen(name, 'r');
while 1
    line = fgetl(ID);
    if ~ischar(line), break, end
    if strncmp(line,'function ',9)
        pars=regexp(line,'[,\)\(]','split');
        pars=pars(2:end-1);
        if length(pars)>2
            gg_grind.pars=pars(3:end);
            p=[];
            break
        end
    end
    k = strfind(line, 'global ');
    if ~isempty(k)
        isspace = 1;
        for i = 1:k - 1
            if line(i) ~= ' '
                isspace = 0;
            end
        end
        if isspace
            k = k + 7;
            while line(k) ~= ';'
                while line(k) == ' '
                    k = k + 1;
                end
                k0 = k;
                while (line(k)~=' ')&&(line(k)~=';')
                    k = k + 1;
                end
                apar = line(k0:k - 1);
                if ~strcmp(apar, 'gg_grind')
                    gg_grind.pars{p} = apar;
                    p = p + 1;
                end
            end
        end
    end
end
if ~isempty(p)
    gg_grind.pars = gg_grind.pars(1:p - 1);
end
fclose(ID);
if gg_grind.statevars.dim > 0
    if gg_grind.statevars.vector
        evalin('base', i_globalstr(gg_grind.statevars.vectnames));
    else
        evalin('base', i_globalstr(gg_grind.statevars.names));
    end
end
evalin('base', i_globalstr(gg_grind.pars));
if isfield(args,'commands')
    gg_grind.commands=args.commands;
end
gg_grind.model = {'%external odefile'};
if nargout==0
    %if there are no outputarguments, really load the new model
    finishgrind;
    evalin('base','initgrind');
    g_grind=gg_grind;
    resetpars;
    i_modelinit;
    if g_grind.statevars.dim == 0
        disp(' ');
        disp('>>>>Please enter names of state variables as list');
    end
    i_parcheck(1);
else
    ggrind=gg_grind;
end
