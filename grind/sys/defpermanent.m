%DEFPERMANENT   Define permanent variables. 
%   "Permanent variables" are variables that may change during a run. They are not state variables
%   as they are not part of the differential or difference equation. They offer for instance an easy
%   way to monitor the maximum of a state variable. The initial value can be set in the same 
%   way as with state variables. You can plot the permanent variables in <a href="matlab:help time">time</a> plots etc. Their 
%   values are stored in the g_permanent global variable.
%
%   Usage:
%   DEFPERMANENT VAR - Defines the permaneent variable VAR.
%   DEFPERMANENT VAR VALUE - Adds the initial value VALUE.
%   DEFPERMANENT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     't' [number] - time or time span (used in combination with -p and -initiate)
%     'value' [number] - initial value.
%     'var' [identifier] - name of the permanent variable
%     'varno' [integer>0] - the number of the variable (used in combination with -g or -e)
%   DEFPERMANENT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-activate' - activate the variable
%     '-deactivate' - deactivate the variable
%     '-e' varno - get the end values after a run.
%     '-g' varno - get the output after a run
%     '-initiate' - initiate the permanent variable (called automatically)
%     '-i_makemodel' - used internally
%     '-l' - list the settings.
%     '-nodialog' - suppresses warning dialog for resetting current run
%     '-p' t - get the values at time t
%     '-s' values - sets the values of the permanent variables.
%     '-updatevars' - used internally
%
%
%   See also model, definepars, defextern, outfun, time
%
%   Reference page in Help browser:
%      <a href="matlab:commands('defpermanent')">commands defpermanent</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res = defpermanent(varargin)
%(avar, avalue, apar)
global g_permanent g_grind g_t;
if nargin==0
    error('grind:defpermanent','Not enough parameters')
end
fieldnams={'var', 'U1', 'name of the permanent variable','';...
   'value', 'n', 'initial value.',0;...
   't', 'n', 'time or time span (used in combination with -p and -initiate)',[];...
   'varno', 'i>0', 'the number of the variable (used in combination with -g or -e)',1}';
args=i_parseargs(fieldnams,...
    'if(hasoption(''-initiate,-p'')),deffields=''t'';elseif(hasoption(''-g,-e'')),deffields=''varno'';elseif(hasoption(''-s,-a'')),deffields=''value'';else,deffields=''var,value'';end;','-initiate,-updatevars,-deactivate,-activate,-i_makemodel,-s,-g,-l,-p,-e,-nodialog',varargin,...
    true,{@i_isid});
if strcmp(g_grind.model{1},'%external odefile')
    error('grind:defpermanent:externalode','Cannot use this command when an external odefile is used');
end
if any(strcmp(args.opts, '-initiate'))
    %initiate before run args.t = tspan
    
    %   in  curr_ode.m:
    %   remove defpermanent lines
    %   replace "avar" by "g_grind.permanent{avarno}.currvalue" with i_changevar
    %   update all at the end
    %   i_updatepermanent(at);
    %
    tspan = args.t;
    offs = 1;
    for i = 1:length(g_grind.permanent)
        g_grind.permanent{i}.currvalue = evalin('base', g_grind.permanent{i}.name);
        if isempty(g_grind.permanent{i}.currvalue)
            g_grind.permanent{i}.dims = [1,1];
        else
            g_grind.permanent{i}.dims = size(g_grind.permanent{i}.currvalue);
        end
        g_grind.permanent{i}.from = offs;
        offs = offs + prod(g_grind.permanent{i}.dims);
        g_grind.permanent{i}.to = offs - 1;
    end
    if length(tspan) == 2
        tspan(2)=tspan(1)+tspan(2);
        if isnan(g_grind.tstep)
            tspan = tspan(1):tspan(2);
        else
            tspan = tspan(1):(tspan(2) - tspan(1)) / g_grind.tstep:tspan(2);
        end
    end
    g_permanent.t = tspan(:);
    g_permanent.Y = zeros(length(g_permanent.t), offs - 1) + NaN;
    g_permanent.nextt = 1;
    g_permanent.active = 1;
    g_permanent.lasti = 1;
    g_permanent.lastt(g_permanent.lasti) = tspan(1) - 1;
    g_permanent.lastY=cell(10,1);
    i_updatepermanent(tspan(1)); %set t = 0 to initial values
elseif any(strcmp(args.opts, '-updatevars'))
    for i = 1:length(g_grind.permanent)
        g_grind.permanent{i}.currvalue = evalin('base', g_grind.permanent{i}.name);
    end
elseif any(strcmp(args.opts, '-deactivate'))
    if ~isempty(g_permanent)&& g_permanent.active&&(g_permanent.nextt <= length(g_permanent.t))
        i_updatepermanent;
        g_permanent.nextt = g_permanent.nextt - 1;
    end
    g_permanent.active = 0;
elseif any(strcmp(args.opts, '-activate'))
    g_permanent.active = 1;
    if isfield(args,'value')&&~isempty(args.value)
        defpermanent('-s', args.value);
    end
elseif any(strcmp(args.opts, '-g')) %get the value of variable args.varno
    %output
    defpermanent('-deactivate');
    if ~isfield(args,'varno')
        args.varno=[];
    end
    if isempty(g_permanent) || isempty(g_permanent.Y) || (~isempty(args.varno)&&(args.varno > length(g_grind.permanent)))
        res = [];
    else
        if isempty(args.varno)
            res = g_permanent.Y;
        else
            res = g_permanent.Y(:, g_grind.permanent{args.varno}.from:g_grind.permanent{args.varno}.to);
        end
        if length(g_t) ~= length(g_permanent.t)
            ndx=~isnan(sum(res,2));
            sumndx=sum(ndx);
            if sumndx>2
               res = interp1(g_permanent.t(ndx), res(ndx,:), g_t);
            elseif sumndx==1
                res=res(ndx,:);
            else
                res=res(1,:);
            end
        end
    end
elseif any(strcmp(args.opts, '-l')) %-list
    disp('');
    maxparlen = par('-maxparlen');
    if ~isfield(g_grind,'permanent')||isempty(g_grind.permanent)
        disp('No permanent variables');
    elseif isempty(g_permanent.Y)
        disp('Initial values of permanent variables:')
        s=['%-' num2str(maxparlen) 's = %0.6g\n'];
        for i = 1:length(g_grind.permanent)
            p = evalin('base', g_grind.permanent{i}.name);
            fprintf(s, g_grind.permanent{i}.name, p(:));
        end
    else
        disp('Initial and final values of permanent variables in last run:')
        s=['%-' num2str(maxparlen) 's = %0.5g / %0.5g\n'];
        for i = 1:length(g_grind.permanent)
            p = evalin('base', g_grind.permanent{i}.name);
            fprintf(s, g_grind.permanent{i}.name, p(:), defpermanent('-e', i));
        end
    end
elseif any(strcmp(args.opts, '-p'))
    %get the values at t=at
    if ~isfield(args,'t')
        args.t=[];
    end
    if isempty(g_grind.permanent)
        res = [];
    elseif isempty(args.t) %defpermanent('-p') gets the initial values
   %     defpermanent('-deactivate');
        if isfield(g_grind.permanent{end},'to')
            res = zeros(g_grind.permanent{end}.to, 1);
            for i = 1:length(g_grind.permanent)
                res(g_grind.permanent{i}.from:g_grind.permanent{i}.to) = g_grind.permanent{i}.currvalue(:);
            end
        else
            res=zeros(size(g_grind.permanent));
            for i = 1:length(g_grind.permanent)
                if ischar(g_grind.permanent{i}.currvalue)
                    curr=evalin('base',g_grind.permanent{i}.currvalue);
                else
                    curr=g_grind.permanent{i}.currvalue;
                end
                if numel(curr)==1
                    res(i) = curr;
                end
            end
        end
    elseif ~(isempty(g_permanent) || isempty(g_permanent.Y))
        f=find(g_permanent.t >= args.t);
        if isempty(f)
            args.t = length(g_permanent);
        else
            args.t = f(1);
            if args.t == g_permanent.nextt %add the last values without increasing nextt
                for i = 1:length(g_grind.permanent)
                    g_permanent.Y(args.t, g_grind.permanent{i}.from:g_grind.permanent{i}.to) = transpose(g_grind.permanent{i}.currvalue(:));
                end
            end
        end
        if args.t >= 1
            res = transpose(g_permanent.Y(args.t, :));
        else
            res = [];
        end
    else
        res=[];
    end
elseif any(strcmp(args.opts, '-s'))
    %set the current values to avalue
    if ~isfield(args,'value')||isempty(args.value)
        args.value = defpermanent('-p', length(g_permanent.t));
    end
    if  ~isempty(args.value)
        offset=1;
        for i = 1:length(g_grind.permanent)
            if isfield(g_grind.permanent{end},'to')
                p = args.value(g_grind.permanent{i}.from:g_grind.permanent{i}.to);
                p = reshape(p, g_grind.permanent{i}.dims);
                assignin('base', g_grind.permanent{i}.name, p);
            else
                afrom=offset;
                ato=afrom+numel(g_grind.permanent{i}.currvalue)-1;
                siz=size(g_grind.permanent{i}.currvalue);
                p=args.value(afrom:ato);
                offset=offset+numel(p);
                g_grind.permanent{i}.currvalue=reshape(p,siz);
            end
        end
    end
elseif any(strcmp(args.opts, '-e'))
    %last value of variable varno
    defpermanent('-deactivate');
    if isempty(g_permanent) || isempty(g_permanent.Y)|| (isfield(args,'varno')&&args.varno > length(g_grind.permanent))
        res = [];
    elseif ~isfield(args,'varno')||isempty(args.varno)
        res = g_permanent.Y(end,:);
    else
        res = g_permanent.Y(end, g_grind.permanent{args.varno}.from:g_grind.permanent{args.varno}.to);
    end
elseif any(strcmp(args.opts,'-i_mmodel'))
    %define permanent variable (internally used by i_mmodel)
    %defpermanent(avar)
    args.var = strtrim(args.var);
    f =  strfind(args.var, ';');
    if ~isempty(f)
        args.var = args.var(1:f(1) - 1);
    end
    f =  strfind(args.var, ' ');
    if ~isempty(f)
        avalue = str2num(args.var(f(1) + 1:end));  %#ok<ST2NM>
        avalue =  i_checkstr(avalue);
        args.var = args.var(1:f(1) - 1);
        evalin('base',sprintf('global %s g_permanent',args.var));
        assignin('base', args.var, avalue);
    else
        evalin('base',sprintf('global %s g_permanent',args.var));
        assignin('base', args.var, NaN);
    end
    if ~isempty(g_grind.permanent)
        no = length(g_grind.permanent) + 1;
    else
        no = 1;
    end
    g_permanent.Y = [];
    g_permanent.t = [];
    g_permanent.nextt = 1;
    g_permanent.active = 0;
    g_permanent.lastt = -9999+zeros(10,1);
    g_permanent.lastYs = cell(10,1);
    g_permanent.lasti = 1;
    g_grind.permanent{no}.name = args.var;
    if isempty(avalue)
        g_grind.permanent{no}.currvalue = 0;
    else
        g_grind.permanent{no}.currvalue = avalue;
    end
    g_grind.permanent{no}.initiate=1;
    %       ppars = {};
    %        for i = 1:length(g_grind.pars)
    %          if ~strcmp(g_grind.pars{i}, args.var)
    %            ppars = [ppars g_grind.pars(i)];
    %          end
    %         end
    g_grind.pars = g_grind.pars(~strcmp(args.var,g_grind.pars));
    for i=1:length(g_grind.permanent)
        g_grind.funcnames.names(strcmp(g_grind.permanent{i}.name,g_grind.funcnames.names))=[];
    end
elseif any(strcmp(args.opts, '-i_makemodel'))
    %define permanent variable (internally used by i_mmodel)
    %defpermanent(args.var)
    if isfield(args,'value')
     %   avalue = str2num(apar); 
        evalin('base',sprintf('global %s g_permanent',args.var));
        assignin('base', args.var, args.value);
    else
        evalin('base',sprintf('global %s g_permanent',args.var));
        assignin('base', args.var, NaN);
        args.value=[];
    end
    if ~isempty(g_grind.permanent)
        no = length(g_grind.permanent) + 1;
    else
        no = 1;
    end
    g_permanent.Y = [];
    g_permanent.t = [];
    g_permanent.nextt = 1;
    g_permanent.active = 0;
    g_permanent.lastt = -9999+zeros(10,1);
    g_permanent.lastYs = cell(10,1);
    g_permanent.lasti = 1;
    g_grind.permanent{no}.name = args.var;
    if isempty(args.value)
        g_grind.permanent{no}.currvalue = 0;
    else
        g_grind.permanent{no}.currvalue = args.value;
    end
    g_grind.permanent{no}.initiate=1;
    %       ppars = {};
    %        for i = 1:length(g_grind.pars)
    %          if ~strcmp(g_grind.pars{i}, args.var)
    %            ppars = [ppars g_grind.pars(i)];
    %          end
    %         end
    g_grind.pars = g_grind.pars(~strcmp(args.var,g_grind.pars));
    for i=1:length(g_grind.permanent)
        g_grind.funcnames.names(strcmp(g_grind.permanent{i}.name,g_grind.funcnames.names))=[];
    end
else
    i_parcheck;
    for i = 1:length(g_grind.permanent)
        if strcmp(args.var, g_grind.permanent{i}.name)
            fprintf('"%s" is already defined as permanent variable\n', args.var);
            return;
        end
    end
    if ~any(strcmp(args.opts,'-nodialog'))
    ButtonName = questdlg('Defining a permanent variable requires a reset of the current run, OK to continue?', ...
        sprintf('defpermanent %s',args.var),'OK','Cancel','OK') ;
    else
        ButtonName='OK';
    end
    if strcmp(ButtonName, 'OK')
        if ~isfield(args,'value') && evalin('base',sprintf('exist(''%s'',''var'')',args.var))
            args.value = evalin('base', args.var);
        elseif ~isfield(args,'value')
            args.value = [];
        end
        finishgrind;
        g_grind.model=[{strtrim(sprintf('defpermanent %s %g',args.var,args.value))}, g_grind.model];
        i_makemodel(g_grind.model, g_grind.commands, g_grind.inifile, g_grind.scheme);
    end
end

