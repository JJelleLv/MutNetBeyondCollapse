%OUT   Select variables for time plots
%   Specify the variables/functions for the time plots. You can define different
%   time plots that are updated simultaneously.
%
%   Usage: 
%   OUT - opens a <a href="matlab:commands outdlg">dialog box</a> to add variables to time plots 
%   OUT FUN1 FUN2 FUN3 - plots the variables FUN1,FUN2 and FUN3 (replaces current
%   list). FUN1 may be a state variable, a parameter, an auxiliary variable
%   (see FUNCS), or an MATLAB expression with these variables;
%   OUT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'fun' [equation or shortcut] - list of variables-functions for a time plot
%     'no' [integer>0] - the number of the time plot (see '-1')
%     'silent' [logical] - silent mode, no output
%     'time' [string] - the variable for the x-axis.
%   OUT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-1 or -2 or -value' - define the variable/functions for time figure number value (default value=1, can be in 
%    combination with other options, but this should then be the first option)
%     '-?' - show current settings
%     '-add' FUN1 FUN2 - adds the    FUN1 FUN2 to the current list
%     '-ch' - remove numbers of time plots that are empty
%     '-cl' - clear all time plots. To remove one plot, use in combination with the first option.
%     '-cleantimevars' - used internally to remove empty plots
%     '-remove' FUN1 FUN2 - remove the variables/functions FUN1 FUN2 from the current list
%     '-silent' - do no show the new settings
%     '-time' FUN - replaces the time axis with the function FUN
%
%   Short hands (can be used in combinations and with previous options)
%     -defaults, all state variables (including data) and external variables.
%     -data, all state variables for with data is defined.
%     -extern, all external variables.
%     -sum, sets the sum of all vector state variables.
%
%   Examples:
%   OUT A  - only the state variable A on the axis
%   OUT -add cons - add the intermediate variable cons
%   OUT -2 A - make a second plot with state variable A
%   OUT -3 -mean -sum -max - mean sum and maximum of vector state variables in plot 3.
%   OUT -add A/K - add an expression with parameters and state variables
%   OUT -time mod(t,365) - replaces the time axis with day number
%
%   See also time, funcs, outf, timesens
%
%   Reference page in Help browser:
%      <a href="matlab:commands('out')">commands out</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function out(varargin)
global g_grind;
if isempty(g_grind)
    error('GRIND:out:NoModel','No model selected');
end
if nargin==0
    %  outdlg;
    i_outdlg;
    return;
end
fieldnams={'fun', 'U1', 'list of variables-functions for a time plot','';...
   'silent', 'l', 'silent mode, no output',false;...
   'time', 's', 'the variable for the x-axis.','t';...
   'no', 'i>0', 'the number of the time plot (see ''-1'')',1}';
args=i_parseargs(fieldnams,'fun(+)','-remove,-?,-add,-silent,-cleantimevars,-cl,-ch,-time,(-[1-9]+[0-9]*)',varargin,false,{@i_is_out_equation});

if ~isfield(args,'silent')
    args.silent = 0;
end
No = 1;
hasNo = 0;
if isfield(args,'no')
    hasNo=1;
    No=args.no;
end
if any(strcmp(args.opts, '-silent'))
    args.silent=1;
end
numopts=~cellfun('isempty', regexp(args.opts,'(-[1-9]+[0-9]*)','match','once'));
if any(numopts)
    No = str2double(args.opts{numopts});
    if isempty(No)
        No = 1;
    else
        hasNo = 1;
        No = -No;
        args.no=No;
    end
end
if any(strcmp(args.opts, '-ch'))
    if hasNo
        g_grind.timevars{No} = checklst(g_grind.timevars{No});
    else
        for iNo = length(g_grind.timevars):-1:1
            if isempty(g_grind.timevars{iNo})
                g_grind.timevars =  {g_grind.timevars{1:iNo - 1} g_grind.timevars{iNo + 1:end}};
                g_grind.outt =  {g_grind.outt{1:iNo - 1} g_grind.outt{iNo + 1:end}};
            end
        end
        for iNo = 1:length(g_grind.timevars)
            g_grind.timevars{iNo} = checklst(g_grind.timevars{iNo});
        end
    end
end
if any(strcmp(args.opts, '-cleantimevars'))
    ndx=cellfun('isempty',g_grind.timevars);
    g_grind.timevars=g_grind.timevars(~ndx);
    g_grind.outt=g_grind.outt(~ndx);
end
if any(strcmp(args.opts, '-cl'))
    if hasNo
        g_grind.timevars{No} = {};
        g_grind.outt{No} = '';
        if ~args.silent
            fprintf('Time plot %d cleared\n', No);
        end
    else
        g_grind.timevars = {};
        g_grind.outt = {};
    end
end

if any(strcmp(args.opts, '-remove'))
    if ~isfield(args,'no')
       No=1:length(g_grind.timevars);
    end
    varlist = getvarlist(args.fun);
    for k=length(No):-1:1
        for j = 1:length(varlist)
            if strcmp(g_grind.outt{No(k)}, varlist)
                g_grind.timevars(No(k))=[];
                g_grind.outt(No(k))=[];
            elseif ~isempty(g_grind.timevars{No(k)})
               varndx=~strcmp(g_grind.timevars{No(k)}, varlist{j});
               g_grind.timevars{No(k)}=g_grind.timevars{No(k)}(varndx);
            end
        end
    end
    if ~args.silent
        out('-?');
    end
    return;
end
% if any(strcmp(args.opts, '-remove'))&&any(strcmp(args.opts, '-data'))
%     args.opts{end+1}= '-removedata';
%     ntime = size(g_grind.timevars{No}, 2);
%     varlist = getvarlist(args.fun);
%     for j = 1:length(varlist)
%         i = 1;
%         while (i <= ntime) && ~strcmp(g_grind.timevars{No}{i}, varlist{j})
%             i = i + 1;
%         end
%         [g_grind.timevars{No},ntime] = removei(g_grind.timevars{No},i, ntime);
%     end
% end
% if any(strcmp(args.opts, '-removeextern'))
%     out('-silent',num2str(-No),'-remove','-extern');
% end
% if any(strcmp(args.opts, '-removedata'))
%     ntime = size(g_grind.timevars{No}, 2);
%     for k = length(g_grind.timevars{No}):-1:1
%         if length(g_grind.timevars{No}{k}) >= 9
%             if strcmpi(g_grind.timevars{No}{k}(1:9), 'observed ')||~isempty(strfind(g_grind.timevars{No}{k},'observed('))
%                 [g_grind.timevars{No},ntime] = removei(g_grind.timevars{No},k, ntime);
%             end
%         end
%     end
% end
if any(strcmp(args.opts, '-add'))
    addtovars(getvarlist(args.fun), No);
    if ~args.silent
        out('-?');
    end
    return;
end
if any(strcmp(args.opts, '-time'))
    if isfield(args,'fun')&&~isempty(args.fun)
        args.time=args.fun{1};
        args.fun={};
    else
        args.time='t';
    end
end
if isfield(args,'time')
    g_grind.outt{No} = args.time;
end
if isfield(args,'fun')&&~isempty(args.fun)
    g_grind.timevars{No} =getvarlist(args.fun);
end
if g_grind.statevars.vector
    done = 0;
    timvars = {};
    if ~isempty(g_grind.timevars)
        for i = size(g_grind.timevars{No}, 2):-1:1
            if (i <= length(g_grind.timevars{No})) && isempty(strfind(g_grind.timevars{No}{i}, '('))
                for j = 1:size(g_grind.statevars.vectnames, 2)
                    if (i <= length(g_grind.timevars{No})) && strcmp(g_grind.timevars{No}{i}, g_grind.statevars.vectnames{j})
                        if g_grind.statevars.dims{j}.dim2 == 1
                            done = 1;
                            g_grind.timevars{No}(i) = [];% removei(g_grind.timevars{No}, i, length(g_grind.timevars{No}));
                            for k = 1:g_grind.statevars.dims{j}.dim1
                                addtimvars = {sprintf('%s(%d)', g_grind.statevars.vectnames{j}, k)};
                                if isempty(timvars)
                                    timvars = addtimvars;
                                else
                                    timvars = [timvars, addtimvars];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if done
        g_grind.timevars{No} = [g_grind.timevars{No}, timvars];
    end
end
if (length(varargin) ~= 1) || ~strcmp(varargin{1}, '?')
    if No <= length(g_grind.timevars)&&~isempty(g_grind.timevars{No})
        fsub=find(~cellfun(@isempty,regexp(g_grind.timevars{No},'(?<![A-Za-z0-9])_[A-Za-z]*')));
        for i = 1:length(fsub)
            g_grind.timevars{No}{fsub(i)} = outf('changeshortcut', g_grind.timevars{No}{fsub(i)});
        end
        if No > length(g_grind.outt) || isempty(g_grind.outt{No})
            g_grind.outt{No} = 't';
        end
    end
end
if ~args.silent
    if isempty(g_grind.timevars)
        disp('No time variables');
    end
    for iNo = 1:length(g_grind.timevars)
        if ~isempty(g_grind.timevars{iNo}) && (~hasNo || (iNo == No))
            if length(g_grind.timevars{iNo}) > 10
                s1 = i_cell2str({g_grind.timevars{iNo}{1:9} '...' g_grind.timevars{iNo}{end}});
            else
                s1 = i_cell2str(g_grind.timevars{iNo});
            end
            fprintf('Time (%s) plot %d:  %s\n',g_grind.outt{iNo}, iNo, s1);
        end
    end
end

function reslist = checklst(varlist)
global g_data;
k = 1;
reslist = cell(1, length(varlist));
for i = 1:length(varlist)
    OK2add = ~isempty(strtrim(varlist{i}));
    if OK2add
        % check for 'observed ' without external data
        if strncmpi(varlist{i}, 'observed ',9)
            if ~isempty(g_data) && ~isempty(g_data.obs)
                avar = strtrim(varlist{i}(9:end));
                iX = i_getno(avar);
                OK2add = iX.isvar & ~min(isnan(g_data.obs(:, iX.no)));
            else
                OK2add = 0;
            end
        end
        %check for double entries
        if OK2add
            for j = 1:i - 1
                if strcmp(varlist{j}, varlist{i})
                    OK2add = 0;
                    break;
                end
            end
        end
    end
    if OK2add
        reslist{k} = varlist{i};
        k = k + 1;
    end
end
if k > 1
    reslist = reslist(1:k - 1);
else
    reslist = {};
end

function  varlist1 = getvarlist(varlist)
global g_grind;
[varlist1, addition] = i_getoutlist(varlist);
k = 1;
for i = 1:length(addition)
    if strncmpi(addition{i}, '-def',4)
        if g_grind.statevars.vector
            for h = 1:length(g_grind.statevars.vectnames)
                varlist1{k} = sprintf('_mean(%s)', g_grind.statevars.vectnames{h});
                k = k + 1;
            end
        else
            for h = 1:g_grind.statevars.dim
                varlist1{k} = i_statevars_names(h);
                k = k + 1;
            end
            v1 = i_getoutlist( {'-extern'});
            for h = 1:length(v1)
                varlist1{k} =  v1{h};
                k = k + 1;
            end
            v1 = i_getoutlist({'-data'});
            for h = 1:length(v1)
                varlist1{k} =  v1{h};
                k = k + 1;
            end
        end
    end
end




function addtovars(vars, No)
global g_grind;
if length(g_grind.timevars) < No
    g_grind.timevars{No} = {};
    g_grind.outt{No} = 't';
end
oudtimevars = g_grind.timevars{No};
k = length(g_grind.timevars{No}) + 1;
newvars=setdiff(vars,oudtimevars);
for j = 1:length(newvars)
    g_grind.timevars{No}{k} =  newvars{j};
    k = k + 1;
end


function res = isnumstring(s) %#ok<DEFNU>
res = 1;
for i = 1:length(s)
    if ~strcontains('0123456789-', s(i))
        res = 0;
        return;
    end
end






