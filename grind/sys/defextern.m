%DEFEXTERN   Define external variables. 
%   External variables are parameters that are variable in time. With the command <a href="matlab:help setdata">setdata</a> or 
%   <a href="matlab:help loaddata">loaddata</a> the data are entered. This command can be used in the model definition, but can 
%   also be entered as a GRIND command to change a parameter in an external variable. If you want to shift in time, 
%   you can use the function <a href="matlab:help externlag">externlag</a>.
%
%   Usage:
%   DEFEXTERN VAR DEFAULT - Defines the variable VAR. DEFAULT is the default value used for 
%   simulation outside the scope of the data.
%   DEFEXTERN NAME DEFAULT DATA - You can also enter the data matrix directly.
%   DEFEXTERN NAME DEFAULT -cycle - option cycle reuses the data outside the range
%   DEFEXTERN('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'cycle' [logical] - option cycle reuses the data outside the range
%     'data' [number] - You can also enter the data matrix directly.
%     'default' [string or number] - The default value that is used if there is no data.
%     'floor' [logical] - floor, do not interpolate within time steps.
%     'name' [identifier] - Name of the external variable.
%   DEFEXTERN('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - reactivate external variables.
%     '-c' - cycle: option cycle reuses the data outside the range
%     '-d' - deactivate: external variable are (temporarily) considered to be parameters so all 
%    data is neglected.
%     '-f' - floor, do not interpolate within time steps.
%     '-n' - nocycle: outside the range of the data the default value is used (default behaviour).
%
%   See also model, definepars, defpermanent, setdata, loaddata, externlag
%
%   Reference page in Help browser:
%      <a href="matlab:commands('defextern')">commands defextern</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function g_result = defextern(varargin)
%(name, default, data, opt1, opt2)
global g_grind g_Y;
fieldnams={'name','U1','Name of the external variable.','';...
   'default','s#n','The default value that is used if there is no data.',0;...
   'data','n','You can also enter the data matrix directly.',[];...
   'cycle','l','option cycle reuses the data outside the range',true;...
   'floor','l','floor, do not interpolate within time steps.',false}';
args=i_parseargs(fieldnams,'name,default,data','-c,-n,-f,-d,-a',varargin,false,{@i_isid});
if strcmp(g_grind.model{1},'%external odefile')
    error('grind:defextern:externalode','Cannot use this command when an external odefile is used');
end
if any(strcmp(args.opts,'-c'))
    args.cycle=true;
elseif any(strcmp(args.opts,'-n'))
    args.cycle=false;
elseif ~isfield(args,'cycle')
    args.cycle=false;
end
if any(strcmp(args.opts,'-f'))
    args.floor=true;
elseif ~isfield(args,'floor')
    args.floor=false;
end

if nargin < 1
    if ~isfield(g_grind, 'externvars') || (isempty(g_grind.externvars))
        disp('No external variables defined');
    else
        disp('External variables (default values):');
        maxparlen=par('-maxparlen');
        for i = 1:length(g_grind.externvars)
            if isempty(g_grind.externvars{i}.data)
                fprintf(['%-' num2str(maxparlen) 's = %s  [No data] %s\n'], ...
                    g_grind.externvars{i}.name, g_grind.externvars{i}.default,getoption(i));
            else
                fprintf(['%-' num2str(maxparlen) 's = %s  [datasize:%dx%d] %s\n'], ...
                    g_grind.externvars{i}.name, g_grind.externvars{i}.default, ...
                    size(g_grind.externvars{i}.data), getoption(i));
            end
        end
    end
    return;
    % error('Cannot define external variable');
end
if any(strcmp(args.opts,'-d'))
    %deactivate
    for i=1:length(g_grind.externvars)
        g_grind.externvars{i}.options.active=0;
    end
    g_grind.lastsettings={};
    disp('All external variables are now fixed');
    return;
end
if any(strcmp(args.opts,'-a'))
    %activate
    for i=1:length(g_grind.externvars)
        g_grind.externvars{i}.options.active=1;
    end
    g_grind.lastsettings={};
    disp('All external variables are activated');
   return;
end
%     if isempty(strfind(' ',name))&evalin('base',sprintf('exist(''%s'',''var'')',name))
%        default=num2str(evalin('base',name));
%     else
if ~isfield(args,'default')
    args.default = '0';
elseif isnumeric(args.default)
    args.default=mat2str(args.default);
    %     end
end
if isfield(args,'data')&&iscell(args.data)
    args.data=args.data{1};
end
%    if isfield(args,'name')
%    f = strfind(args.name, ' ');
%    f4=strfind(args.name,'[');
%    f3= strfind(args.name,']');
%    if ~isempty(f3)
%        f=f(f<f4(1)|f>f3(1));
%    end
%    if ~isempty(f)
%       f2 = strfind(name, ';');
%       if ~isempty(f2)
%          if ~isempty(f3)
%             f2=f2(f2<f4(1)|f2>f3(1));
%          end
%          if ~isempty(f2)
%            name = name(1:f2(1) - 1);
%          end
%       end
%       f2 = strfind(name, '%');
%       if ~isempty(f2)
%          name = name(1:f2(1) - 1);
%       end
%       n = name;
%       if (length(f) == 1)
%          name = n(f(1) + 1:length(n));
%          default = '0';
%          data = NaN;
%       else
%          name = n(f(1) + 1:f(2) - 1);
%          f1=strfind(name,'''');
%          if length(f1) == 2
%             name = name(f1(1) + 1:f1(2) - 1);
%          end
%          if length(f) == 2
%             default = n(f(2) + 1:length(n));
%             f1=strfind(default,'''');
%             if length(f1) == 2
%                default = default(f1(1) + 1:f1(2) - 1);
%             end
%             data = NaN;
%          else
%             default = n(f(2) + 1:f(3) - 1);
%             if length(f) == 3
%                opt1 = n(f(3) + 1:length(n));
%                opt2 = ' ';
%             else
%                opt1 = n(f(3) + 1:f(4) - 1);
%                opt2 = n(f(4) + 1:length(n));
%             end
%             if strncmpi('-c',opt1, 2)||strncmpi('''-c',opt1, 2)||strncmpi('-c',opt2, 2)||strncmpi('''-c',opt2, 2)
%                cycle = 1;
%                data = NaN;
%             elseif strncmpi('-n',opt1, 2)||strncmpi('''-n',opt1, 2)||strncmpi('-n',opt2, 2)||strncmpi('''-n',opt2, 2)
%                cycle = 0;
%                data = NaN;
%             end
%             if strncmpi('-f',opt1, 2)||strncmpi('''-f',opt1, 2)||strncmpi('-f',opt2, 2)||strncmpi('''-f',opt2, 2)
%                tofloor = 1;
%                data = NaN;
%             end
%          end
%       end
%    else
%  %     if nargin < 2
%  %       default = '0';
%   %    end
%       if iscell(default)
%          default = default{1};
%       end
%       if ~ischar(default)
%          default = mat2str(default);
%       end
%       if (nargin == 3)
%          if ischar(data)
%             if strncmpi('-c',data, 2)
%                cycle = 1;
%                data = NaN;
%             elseif strncmpi('-n',data,  2)
%                cycle = 0;
%                data = NaN;
%             elseif strncmpi('-f',data,  2)
%                tofloor = 1;
%                data = NaN;
%             else
%                data = eval(data);
%             end
%          end
%       elseif nargin<3
%          data = NaN;
%       end
%       if (nargin == 4)
%          if ischar(opt1)
%             if strncmpi(opt1, '-c', 2)
%                cycle = 1;
%             elseif strncmpi(opt1, '-n', 2)
%                cycle = 0;
%             elseif strncmpi(opt1, '-f', 2)
%                tofloor = 1;
%             end
%          end
%       end
%       if (nargin == 5)
%          if ischar(opt2)
%             if strncmpi(opt2, '-c', 2)
%                cycle = 1;
%             elseif strncmpi(opt2, '-n', 2)
%                cycle = 0;
%             elseif strncmpi(opt2, '-f', 2)
%                tofloor = 1;
%             end
%          end
%       end
%
%    end
% end
ivar = 0;
for i = 1:length(g_grind.externvars)
    if strcmp(args.name, g_grind.externvars{i}.name)
        ivar = i;
    end
end
if ivar == 0
    ivar = length(g_grind.externvars) + 1;
end
ndx= ~strcmp(args.name,g_grind.pars);
ppars=g_grind.pars(ndx);
%if length(g_grind.pars) > length(ppars)
eval(['global ' args.name]);
eval(['clear ' args.name]);
if length(args.default)<20000
   evalin('base',sprintf('%s=%s;',args.name,args.default));
end
g_grind.pars = ppars;
g_grind.externvars{ivar}.name = args.name;
g_grind.externvars{ivar}.default = args.default;
n=i_checkstr(args.default);
g_grind.externvars{ivar}.dim1 = size(n,1);
g_grind.externvars{ivar}.dim2 = size(n,2);
g_grind.externvars{ivar}.options.cycle = args.cycle;
g_grind.externvars{ivar}.options.tofloor = args.floor;
g_grind.externvars{ivar}.options.active =1;
if isfield(args,'data')
    if ~isfield(g_grind.externvars{ivar}, 'data') || xor(isempty(g_grind.externvars{ivar}.data), isempty(args.data)) ...
            || (min(size(g_grind.externvars{ivar}.data) ~= size(args.data)) == 0) ...
            ||       (min(min(g_grind.externvars{ivar}.data == args.data)) == 0)
        g_Y = [];
    end
    g_grind.externvars{ivar}.data = args.data;
else
    if ~isfield(g_grind.externvars{ivar}, 'data')
        g_grind.externvars{ivar}.data = [];
    end
end

modelchange = 1;
opts=sprintf(' %s',args.opts{:});
%remove comments
g_model=regexp(g_grind.model, '^[^%]*','match','once');

ndx1=~cellfun('isempty',regexp(g_model,'\<defextern ','once'));
ndx2=~cellfun('isempty',regexp(g_model,'\<defextern(','once'));
for i = 1:length(g_grind.model)
    %parses defextern name value;%ff
    %or defextern(name,value);%DD
    if ndx1(i)
        d1=parsed_equation(g_grind.model{i}).fs;
         %d1=regexp(g_grind.model{i},'(?<!%.*)[ ;]','split');
    elseif ndx2(i)
         d1=parsed_equation(g_grind.model{i}).fs;
         %d1=regexp(g_grind.model{i},'(?<!%.*)[\(\),]','split');
         %not perfect for nested functions
    else
        d1={};
    end
    d1=d1(~strcmp(d1,' '));
    if length(d1)>1&&strcmpi(d1{1},'defextern')
        if length(d1)==2
            d1{3}='0';
        end
        if strcmp(d1{2},args.name)&&~strcmpi(d1{3},args.default)
            g_grind.model{i} = sprintf('defextern %s %s%s',args.name,args.default,opts);
            i_makemodel;
            g_Y = [];
            if nargout == 1
                g_result= sprintf('%s=externvar(%d,%s,t);',args.name,ivar,args.default);
            end
            return;
        elseif strcmp(d1{2},args.name)
            modelchange = 0;
        end
    end
%     if strncmpi(d2, g_grind.model{i}, length(d2))
%         d3=[d2 ', ''' args.default ''''];
%         if ~strncmpi(d3, g_grind.model{i}, length(d3))
%             g_grind.model{i} = [d3 ');'];
%             i_makemodel;
%             g_Y = [];
%             if nargout == 1
%                 g_result= sprintf('%s=externvar(%d,%s,t);',args.name,ivar,args.default);
%             end
%             return;
%         else
%             modelchange = 0;
%         end
%     end
end
if modelchange
    g_grind.model = [sprintf('defextern %s %s%s',args.name,args.default,opts), g_grind.model];
    starting=isfield(g_grind,'starting');
    i_makemodel;
    if starting
        g_grind.starting=1;
    end
    g_Y = [];
end
if nargout == 1
    g_result= sprintf('%s=externvar(%d,%s,t);',args.name,ivar,args.default);
end
function s=getoption(ivar)
global g_grind;
s='';
if g_grind.externvars{ivar}.options.cycle
    s=[s '-cycle '];
end
if g_grind.externvars{ivar}.options.tofloor
    s=[s '-floor '];
end
