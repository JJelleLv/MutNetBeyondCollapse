%PAR   Show parameters and model
%   Show the model equations and the values of parameters
%  
%   Usage:
%   PAR - shows the parameters and model
%   PAR MODELONLY - shows only the model
%   PAR('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'output' [string] - determines which part of the model is shown (modelonly | paronly | full)
%     'pargroup' [string] - name of a group of parameters (useful for complex models)
%     'pars' [parameter or variable] - a parameter/state variable/external variable or a list of variables
%     'statevars' [logical] - set to true to get a structure including the state variables (see '-v')
%     'xtra' [general] - extra field for some options
%   PAR('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-d' PARS - get the description of the parameters PARS (if omitted, all parameters are used)
%     '-e' PARS - edit the parameters of the model
%     '-g' PARGROUP PARS - group the parameters PARS in a group named PARGROUP
%     '-la' - shows the equations in LateX
%     '-makeshortref' - used internally
%     '-maxparlen' - used internally to get the maximum length of a parameter
%     '-set' PAR XTRA - change the description of PAR to EXTRA (if there is comma in the description the remaining part
%    is assumed to be a unit)
%     '-setu' PAR XTRA - change the unit of PAR to XTRA
%     '-size' - list the sizes of parameters
%     '-u' - get the units of the parameters PARS (if omitted, all parameters are used)
%     '-v' - get a structure with all values of the parameters (fields are parameter names). If 
%    'statevars'=true the state variables are added to the structure
%     '-vector' - get all parameters in a single vector (needed for MatCont)
%
%   See also gstat
%
%   Reference page in Help browser:
%      <a href="matlab:commands('par')">commands par</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [p,p2] = par(varargin)
global g_grind t;
if nargin>=1&&strcmp(varargin{1},'group')
    %old format group should be replaced by -group (for consistency)
    varargin{1}='-g';
end
fieldnams={'output', 's', 'determines which part of the model is shown (modelonly | paronly | full)','full';...
   'pars', 'U1', 'a parameter/state variable/external variable or a list of variables','';...
   'pargroup', 's', 'name of a group of parameters (useful for complex models)','';...
   'xtra', '', 'extra field for some options',[];...
   'statevars', 'l', 'set to true to get a structure including the state variables (see ''-v'')',false}';
args=i_parseargs(fieldnams, 'if(hasoption(''-d,-u,-set,-setu'')),deffields=''pars,xtra(+)'';elseif(hasoption(''-g'')),deffields=''pargroup,pars(+)'';elseif(hasoption(''-v'')),deffields=''statevars'';else,deffields=''output,xtra'';end;','-size,-setu,-set,-makeshortref,-maxparlen,-vector,-v,-la,-e,-g,-d,-u'...
    ,varargin,false,{@i_isparvar});
if ~isfield(args,'output')
    args.output='';
end
% elseif ~isempty(fullout) && strcmp(fullout, '?')
%    e = par('');
%    [str,ok]=listdlg('liststring',e,'Name','Parameters',...
%       'PromptString','Select parameter','SelectionMode',...
%       'single','uh',40,'ffs',10,'fus',10);
%    if ok
%       f = strfind(e{str}, '=');
%       if ~isempty(f)
%          parname = strtrim(e{str}(1:f(1) - 1));
%          val = strtrim(e{str}(f(1) + 1:length(e{str}) - 1));
%          val=inputdlg(['Change parameter ' parname],'Change parameter',1,{val});
%          assignin('base', parname, str2num(val{1})); 
%       end
%    end
if any(strcmp(args.opts, '-vector'))
    %fast method to get all parameters as one vector
    s = sprintf('%s(:);', g_grind.pars{:});
    p = evalin('base',sprintf('[%s];',s));
    return
end
if any(strcmp(args.opts, '-v')) %values
    p.t = t;
    if ~isempty(g_grind)&&isfield(g_grind,'pars')
        for i = 1:length(g_grind.pars)
            p.(g_grind.pars{i}) =  evalin('base', char(g_grind.pars{i}));
        end
        if isfield(args,'statevars')&&args.statevars
            if g_grind.statevars.vector
                for i = 1:length(g_grind.statevars.vectnames)
                    p.(g_grind.statevars.vectnames{i}) =  evalin('base', char(g_grind.statevars.vectnames{i}));
                end
            else
                for i = 1:length(g_grind.statevars.names)
                    nam = g_grind.statevars.names{i};
                    f = strfind(nam, '(');
                    if ~isempty(f)
                        nam = nam(1:f(1) - 1);
                    end
                    p.(nam) =  evalin('base', char(nam));
                end
            end
        end
    end
    return;
end
if any(strcmp(args.opts, '-la')) %latex
    if ~isfield(g_grind,'model')||isempty(g_grind.model)
        error('grind:par','No model selected');
    end
    mod=g_grind.model;
    for i=1:length(mod)
        f=strfind(mod{i},'%');
        if isempty(f)||f(1)>3
            
            if (length(mod{i})>1)&&(mod{i}(end)==';')
                mod{i}=mod{i}(1:end-1);
            end
            mod{i}=strrep(mod{i},'(t+1)','_#');
            mod{i}=strrep(mod{i},'(t)','_t');
            mod{i}=strrep(mod{i},sprintf('...\n'), '');
            %   mod{i}=strrep(mod{i},'(t)','_{t}');
        end
    end
    str2latex(mod)
    set(gcf,'Name',sprintf('Model equations %s',g_grind.inifile));
    %    ch = get(0, 'children');
    %    latexview=strcmp(get(ch,'tag'),'equationsdlg');
    return;
end
if any(strcmp(args.opts, '-e'))
    %use options -e and -d for debug mode of analunits
    debug=any(strcmp(args.opts, '-d'));
    i_parcheck;
    i_pardlg(debug);
    return;
end
if any(strcmp(args.opts, '-g')) % || strncmpi(args.output, 'group',5)
    %Create/list parameter groups
    if nargin == 1
        if isfield(g_grind, 'pargroups')
            pg = sort(unique(g_grind.pargroups));
            maxlen = 0;
            pp = {};
            for i = 1:length(pg)
                if length(pg{i}) > maxlen
                    maxlen = length(pg{i});
                end
            end
            if length(pg) > 1
                pp = cell(size(pg));
                for i = 1:length(pg)
                    pars = g_grind.pars(strcmp(pg{i}, g_grind.pargroups));
                    pp{i}=sprintf(['par group %' num2str(maxlen+2) 's  %s'],['''' pg{i} ''''],strtrim(sprintf('%s ',pars{:})));
                end
            end
            if nargout >= 1
                p = pg;
                if nargout > 1
                    p2 = pp;
                end
            else
                fprintf('%s\n', pp{:});
            end
        else
            if nargout >= 1
                p = {};
                p2 = {};
            end
        end
        return;
    end
    defgroup = '-';
    if ~isfield(g_grind, 'pargroups')
        g_grind.pargroups = cell(size(g_grind.pars));
        for i = 1:length(g_grind.pargroups)
            g_grind.pargroups{i} = defgroup;
        end
    end
    if isfield(args,'pars')
        pars = args.pars;
    elseif nargin == 2
        pars = g_grind.pars(strcmp(defgroup, g_grind.pargroups));
    else
        pars = varargin(2:end);
    end
    selndx = ones(size(pars));
    if isfield(args,'pargroup')
        pargroup = args.pargroup;
    end
    for i = 1:length(pars)
        ndx = strcmp(g_grind.pars, pars{i});
        if ~any(ndx)
            selndx(i) = 0;
            warning('GRIND:par:unknown','Parameter %s is unknown',pars{i});
        else
            g_grind.pargroups{ndx} = pargroup;
        end
    end
    return;
end
if any(strcmp(args.opts, '-makeshortref'))
    %internal function to make a nice short text of a full reference
    p = args.output;
    if isfield(args,'xtra')
        makelink=args.xtra;
    else
        makelink=false;
    end
    f = strfind(p, ' ');
    if ~isempty(f)
        p = ['%' p(f(1) + 1:end)];
        f1 = strfind(p, '; ');
        if length(f1) > 1
            f2 = strfind(p, ' (');
            f3 = strfind(p,',');
            if ~isempty(f2)&&~isempty(f3)
                p = sprintf('%s et al.%s', p(1:f3(1) - 1), p(f2(1) + 1:end));
            end
        elseif length(f1) == 1
            f3 = strfind(p,',');
            f2 = strfind(p, ' (');
            if length(f3) > 1
                p = sprintf('%s &%s%s', p(1:f3(1) - 1), p(f1(1)+1:f3(2) - 1), p(f2(1)+1:end));
            end
        end
        f1 = strfind(p, 'http://');
        if makelink && ~isempty(f1)
            url = p(f1:end);
            f2 = strfind(p, ')');
            if ~isempty(f2)
                p=sprintf('%%<a href="matlab: web(''%s'',''-browser'')">%s</a>%s',url,p(2:f2),p(f2+1:f1-1));
            else
                p=sprintf('%s\n%%[<a href="matlab: web(''%s'',''-browser'')">full text</a>]',p(1:f1-1),url);
            end
        end
    end
    return;
end
if any(strcmp(args.opts, '-maxparlen'))
    p = 1;
    if isfield(g_grind, 'model')
        for i = 1:length(g_grind.pars)
            if length(g_grind.pars{i}) > p
                p = length(g_grind.pars{i});
            end
        end
    end
    return;
end
if any(strcmp(args.opts, '-d')) %descriptions
    if ~isfield(args,'pars')
        args.pars = g_grind.pars;
    end
    if ischar(args.pars)
        args.pars={args.pars};
    end
    p='';
    p2='';
    pc = cell(1, length(args.pars));
    pc2=pc;
    for i = 1:length(args.pars)
        p='';
        pc2='';
        s = findparcomm(args.pars{i});
        if ~isempty(s)
            f = strfind(s,',');
            if ~isempty(f)
                p = s(1:f(end) - 1);
                p2=  s(f(end) + 1:end);
            else
                p=s;
                p2='';
            end
        end
        pc{i} = p;
        pc2{i} = p2;
    end
    if length(pc) > 1
        p = pc;
        p2 = pc2;
    end
    return;
end
if any(strcmp(args.opts, '-setu'))
    descr = par('-d',args.pars);
    comm = sprintf('%%%s,%s',descr,args.xtra{1});
    p=replaceparcomm(args.pars, comm);
    return
end
if any(strcmp(args.opts, '-set'))
    %you may alternatively give 4 arguments 1)parameter 
    %2)default value, 3) description,4)unit 
    valu=[];
    if iscell(args.xtra)
        if length(args.xtra)==3
            valu=args.xtra{1};
            args.xtra=sprintf('%s,',args.xtra{2:end});
        else
            args.xtra=sprintf('%s,',args.xtra{:});
        end
        args.xtra=args.xtra(1:end-1);
    end
    if isempty(strfind(args.xtra,','))
        comm=sprintf('%%%s,%s',args.xtra,par('-u',args.pars));
    else
        comm = ['%' args.xtra];
    end
    p= replaceparcomm(args.pars, comm, valu);
    if ~isempty(valu)
       assignin('base',args.pars,evalin('base',valu));
    end
    return;
end
if any(strcmp(args.opts, '-u')) %units
    if ~isfield(args,'pars')
        args.pars = g_grind.pars;
    end
    if ischar(args.pars)
        args.pars={args.pars};
    end
    p='';
    pc = cell(1, length(args.pars));
    for i = 1:length(args.pars)
        p = findparcomm(args.pars{i});
        if ~isempty(p)
            f = strfind(p,',');
            if isempty(f)
                p = '';
            else
                p = p(f(end) + 1:end);
            end
        end
        pc{i} = p;
    end
    if length(pc) > 1
        p = pc;
    end
    return;
end
if ~strcmpi(args.output, 'paronly')
    % disp(' ');
    if ~isfield(g_grind, 'model') || isempty(g_grind.model)
        disp('No model defined, use <a href="matlab: model">model</a> or <a href="matlab:vismod">vismod</a> to define a model');
    else
        [~,ininame,iniext]=fileparts(g_grind.inifile);
        if strcmpi(args.output, 'modelonly')
            fprintf('<a href="matlab:par(''-la'')">Model equation(s)</a> [%s%s]\n',ininame,iniext);
        elseif any(strcmp(args.opts, '-size'))
            disp('Sizes of parameters/variables');
            if g_grind.statevars.vector
                dispsiz(g_grind.statevars.vectnames);
            else
                dispsiz(g_grind.statevars.names);
            end
            dispsiz(g_grind.pars);
            return
        else
            fprintf('<a href="matlab:par(''-la'')">Model equation(s)</a> and <a href="matlab:par(''-e'')">parameters</a> [%s%s]\n',ininame,iniext);
        end
        if ~isempty(g_grind.model) && strcmp(g_grind.model{1}, '%external odefile')
            f=~cellfun('isempty',strfind(g_grind.model,'%'));
            fprintf('%s\n',g_grind.model{f});
            type(g_grind.odefile)
        else
            for i = 1:size(g_grind.model, 2)
                if strncmpi(g_grind.model{i}, '%#ref', 5)
                    disp(par('-makeshortref', g_grind.model{i}, 1));
                else
                    disp(char(g_grind.model{i}));
                end
            end
        end
        if isfield(g_grind.solver,'polar')
            makecartesian('-?')
        end
        if ~isempty(g_grind.scheme) && any(strncmp(g_grind.scheme, '%sym=', 4))
            disp('<a href="matlab: vismod">scheme</a>');
        end
    end
end
if isfield(g_grind, 'space')&&~strcmpi(args.output, 'modelonly')
    definespace('-listshort');
end
maxparlen = par('-maxparlen');

if isfield(g_grind, 'pars')&&~strcmpi(args.output, 'modelonly')
    pp = cell(1, numel(g_grind.pars));
    for i = 1:numel(g_grind.pars)
        ppar = evalin('base', g_grind.pars{i});
        if issparse(ppar)
            ppar=full(ppar);
        end
        if (size(ppar, 1) > 1) || (size(ppar, 2) > 1)
            minppar = min(min(ppar));
            if minppar == max(max(ppar))
                if minppar == 0
                    s=sprintf('%s = zeros(%d,%d);',g_grind.pars{i},size(ppar));
                elseif minppar == 1
                    s=sprintf('%s = ones(%d,%d);',g_grind.pars{i},size(ppar));
                else
                    s=sprintf('%s = %0.5g+zeros(%d,%d);',g_grind.pars{i},minppar,size(ppar));
                end
            else
                s=[g_grind.pars{i} ' = [' ];
                siz = size(ppar);
                sh = 0;
                if ~strcmpi(args.output, 'full')
                    if (siz(1) > 10)
                        siz(1) = 10;
                        sh = 1;
                    end
                    if siz(2) > 10
                        siz(2) = 10;
                        sh = 1;
                    end
                end
                
                for j = 1:siz(1)
                    s1=sprintf('%g, ',ppar(j,:));
%                     for k = 1:siz(2)
%                         s = sprintf('%s%g, ',s,ppar(j,k));
%                     end
                    if j < size(ppar, 1)
                        s1(length(s1) - 1) = ';';
                        s = sprintf('%s%s...\n    ', s,s1);
                    else
                        s=[s s1];                       
                    end
                end
                s(length(s) - 1) = ']';
                s(length(s)) = ';';
                if sh
                    s =  sprintf('%s\n%%(first 10 elements shown)', s);
                end
            end
            pp{i} = s;
        elseif ~isempty(ppar)
            if isreal(ppar)
               pp{i}=sprintf(['%-' num2str(maxparlen) 's = %0.5g;'], g_grind.pars{i}, ppar);
            else
               pp{i}=sprintf(['%-' num2str(maxparlen) 's = %0.5g%+0.5gi;'], g_grind.pars{i}, real(ppar),imag(ppar));
            end
        else
            pp{i}=sprintf(['%-' num2str(maxparlen) 's = [];'], g_grind.pars{i});
        end
        [descr,aunit] = par('-d', g_grind.pars{i});
        %  aunit = par('-u', g_grind.pars{i});
        if ~isempty(descr)||~isempty(aunit)
            pp{i}=sprintf(['%-' num2str(maxparlen+25-8) 's %%%s,%s'],pp{i},descr,aunit);
        end
    end
    if g_grind.solver.hasevents
        disp('Discrete events');
        setevent  '-list';
    end
    if nargout == 1
        p = pp;
    else
        if isfield(g_grind, 'pars')
            pg = par('-groups');
            if length(pg)<=1
                fprintf('%s\n', pp{:});
            else
                for j = 1:length(pg)
                    fprintf('%%%s:\n', pg{j});
                    ndx = strcmp(pg{j}, g_grind.pargroups);
                    fprintf('%s\n', pp{ndx});
                end
            end
        end
        if ~isempty(g_grind.externvars)
            defextern();
        end
        i_parcheck(1);
    end
end
function dispsiz(v)
for i = 1:size(v, 2)
    siz = evalin('base', sprintf( 'size(%s)',v{i}));
    fprintf('%s = [%dx%d]\n', v{i}, siz);
end
function [parcomm] = findparcomm(p)
global g_grind;
parcomm = '';
lenp = length(p);
%f=strfind(g_grind.commands,p);
fndx=i_findparcom(p);
% f=regexp(g_grind.commands,sprintf('\\<%s\\>',p));
% fndx=find(~cellfun(@isempty,f));
for i = 1:length(fndx)
    s = strtrim(g_grind.commands{fndx(i)});
    f1 = strfind(s, '%');
    if ~isempty(f1) && (f1(1) < length(s))
        p1=strfind(s, '=');
        if ~isempty(p1)
            par1 = strtrim(s(1:p1(1) - 1));
            if strncmp(p, par1, length(p))
                if (length(s) > lenp) && strcontains(s(lenp + 1:end), '=')
                    parcomm = s(f1 + 1:end);
                    return;
                end
            end
        end
    end
end


function done=replaceparcomm(p, comm, valu)
global g_grind;
fndx=i_findparcom(p);
if ~isempty(fndx)
    s = strtrim(g_grind.commands{fndx});
    p1=strfind(s, '=');
    if ~isempty(p1)
        par1 = strtrim(s(1:p1(1) - 1));
        if strncmp(p, par1, length(p))
            if nargin==2||isempty(valu)
                f1 = strfind(s, '%');
                if isempty(f1)
                    g_grind.commands{fndx} = sprintf('%s    %s', g_grind.commands{fndx}, comm);
                else
                    g_grind.commands{fndx} = [s(1:f1(1) - 1) comm];
                end
            else
                if isnumeric(valu)
                    valu=mat2str(valu);
                end
                g_grind.commands{fndx} = sprintf('%s=%s;    %s',p,valu, comm);
            end
            done=true;
            return;
        end
    end
end
if nargin==2||isempty(valu)
    valu='[]';
end
warning('grind:par','Could not add set parameter attributes as there is no single line assignment, added an extra line')
g_grind.commands=[{sprintf('%s=%s;    %s',p,valu, comm);} g_grind.commands];
done=true;
function [validated_x,errmsg]=i_isparvar(x)
%just the the default check for an equation, but if the equation contains
%short cuts (as we can use in setdata) they are handled
errmsg='';
validated_x='';
if nargin==0
    validated_x='parameter or variable';
    return;
end
if iscell(x)
    for i=1:length(x)
        [~,errmsg]=i_isparvar(x{i});
        if ~isempty(errmsg)
            return
        end
    end
    validated_x=x;
    return;
end
if ischar(x)
    iX=i_getno(x);
    if iX.ispar||iX.isext||iX.isvar
        validated_x = x;
    else
        errmsg = 'argument should be parameter, state variable or external variable';
    end
else
    errmsg = 'argument should be parameter, state variable or external variable';
end

