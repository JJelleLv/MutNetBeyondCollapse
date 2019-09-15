%SETDATA   Enter/edit data for external or state variables
%   Command to enter data for parameter optimizing or for external variables. 
%   You can either enter the data manually or as a file. (see for formats: 
%   <a href="matlab:help loaddata">loaddata</a>)
%
%   Usage:
%   SETDATA - Enter/edit data matrix
%   SETDATA MATRIX - Use MATRIX as data matrix.
%   SETDATA(MATRIX,VARLIST) MATRIX=matrix with data.
%   VARLIST is list with variables (t = time) (if omitted,
%   it is assumed that the first column is time and the
%   other columns are in same order as g_grind.statevars).
%   SETDATA ADATASET - if the <a href="matlab:commands toolboxes">statistics toolbox</a> is installed, you can also
%   use the dataset class to import data (the variable names should be consistent).
%   SETDATA ATABLE - in the newer MATLAB versions you can also
%   use the table class to import data (the variable names should be consistent).
%   SETDATA('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'matrix' [general] - variable  with the data (matrix dataset or table)
%     'varlist' [equation] - list of variables in cell of strings
%   SETDATA('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-clear' - clear the data
%     '-current' - returns the current matrix and varlist in a struct
%     '-list' - list a summary of the data
%
% 
%    See also loaddata, optimpars, defextern, g_data 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setdata')">commands setdata</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function outlist = setdata(varargin)
%(amatrix, varlist)
global g_data g_grind;
deflist = cell(1, length(g_grind.externvars));
k=1;
for i = 1:length(g_grind.externvars)
    if g_grind.externvars{i}.dim1*g_grind.externvars{i}.dim2>1
        d = elemstr(g_grind.externvars{i}.name,g_grind.externvars{i}.dim1,g_grind.externvars{i}.dim2);
        deflist(k:numel(d))=d;
        k=k+numel(d);
    else
        deflist{k} = g_grind.externvars{i}.name;
        k=k+1;
    end
end
if isempty(g_data)||~isfield(g_data,'varlist')
    deflist = [{'t'} i_statevars_names deflist];
else
    deflist = [{'t'} g_data.varlist deflist];
end
fieldnams={'matrix', '', 'variable  with the data (matrix dataset or table)',[];...
    'varlist', 'U1', 'list of variables in cell of strings',deflist}';
args=i_parseargs(fieldnams,'matrix,varlist','-clear,-list,-current',varargin,false,{@i_is_setdata_equation});
if ~isempty(g_data)&&isfield(g_data,'varlist')
    oldvarlist=g_data.varlist;
else
    oldvarlist={};
end
if any(strcmp(args.opts, '-clear'))
    out '-silent' '-remove' '-data';
    out '-silent' '-remove' '-extern';
    g_data.obs = [];
    g_data.pred = [];
    g_data.varlist={};
    g_data.t = [];
    g_data.pars = [];
    g_data.minobs = [];
    g_data.maxobs = [];
    g_grind.checks.lastsettings = [];
    for i = 1:length(g_grind.externvars)
        g_grind.externvars{i}.data =  [];
        g_grind.externvars{i}.nodatawarning=0;
    end
    return;
end
if any(strcmp(args.opts, '-list'))
    outlist=deflist;
    return;
end
if any(strcmp(args.opts, '-current'))
    outlist.varlist = deflist;
    amatrix = [];
    if ~isempty(g_data)
        if isfield(g_data, 'obs')
            amatrix = g_data.obs;
        end
        if isfield(g_data, 't')
            tr = g_data.t;
        else
            tr = [];
        end
    else
        amatrix=nan(1,g_grind.statevars.dim);
        tr=0;
    end
    for i = 1:length(g_grind.externvars)
        g_grind.externvars{i}.nodatawarning=0;
        if ~isempty(g_grind.externvars{i}.data)
            if size(g_grind.externvars{i}.data,2)==1
                [tr, amatrix] = i_concatdata(tr, amatrix, [], g_grind.externvars{i}.data(:, 1));
            else
                [tr, amatrix] = i_concatdata(tr, amatrix, g_grind.externvars{i}.data(:, 1), g_grind.externvars{i}.data(:, 2:end));
            end
        end
    end
    if isempty(tr)
        tr=transpose(1:size(amatrix,1));
    end
    outlist.matrix = [tr amatrix];
    return;
end
i_parcheck;
if isfield(args,'matrix')
    if ischar(args.matrix)
        args.matrix = i_checkstr(args.matrix);
    elseif iscell(args.matrix)
        args.matrix = i_checkstr(char(args.matrix));
    end
    if isa(args.matrix,'dataset') %statistics toolbox
        args.varlist=args.matrix.Properties.VarNames;
        args.matrix=double(args.matrix);
        setdata(args.matrix,args.varlist);
        return;
    elseif isa(args.matrix,'table')
        args.varlist=args.matrix.Properties.VariableNames;
        setdata(args.matrix{:,:},args.varlist);
        return;
    end
end

if isfield(args,'varlist')
    if ischar(args.varlist)
        eoln=sprintf('\n'); %#ok<SPRINTFN>
        v=strrep(args.varlist,sprintf('\t'),eoln);
        v=strrep(v,',',eoln);
        args.varlist=transpose(str2cell(v));
    end
    if size(args.varlist,1)>1
        args.varlist=transpose(args.varlist);
    end
    if g_grind.statevars.vector && (size(args.matrix,2)>length(args.varlist))
        vecvarlist=args.varlist;
        k=1;
        for i=1:length(vecvarlist)
            if isempty(strfind(vecvarlist{i},'('))
                no=i_getno(vecvarlist{i});
                if no.isvar
                    vars=elemstr(vecvarlist{i},g_grind.statevars.dims{no.no}.dim1,g_grind.statevars.dims{no.no}.dim2);
                elseif no.isext && (g_grind.externvars{no.no}.dim1*g_grind.externvars{no.no}.dim2>1)
                    vars=elemstr(vecvarlist{i},g_grind.externvars{no.no}.dim1,g_grind.externvars{no.no}.dim2);
                else
                    vars=vecvarlist(i);
                end
                for l=1:length(vars)
                    args.varlist{k+l-1}=vars{l};
                end
                k=k+length(vars);
            else
                args.varlist{k}=vecvarlist{i};
                k=k+1;
            end
        end
    end
end


TAB = sprintf('\t');
if ~isfield(args,'varlist')
    disp('Assumed order of columns:');
    args.varlist = deflist;
    disp(sprintf('%s ', args.varlist{:}));  %#ok<DSPS>
end
if ~isfield(args,'matrix')
    outvar = setdata('-current');
    [args.matrix, args.varlist] = i_setdatadlg( outvar.matrix, outvar.varlist);
end
if isfield(args,'varlist')
    if ischar(args.varlist)
        l = length(args.varlist) + 1;
        l2 = l - 1;
        while l2 < l
            args.varlist = strrep(args.varlist, [TAB TAB], [TAB '###skip###' TAB]);
            args.varlist = strrep(args.varlist, TAB, ' ');
            args.varlist=strrep(args.varlist,'  ',' ');
            l = l2;
            l2 = length(args.varlist);
        end
        args.varlist=['{''' strrep(strtrim(args.varlist),' ',''' ''') '''}'];
        args.varlist = memo2str(args.varlist, 0);
    end
    args.varlist = i_checkstr(args.varlist);
end
if isempty(args.matrix)
    disp('No data entered or cancelled');
    return;
end
i_savedata(args.varlist, args.matrix);
if ~isempty(g_data)&&~isempty(setdiff(g_data.varlist,oldvarlist))
    out('-silent','-add','-data');
end


function s = memo2str(answer, addmatrixthings)
s1 = i_memo2cell(answer);
%s = '';
f=strfind(s1{1}, '=');
if ~isempty(f)
    s1{1} = s1{1}(f(1) + 1:length(s1{1}));
end
for i = 1:length(s1)
    f = strfind(s1{i}, '...');
    if ~isempty(f)
        s1{i} = s1{i}(1:f(1) - 1);
    end
    if addmatrixthings
        f = strfind(s1{i}, ';');
        if isempty(f)
            s1{i} = [s1{i} ';'];
        end
    end
    %    if addmatrixthings
    %       f = strfind(s1{i}, ';');
    %       if isempty(f)
    %          s = [s ';' s1{i}];
    %       else
    %          s = [s s1{i}];
    %       end
    %    else
    %       s = [s s1{i}];
    %    end
end
s = sprintf('%s', s1{:});
if addmatrixthings && ~strcontains(s, '[')
    s=['[' s ']'];
end

function   An=elemstr(var,dim1,dim2)
An=cell(1,dim1*dim2);
if dim2==1
    for i=1:dim1
        An{i}=sprintf('%s(%d)',var,i);
    end
else
    for i=1:dim1
        for j=1:dim2
            An{(j-1)*dim1+i}=sprintf('%s(%d,%d)',var,i,j);
        end
    end
end
function [validated_x,errmsg]=i_is_setdata_equation(x)
%just the the default check for an equation, but if the equation contains
%short cuts (as we can use in setdata) they are handled
errmsg='';
validated_x='';
if nargin==0
    validated_x='equation';
    return;
end
if iscell(x)
    for i=1:length(x)
        [~,errmsg]=i_is_setdata_equation(x{i});
        if ~isempty(errmsg)
            return
        end
    end
    validated_x=x;
    return;
end
if strncmp(x,'###',3)||isempty(x)
    validated_x=x;
    errmsg='';
else
    [validated_x,errmsg]=i_validat(x,'q',{});
end