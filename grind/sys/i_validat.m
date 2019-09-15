function [validated_x, errormsg] = i_validat(x, cond,validate_funs)
% %
%  c = cell of strings (or string: 'a b c')
%  n = number or equation (or string: '1E-8')
%  i = integer (or string: '100')
%  s = string
%  r = struct
%  l = logical (or 'on' 'off' 'true' 'false' '0' '1')
%  E = empty use only for instance as l#E logical or empty
%  v = exisiting state variable (or list)
%  p = existing parameter (or list)
%  f = function handle
%  U = userfunction
%  F = file name
%  e = enumerated e[a|b|c]
%

global g_grind;
errormsg = '';
if ischar(x)&&strcmp(x,'-make_help_text')
    %use this option to translate the condition text to help text
    if nargin<3
       usertext='';
    else
       usertext=validate_funs;
    end
    if strncmp('e[',cond,2)
        enums=regexp(cond,'[\|\[\]]','split');
        enums=enums(2:end-1);
        for k=1:length(enums)
            if isempty(enums{k})
                enums{k}='''''';
            end
        end
        cond=sprintf(' | %s',enums{:});
        cond=cond(4:end);
        %+ is removed from the text, ++ is kept as "+" 
        cond=regexprep(cond,'(?<=[^+])[+]','');
    else
        reps={'\<c\>','cell';...
        '\<n\>','number';...
        '\<d\>','double';...
       '\<f\>','function_handle';...
       '\<h\>','handle';...
       '\<i\>','integer';...
        '\<s\>','string';...
        '\<r\>','struct';...
       '\<l\>','logical';...
       '\<v\>','state variable';...
        '\<p\>','parameter';...
       '\<q\>','equation';...
        '\<F\>','file name';...
        '\<E\>','empty';...
        '&',' and ';...
        '|',' or '};
        for k=1:length(usertext)
            reps(end+1,:)={sprintf('\\<U%d\\>',k),usertext{k}};
        end
        cond=regexprep(cond,reps(:,1),reps(:,2));
    end
    if isempty(cond)
        cond = 'general';
    end
    validated_x=cond;
    return;
end
if nargin<3
    validate_funs={};
end
if isempty(cond) %do nothing
    validated_x = x;
    return;
end
f=cond=='#';
if any(f)
    f=[0 find(f) length(cond)+1];
    %   conds=regexp(cond,'[#]','split');
    errormsg='';
    for i=2:length(f)
        [validated_x, errormsg1] = i_validat(x, cond(f(i-1)+1:f(i)-1),validate_funs);
        errormsg=sprintf('%s or %s',errormsg,errormsg1);
        if isempty(errormsg1)
            errormsg='';
            return;
        end
    end
    return;
end
vartype= regexp(cond,'\<[rlceshdnivpqfEF]\>|\<U(?<=[0-9]*)','match','once');
if isempty(vartype)
    error('grind:parseargs','unknown type of argument');
end
validated_x = [];
switch vartype
    case 'U' %userfunction for instance i_isid
        no=round(str2double(cond(2:end)));
        [validated_x, errormsg]=validate_funs{no}(x);
    case 'f'
        if iscell(x)
            for i=1:length(x)
                [v1,err]=i_validat(x{i}, cond,validate_funs);
                x{i}=v1;
                if ~isempty(err)
                    errormsg=err;
                end
            end
            validated_x=x;
            if length(validated_x)==1
                validated_x=validated_x{1};
            end
            return;
        end
        if ischar(x)
            if strcontains(x,'@')
                x=evalin('base',x);
            elseif length(regexp(x,'[a-zA-Z_][a-zA-Z_0-9]*','match','once'))==length(x)
                x=str2func(x);
            end
        end
        if ~isa(x,'function_handle')&&~isempty(x)
            errormsg = 'argument should be a function_handle';
        else
            validated_x=x;
        end
    case 'r'
        if ischar(x)
            try
                x=evalin('base',x);
            catch
            end
        end
        if ~isstruct(x)&&~isempty(x)
            errormsg = 'argument should be struct';
        else
            validated_x=x;
        end
    case 'E' %emtpy (only use in combination with other)
        if (ischar(x)&&strcmp(x,'[]'))||isempty(x)
            validated_x=[];
            return;
        else
            errormsg = 'not empty';
        end
    case 'l'
        %         if (ischar(x)&&strcmp(x,'[]'))||isempty(x)
        %             v=[];
        %             return;
        %         end
        if ischar(x)
            try
                x=checklogical(x);
                %                 siz=[];
                %                 if length(cond)>1
                %                     siz=regexp(cond,'l[\[][0-9 :]*[\]]','match','once');
                %                     if ~isempty(siz)
                %                         siz=eval(siz(2:end));
                %                     end
                %                 end
                
                if any(x > 1)||any(x < 0)
                    errormsg = 'value should be logical';
                else
                    validated_x = logical(x);
                    %                     if numel(siz)==1&&numel(v)~=siz
                    %                         errormsg = sprintf('value should be vector of size %d',siz);
                    %                     elseif numel(siz)==2&&any(size(v)~=siz)
                    %                         errormsg = sprintf('value should be matrix of size %dx%s',siz(1),siz(2));
                    %                     end
                end
            catch
                validated_x = [];
                errormsg = 'value should be logical';
            end
        else
            if any(x > 1)||any(x < 0)
                errormsg = 'value should be logical';
            else
                validated_x = logical(x);
            end
        end
        if length(cond) > 1
            l=x;
            if strncmp(cond,'l&',2)
                cond=cond(3:end);
            end
            ok = eval(cond);
            if ~all(ok)
                errormsg = sprintf('value should be logical l with condition: [%s]', cond);
                validated_x=[];
            end
        end
    case 'c'
        if ischar(x)
            if ~isempty(strfind(x,'{')) %#ok<STREMP>
                %from string to cell?
                x = eval(x);
                %         else
                %             %alternatively if a string is entered, we just make one cell of
                %             %it. (\n is used as delimiter)
                %             x=str2cell(x);
            end
        end
        c=x;
        if ~iscell(c)
            errormsg = 'value should be cell';
        else
            if length(cond) > 1
                if strncmp(cond,'c&',2)
                    cond=cond(3:end);
                end
                ok = eval(cond);
                if ~all(ok)
                    errormsg = sprintf('value should be cell c with condition: [%s]', cond);
                end
            else
                validated_x = x;
            end
        end
    case 'e'
        enums=regexp(cond,'[\|\[\]]','split');
        if ischar(x)
            for i=2:length(enums)-1
                s=enums{i};
                f=strfind(s,'+');
                if isempty(f)
                    found=strcmpi(x,s);
                elseif length(f)==2
                    %double ++ is kept as "+" 
                    s=strrep(s,'++','+');
                    found=strcmpi(x,s);
                else
                    %single + is removed from the text (defines how it can
                    %be shortened)
                    found=strncmpi(x,s,f-1);
                    s=s(s~='+'); 
                end
                if found
                    f2=strfind(s,':');
                    if ~isempty(f2)
                        f3=strfind(x,':');
                        if ~isempty(f3)&&f3<length(x)
                            s=[s(1:f2) x(f3+1:end)];
                        end
                    end
                    validated_x=s;
                    return;
                end
            end
        end
        fullcond=sprintf('''%s'' | ',enums{2:end-1});
        fullcond=fullcond(1:end-3);
        fullcond=fullcond(fullcond~='+');
        errormsg = sprintf('value should be %s',fullcond);
    case 'F'
        if ~ischar(x)&&~isempty(x)
            errormsg = 'value should be file name';
        else
            validated_x = x;
        end
    case 's'
        if ~ischar(x)&&~isempty(x)
            errormsg = 'value should be string';
        else
            validated_x = x;
        end
    case 'h'
        if (ischar(x)&&strcmp(x,'[]'))||(isnumeric(x)&&isempty(x))
            validated_x=[];
            return;
        end
        n = checkstr(x);
        if isempty(n)||any(~ishandle(n))
            errormsg = 'value should be a graphical handle';
        else
            validated_x = n;
        end
%             case 'd'
%                 if (ischar(x)&&strcmp(x,'[]'))||(isnumeric(x)&&isempty(x))
%                     v=[];
%                     return;
%                 end
%                 if ischar(x)
%                     d = str2num(x); %#ok<ST2NM>
%                 else
%                     d=x;
%                 end
%                 if isempty(d)
%                     errormsg = 'value should be a double';
%                 elseif isnumeric(d)
%                     if length(cond) > 1
%                         if strncmp(cond,'d&',2)
%                             cond=cond(3:end);
%                         end
%                         ok = eval(cond);
%                         if ~all(ok)
%                             errormsg = sprintf('value should be double d with condition: [%s]', cond);
%                         else
%                             v = d;
%                         end
%                     else
%                         v = d;
%                     end
%                 else
%                     errormsg = sprintf('value should be double d with condition: [%s]', cond);
%                 end
    case {'n','d'}  %d cannot be a parameter
        if (ischar(x)&&strcmp(x,'[]'))||(isnumeric(x)&&isempty(x))
            validated_x=[];
            return;
        end
        if vartype=='d'
            if ischar(x)
                n=str2num(x); %#ok<ST2NM>
            else
                n=x;
            end
        else
            n = checkstr(x);
        end
        if isempty(n)
            errormsg = 'value should be an number';
        elseif isnumeric(n)
            if length(cond) > 1
                if strncmp(cond,'d&',2)
                    d=n;
                    cond=cond(3:end);
                end
                if strncmp(cond,'n&',2)
                    cond=cond(3:end);
                end
                %                 siz=regexp(cond,'n[\[][0-9 :]*[\]]','match','once');
                %                 if ~isempty(siz)
                %                     siz=eval(siz(2:end));
                %                     cond=regexprep(cond,'[\[][0-9 :]*[\]]','');
                %                 end
                if length(cond)>1
                    ok = eval(cond);
                else
                    ok = true;
                end
                if ~all(ok)
                    errormsg = sprintf('value should be number n with condition: [%s]', cond);
                else
                    validated_x = n;
                    %                     if numel(siz)==1&&numel(v)~=siz
                    %                         errormsg = sprintf('value should be vector of size %d',siz);
                    %                     elseif numel(siz)==2&&any(size(v)~=siz)
                    %                         errormsg = sprintf('value should be matrix of size %dx%s',siz(1),siz(2));
                    %                     end
                end
            else
                validated_x = n;
            end
        else
            errormsg = sprintf('value should be number n with condition: [%s]', cond);
        end
    case 'i'
        if (ischar(x)&&strcmp(x,'[]'))||(isnumeric(x)&&isempty(x))
            validated_x=[];
            return;
        end
        if ischar(x)
            i = str2num(x); %#ok<ST2NM>
        else
            i=x;
        end
        if isempty(i)
            errormsg = 'value should be an integer';
        elseif isnumeric(i)&&all(floor(i(:))==i(:))
            if strncmp(cond,'i&',2)
                cond=cond(3:end);
            end
            if length(cond) > 1
                %                 siz=regexp(cond,'i[\[][0-9 :]*[\]]','match','once');
                %                 if ~isempty(siz)
                %                     siz=eval(siz(2:end));
                %                     cond=regexprep(cond,'[\[][0-9 :]*[\]]','');
                %                 end
                if length(cond)>1
                    ok = eval(cond);
                else
                    ok=true;
                end
                if ~all(ok)
                    errormsg = sprintf('value should be integer i with condition: [%s]', cond);
                else
                    validated_x = i;
                    %                     if numel(siz)==1&&numel(v)~=siz
                    %                         errormsg = sprintf('value should be vector of size %d',siz);
                    %                     elseif numel(siz)==2&&any(size(v)~=siz)
                    %                         errormsg = sprintf('value should be matrix of size %dx%s',siz(1),siz(2));
                    %                     end
                end
            else
                validated_x = i;
            end
        else
            errormsg = sprintf('value should be integer i with condition: [%s]', cond);
        end
    case {'v', 'p','q'}
        if vartype=='p'
            if isnumeric(x)&&length(x)==1
                x=g_grind.pars{x};
            end
            msg='value should be a parameter';
        elseif vartype=='v'
            if isnumeric(x)
                x=i_statevars_names(x);
            end
            msg='value should be a state variable';
        else
            msg='value should be a equation';
        end
        if ischar(x)&&strcontains(x,'{')
            x=eval(x);
        end
        if ischar(x)&&(strcontains(x,' ')&&(vartype~='q'))
            x=regexp(x,' ','split');
        end
        if iscell(x)
            errormsg='';
            if (length(cond) > 1)&&~(length(cond)==2&&cond(2)=='+')
                switch vartype
                    case 'p'
                        p=x; %#ok<NASGU>
                    case 'v'
                        v=x; %#ok<NASGU>
                    case 'q'
                        q=x; %#ok<NASGU>
                end
                if strncmp(cond,[vartype '&'],2)
                    cond1=cond(3:end);
                else
                    cond1=cond;
                end
                ok = eval(cond1);
                if ~all(ok)
                    switch vartype
                        case 'p'
                            errormsg = sprintf('value should be parameter p with condition: [%s]', cond1);
                        case 'v'
                            errormsg = sprintf('value should be state variable v with condition: [%s]', cond1);
                        case 'q'
                            errormsg = sprintf('value should be equation q with condition: [%s]', cond1);
                    end
                end
            end
            if isempty(errormsg)
                for i=1:length(x)
                    if ~isempty(x{i})
                        if vartype=='q'
                            [~,errormsg]=checkequation(x{i});
                            if ~isempty(errormsg)
                                break;
                            end
                        else
                            iX=i_getno(x{i});
                            if ~(vartype=='p'&&iX.ispar||vartype=='v'&&iX.isvar)
                                errormsg = msg;
                                break
                            end
                        end
                    end
                end
            end
            if isempty(errormsg)
                validated_x=x;
            end
        elseif ~ischar(x)
            errormsg = msg;
        elseif vartype=='q'
            [validated_x,errormsg]=checkequation(x);
        else
            iX=i_getno(x);
            if (vartype=='p'&&iX.ispar||vartype=='v'&&iX.isvar)
                validated_x = x;
            else
                errormsg = msg;
            end
        end
end
end
function  [validated_x,errormsg]=checkequation(x)
errormsg='';
%shortcut for speed
iX=i_getno(x);
if ~isempty(iX.no)
    validated_x=x;
    return;
end
if strncmpi(x, 'Observed ',9) %obsolete better to use observed('A')
  ivar = strtrim(x(10:end));
  x=sprintf('observed(''%s'')',ivar);
end
validated_x = x;
eq=parsed_equation(x);
try
    eq.structure;
catch err
    errormsg=sprintf('Value should be equation, syntax error %s: %s',x,err.message);
    return;
end
vars=eq.symvar;
for i=1:length(vars)
    iX=i_getno(vars{i});
    if isempty(iX.no)
        errormsg = sprintf('value should be a equation of model variables: "%s" is an unknown variable',vars{i});
        return;
    end
end
end
function x=checklogical(x)
if ischar(x)
    on = true;  %#ok<NASGU>
    off = false;  %#ok<NASGU>
    yes = true;  %#ok<NASGU>
    no = false;  %#ok<NASGU>
    y = true;  %#ok<NASGU>
    n = false;  %#ok<NASGU>
    t = true;  %#ok<NASGU>
    f = false;  %#ok<NASGU>
    x = eval(lower(x));
end
end
function anum=checkstr(astr)
if isempty(astr)
    anum=[];
    return;
end
if ischar(astr)
    anum=str2num(astr);  %#ok<ST2NM>
    if isempty(anum)||~(isnumeric(anum)||islogical(anum))
        try
            anum=evalin('base',astr);
        catch
            anum=[];
        end
    end
else
    anum=astr;
end
end

function res=ndim %#ok<DEFNU> %used in conditions
global g_grind
res=g_grind.statevars.dim;
end
