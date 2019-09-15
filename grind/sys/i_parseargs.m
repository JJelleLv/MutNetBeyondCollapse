function [res,defaults] = i_parseargs(fieldnams, deffields, opts, args, unknownopts,validate_funs)
%validation : s[tring] l[ogical] n[umeric] i[nteger] c[ell] 'n>0' 'i>0&&i<100' 'n<0 v=statevariable p=parameter, r=struct'
%example:
%res=i_parseargs('ndays,steps(i>0||i<5),type(s),test(l),disturbance(n)','ndays,steps,disturbance','-a,-s,-p1',varargin);
%res=i_parseargs('ndays,steps(i>0||i<5),type(s),test(l),disturbance(n)','if (nargs==2),deffields=''steps,disturbance'',end;if(nargs==1&&iscell(arg{1}),deffields=''disturbance'';end','-a,-s,-p1',varargin);
%
%some nested functions that my be used in the deffields definition (for
%convenience
    function n=nargs
        n = length(args);
    end
    function ok=hasoption(anopt)
        %can have a list of options : '-opt1,opt2,opt3'
        ok=false;
        if ~isempty(strfind(anopt,',')) %#ok<STREMP>
            opts1=regexp(anopt,',','split');
            for m=1:length(opts1)
                ok=hasoption(opts1{m});
                if ok
                    return;
                end
            end
            return;
        end
        for g=1:length(args)
            if ischar(args{g})&&strncmp(args{g},anopt,length(anopt))
                ok=true;
                return;
            end
        end
    end
    function ok=argtype(i,aformat)  %#ok<DEFNU>
        %ignore options
        while i<=length(args)&&any(strncmp(args{i},'-',1))
            i=i+1;
        end
        if i<=length(args)
            [~, errormsg] = i_validat(args{i}, aformat,validate_funs);
            ok=isempty(errormsg);
        else
            ok=false;
        end
    end
if nargin < 6
    validate_funs={};
end
if ~iscell(validate_funs)
    validate_funs={validate_funs};
end
if nargin < 5
    unknownopts = false;
end
if nargout==2
    if size(fieldnams,1)<4
        defaults=[];
    else
        s=fieldnams([1 4],:);
        defaults=struct(s{:});
    end
end
res.opts = {};
if isempty(args)||numel(args)==1&&isstruct(args{1})&&isfield(args{1},'opts')
    %all grind commands can also be accessed with the resulting structure (including
    %the field "opts"), then there is no parsing/checks for maximum speed
    %use this option for internal calls.
    if ~isempty(args)
        res=args{1};
    end
    return;
end
if ischar(fieldnams)
    fieldnams=regexp(fieldnams,',','split');
end
if iscell(fieldnams)&&(size(fieldnams,1)~=3&&size(fieldnams,1)~=4)
    fields=regexp(fieldnams,'[a-zA-Z_0-9]*','once','match');
    validation=regexp(fieldnams,'(?<=[\(]).*(?=[\)])','once','match');
    explan=cell(size(fields));
    fieldnams=[fields;validation;explan];
else
    fields=fieldnams(1,:);
    validation=fieldnams(2,:);
    explan=fieldnams(3,:);
end

deffieldcondition='';
if ischar(deffields)
    if ~isempty(strfind(deffields, 'deffields=')) %#ok<STREMP>
        deffieldcondition=deffields;
        eval(deffields);
    end
    if ischar(deffields)
        deffields=regexp(deffields,',','split');
    end
end
if ischar(opts)
    opts=regexp(opts,',','split');
end
% if length(args)==1&&ischar(args{1})&&strcmp(args{1},'??1')
%     %this option is used by help_opt.m to find all options
%     s=sprintf('''%s'',',fieldnams{:});
%     s1=sprintf('res.fields = {%s};\n',s(1:end-1));
%     if ~isempty(deffieldcondition)
%         s=strrep(deffieldcondition,'''','''''');
%         s1=sprintf('%sres.deffields = ''%s'';\n',s1,s);
%     else
%         s=sprintf('%s,',deffields{:});
%         s1=sprintf('%sres.deffields = ''%s'';\n',s1,s(1:end-1));
%     end
%     s=sprintf('''%s'',',opts{:});
%     s1=sprintf('%sres.opts = {%s};\n',s1,s(1:end-1));
%     if unknownopts
%         s1=sprintf('%sres.unknownopts = true;\n',s1);
%     else
%         s1=sprintf('%sres.unknownopts = false;\n',s1);
%     end
%     disp(s1)
%     error('grind:parseargs','Argument "???" shows command line options and stops execution');
% end
splitopts=regexp(opts,'[|]','split');

%validation=regexp(fieldnams,'(?<=[\(])[\<\>a-zE0-9\[\]\(\)\-+.\ |:&~=#]*(?=[\)])','once','match');
used = false(size(args));
if length(args)>=1&&ischar(args{1})
    if any(strcmp(args{1},{'??','???'}))
        gethelp(fields, splitopts, deffields,deffieldcondition, validation, explan, unknownopts,strcmp(args{1},'???'),validate_funs);
        stoprun()
    elseif strncmp(args{1},'??dlg',5)
        aa.opts=args(2:end);
        if ~exist('i_use','file')
            addpath([grindpath filesep 'sys2']);
        end
        if strcmp(args{1},'??dlgnoopt')
            opts={};
        end
        res=i_parseargdlg(getcallingfunction,fieldnams,validate_funs,opts,aa);
        if isempty(res)
            stoprun()
        end
        return;
    elseif strcmp(args{1},'??makecell')
        res=str2cell(evalc(sprintf('help %s',getcallingfunction())));
        fprintf('fieldnams={')
        for i=1:length(fields)
            aa=regexp(res,sprintf(' ''%s''.*',fields{i}),'match','once');
            ndx=find(~cellfun('isempty',aa));
            if isempty(ndx)
                helptxt=' ';
            else
                helptxt=aa{ndx(end)};
                helptxt=strtrim(regexp(helptxt,'(?<=[-]).*','match','once'));
                helptxt=strrep(helptxt,'''','''''');
            end
            if i>1
                fprintf('   ');
            end
            if i==length(fields)
                fprintf('''%s'', ''%s'', ''%s''}'';\n',fields{i},validation{i},helptxt);
            else
                fprintf('''%s'', ''%s'', ''%s'';...\n',fields{i},validation{i},helptxt);
            end
        end
        stoprun();
    end
    %error('grind:parseargs','Argument "??" shows command line options and stops execution');
end

structndx=find(cellfun(@isstruct,args));
for k = 1:length(structndx)
    i=structndx(k);
    if isstruct(args{i})
        sfields = lower(fieldnames(args{i}));
        f=strcmp(sfields,'opts'); %remove the field opts (this makes it easier to use parsed args, and a field cannot be otps)
        values=struct2cell(args{i});
        sfields=sfields(~f);
        values=values(~f);
        sdiff = setdiff(sfields, lower(fields));
        if isempty(sdiff)
            used(i) = true;
            arg1 = transpose([sfields, values]);
            arg1=transpose(arg1(:));
            args = [args, arg1];
            used=[used false(1,length(arg1))];
        elseif length(sdiff) < length(sfields)
            %gethelp(fields, splitopts, deffields, deffieldcondition,validation,explan, unknownopts,false,validate_funs)
            disp(sdiff);
            parseerror('grind:parseargs','structure has unknown fields');
        end
    end
end
for i = length(args):-2:2
    if ~used(i)&&ischar(args{i-1})
        f = strcmpi(args{i-1}, fields);
        if any(f)
            [value, errmsg] = i_validat(args{i}, validation{f},validate_funs);
            if ~isempty(errmsg)
                if i>length(deffields)
                    %gethelp(fields, splitopts, deffields,deffieldcondition, validation,explan, unknownopts,false,validate_funs)
                    parseerror('grind:parseargs','Error in parameter/value pair (%s): %s',args{i-1},errmsg);
                else
                    break;
                end
            end
            res.(fields{f})=value;
            used(i-1:i) = true;
        else
            break;
            %             if i >= length(args)
            %                 error('grind:parseargs','Error in parameter/value pair "%s"',args{i});
            %             else
            %                 gethelp(fields, opts, deffields,deffieldcondition, validation, unknownopts,false,validate_funs)
            %                 error('grind:parseargs','Unknown parameter "%s"',args{i-1});
            %             end
        end
    else
        break;
    end
end

usedfields = fieldnames(res);
if any(used)
    args = args(~used);
    used = false(size(args));
end
%find options with minus
found = false;
%isopt=false(size(args));
for i = 1:length(args)
    if ischar(args{i})
        optndx = strfind(args{i}, '-');
        if ~isempty(optndx)&&optndx(1)==1%&&length(optndx)==1
            %            isopt(i)=true;
            %[~, errormsg] = i_validat(args{i}, 'n',validate_funs);
            % isneg=isempty(errormsg); %is also a valid negative number?
            for j = 1:length(splitopts)
                opt = splitopts{j};
                found = false;
                for k = 1:length(opt)
                    isregexp=~strncmp(opt{k},'-',1);
                    if (~isregexp&&strncmp(opt{k}, args{i}, length(opt{k})))||(isregexp&&~isempty(regexp(args{i},opt{k},'match','once')))
                        f=strfind(args{i}, '=');
                        if ~isempty(f)
                            res.opts = [res.opts; [opt{1} args{i}(f:end)]];
                        elseif isregexp
                            res.opts = [res.opts; regexp(args{i},opt{k},'match','once')];
                        else
                            res.opts = [res.opts; opt{1}];
                        end
                        used(i) = true;
                        found = true;
                        break;
                    end
                end
                if found
                    break;
                end
            end
            if ~found&&unknownopts%&&isempty(regexp(args{i},'[\(\)+*^/]' ,'once'))&&~isneg %still difficult point if there are equations with minus
                res.opts = [res.opts; args{i}];
                used(i) = true;
            end
        end
    end
end
if any(used)
    args = args(~used);
    %   isopt = isopt(~used);
end
i1=1;
inplus=false;
for i = 1:length(args)
    deffield=regexp(deffields{i1},'[a-z_0-9]*','once','match');
    k = find(strcmp(deffield, fields));
    if isempty(k)
        error('grind:parseargs','Error in field definition, bug in GRIND')
    else
        vali = validation{k};
    end
    [aval, errmsg] = i_validat(args{i}, vali,validate_funs);
    if ~isempty(errmsg)&&inplus&&i<=length(args)
        inplus=false;
        i1=i1+1;
        if i1<length(deffields)
        deffield=regexp(deffields{i1},'[a-z_0-9]*','once','match');
        k = strcmp(deffield, fields);
        vali = validation{k};
        [aval, errmsg] = i_validat(args{i}, vali,validate_funs);
        end
    end
    if ~isempty(errmsg)
        parseerror('grind:parseargs','Error in argument %d (%s): %s',i,deffield,errmsg);
    else
        %check if the field is also used with the names, if confict give error
        if any(strcmp(deffield, usedfields))
            val2 = res.(deffield);
            if (ischar(aval)&&ischar(val2)&&~strcmp(aval, val2))||(isnumeric(aval)&&isnumeric(val2)&&(aval~=val2))
                parseerror('grind:parseargs','Conflict in the input argument "%s", twice defined',deffield);
            end
        else
            if ~isempty(strfind(deffields{i1},'(+)'))
                inplus=true;
                if ~iscell(aval)
                    aval={aval};
                end
                if ~isfield(res,deffield)
                    res.(deffield)=aval;
                else
                    res.(deffield)(end+1) = aval;
                end
                i1=i1-1;
            else
                %no errors, no plus
                res.(deffield) = aval;
            end
        end
        i1=i1+1;
    end
end

end

function stoprun(msg)
if nargin==0
    msg='';
end
ms.message=msg;
ms.stack = dbstack;
ms.stack(1:end) = [];
error(ms);
end
function res=getcallingfunction()
ss=dbstack();
ss={ss.name};
f=find(strcmp(ss,'i_parseargs'));
if f==length(ss)
    res='';
else
    res=ss{f(1)+1};
    if strcmp(res,'grind_matcont.set')
        res='conteq';
    end
    f1=strfind(res,'.');
    if ~isempty(f1)
        res=res(f1(1)+1:end);
    end
    if strncmp(res,'i_',2)==1
        res=ss{f(1)+2};
    end
end
end
function gethelp(fields, opts, deffields,deffieldcondition, validation, explan, unknownopts, addhelptext,validate_funs)
%find calling function

usertext=cell(size(validate_funs));
for k=1:length(validate_funs)
    usertext{k}=feval(validate_funs{k});
end
[~,ndx]=sort(lower(fields));
fields=fields(ndx);
validation=validation(ndx);
explan=explan(ndx);
opts1=cell(size(opts));
for i=1:length(opts)
    oo=opts{i};
    ndx=~cellfun('isempty',oo);
    oo=oo(ndx);
    if ~isempty(oo)
        for j=1:length(oo)
            if strcmp(oo{j},'(-[1-9]+[0-9]*)')
                oo{j}='-1 or -2 or -value';
            end
        end
        oopts1=strtrim(sprintf(' or ''%s''',oo{:}));
        opts1{i}=oopts1(4:end);
    else
        opts1{i}='';
    end
end
[~,ndx]=sort(lower(opts1));
opts=opts1(ndx);
ndx=cellfun('isempty',opts);
opts=opts(~ndx);
callingfunction=getcallingfunction();
% ss=dbstack();
% ss={ss.name};
% f=find(strcmp(ss,'i_parseargs'));
% if f==length(ss)
%     callingfunction='';
% else
%     callingfunction=ss{f(1)+1};
%     if strcmp(callingfunction,'grind_matcont.set')
%         callingfunction='conteq';
%     end
%     f1=strfind(callingfunction,'.');
%     if ~isempty(f1)
%         callingfunction=callingfunction(f1(1)+1:end);
%     end
%     if strncmp(callingfunction,'i_',2)==1
%         callingfunction=ss{f(1)+2};
%     end
% end

nofields=true;
if addhelptext
    fprintf('    %s(''argname'',argvalue,...) - Valid argument <a href="func_args.htm">name-value pairs</a> [with type]:<br>\n',upper(callingfunction));
else
    fprintf('<strong>%s</strong> Definition of <a href="matlab:commands func_args">flexible arguments</a>:\n',callingfunction)
end
if ~(isempty(deffields)||numel(deffields)==1&&isempty(deffields{1}))
    if isempty(deffieldcondition)
        disp('   Order of normal arguments:')
        for i = 1:length(deffields)
            fprintf('      %d: ''%s''\n',i,deffields{i});
        end
    else
        %     s1=sprintf(strrep(strrep(deffieldcondition,';',';\n      '),',deffields','\n        deffields'));
        ndx=strfind(deffieldcondition,'deffields=''');
        s1=deffieldcondition;
        for i=length(ndx):-1:1
            s2=s1(ndx(i)+11:end);
            s1=s1(1:ndx(i)+10);
            f1=strfind(s2,'''');
            s3=s2(f1(1)+1:end);
            s2=s2(1:f1(1)-1);
            vars=regexp(s2,',','split');
            for j=1:length(vars)
                vars{j}=sprintf('        %d: %s\n',j,vars{j});
            end
            vars=[{s1},vars,{s3}];
            s1=sprintf('%s',vars{:});
        end
        %        s1= sprintf(strrep(strrep(s1,';','\n      '),',deffields=''','\n'));
        s1= sprintf(strrep(strrep(s1,';','      '),',deffields=''','\n'));
        fprintf('   Order of normal arguments (context dependent):\n      %s\n',s1);
    end
    
    nofields=false;
end
if ~(isempty(fields)||numel(fields)==1&&isempty(fields{1}))
    if ~addhelptext
        %disp('   Valid argument name-value pairs:')
        fprintf('   <strong>%s</strong>(''argname'',argvalue,...) - Valid argument name-value pairs [with type]:\n',callingfunction);
    end
    for i = 1:length(fields)
        valids=regexp(validation{i},'#','split');
        for j=1:length(valids)
            valids{j}=i_validat('-make_help_text',valids{j},usertext);
        end
        vali=sprintf(' or %s',valids{:});
        vali=vali(5:end);
        if isempty(explan{i})
            explan1='';
        else
            explan1=explan{i};
        end
        if addhelptext
            vali=strrep(strrep(vali,'<','&lt;'),'>','&gt;');
            fprintf('   &nbsp;&nbsp; ''%s'' [%s] - %s<br>\n',fields{i},lower(vali),explan1);
        else
            fprintf('      ''%s'' [%s] - %s\n',fields{i},lower(vali),explan1);
        end
    end
    nofields=false;
end

if addhelptext
    fprintf('    %s(''-opt1'',''-opt2'',...) - Valid command line <a href="func_args.htm">options</a>:<br>\n',upper(callingfunction));
    for i = 1:length(opts)
        fprintf('   &nbsp;&nbsp; %s - <br>\n',opts{i});
    end
elseif ~(isempty(opts)||numel(opts)==1&&numel(opts{1})==1&&isempty(opts{1}{1}))
    if unknownopts
        fprintf('   <strong>%s</strong>(''-opt1'',''-opt2'',...) - Valid command line options (not complete):\n',callingfunction);
    else
        fprintf('   <strong>%s</strong>(''-opt1'',''-opt2'',...) - Valid command line options:\n',callingfunction);
    end
    for i = 1:length(opts)
        fprintf('      %s - \n',opts{i});
    end
    nofields=false;
end


if nofields
    if unknownopts
        disp('No arguments allowed (but there can be valid options)')
    else
        disp('No arguments or options allowed')
    end
end
if ~addhelptext
    fprintf('   More help: <a href="matlab:commands %s">commands %s</a> \n',...
        callingfunction,callingfunction)
end

end


function parseerror(id,msg,varargin)
astack=dbstack('-completenames');
fun=getcallingfunction();
fprintf('Error in arguments of grind command "%s", (see <a href="matlab:%s(''??'')">valid arguments</a>)\n', fun,fun)
        %gethelp(fields, splitopts, deffields, deffieldcondition,validation, explan,unknownopts,false,validate_funs);
error(struct('msgid',id,'message',sprintf(msg,varargin{:}),'stack',astack(2:end)));
end
