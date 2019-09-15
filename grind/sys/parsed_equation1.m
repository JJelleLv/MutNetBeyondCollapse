classdef parsed_equation < handle
    % Class parsed_equation implements an equation (MATLAB style) that can be parsed.
    % It uses the matlab <a href="matlab: help precedence">precedence</a> rules by default.
    %
    % Example:
    % obj=parsed_equation('x=a*exp(b/4)-4.*(3^2)+5') returns
    %   x=a*exp(x/4)-4.*(3^2)+5
    %
    % obj.latex
    %   $${x=a\; e^{\frac{b}{4}}-4 \cdot  \left(3^{2} \right)+5}$$
    %
    % obj.minbrackets
    %   a*exp(b/4)-4.*3^2+5
    %
    % obg.maxbrackets
    %   (((a*(exp(b/4)))-(4.*(3^2)))+5)
    %
    % obj.fs
    %   {'x' '=' 'a' '*' 'exp' '(' 'b' '/' '4' ')' '-' '4' '.*' '(' '3' '^' '2' ')' '+' '5'}
    %
    % obj.expanded
    %  'tmp001 = b/4;
    %  tmp002 = exp(tmp001);
    %  tmp004 = 3^2;
    %  tmp003 = a*tmp002;
    %  tmp005 = 4.*tmp004;
    %  tmp006 = tmp003-tmp005;
    %  x = tmp006+5;'
    %
    
    properties
        equation = ''
    end
    
    properties (SetAccess = public) %in fact dependent properties, but stored to avoid recalculation
        fs = []
        types = []
        structure =  []
    end
    
    properties (Constant)
        vr = struct('empty', - 99 ...
            ,'comment',  - 1 ...
            ,'number', 1 ...
            ,'field', 2 ...
            ,'brack1', 3 ...
            ,'brack2', 4 ...
            ,'accent', 5 ...
            ,'comma', 6 ...
            ,'colon', 7 ...
            ,'diff', 8 ...
            ,'t', 11 ...
            ,'space', 12 ...
            ,'semicolon', 13 ...
            ,'sqbrack1', 14 ...
            ,'sqbrack2', 15 ...
            ,'if', 16 ...
            ,'else', 17 ...
            ,'oper', 18 ...
            ,'unary_oper_pre', 19 ...
            ,'unary_oper_post', 20 ...
            ,'defextern', 21 ...
            ,'definepars', 22 ...
            ,'lag', 23 ...
            ,'dwiener', 24 ...
            ,'rednoise', 25 ...
            ,'block1', 26 ...
            ,'block2', 27 ...
            ,'boxcartrain', 28 ...
            ,'setevent', 29 ...
            ,'implicitdisperse', 30 ...
            ,'kw', 31 ...
            ,'assign', 32 ...
            ,'defpermanent', 33 ...
            ,'function', 34 ...
            ,'string', 35 ...
            ,'definespace', 36 ...
            ,'tdiff', 37 ...
            ,'statevar_dim1', 38 ...
            ,'statevar_dim2', 39 ...
            ,'djump',40 ...
            ,'externlag', 41 ...
            ,'funcname', 100 ...
            ,'locfun', 101 ...
            ,'stochfun', 102 ...
            ,'constant', 103 ...
            ,'function_handle', 104 ...
            ,'variable', 200 ...
            ,'statevar', 201 ...
            ,'auxil', 202 ...
            ,'permanent', 203 ...
            ,'extern', 204 ...
            ,'parameter', 205 ...
            ,'statevar_p',206 ...
            ,'implicitvars',207....
            ,'mindiff', 50000);
    end
    
    % Class methods
    methods
        function obj = parsed_equation(varargin)
            if nargin > 0
                if nargin == 1
                    s = varargin{1};
                else
                    s = varargin;
                end
                % Construct a parsed_equation object
                if isa(s, 'parsed_equation')
                    obj.equation = s.equation;
                    obj.fs = s.fs;
                    obj.types = s.types;
                    obj.structure = s.structure;
                elseif ischar(s)
                    obj.equation = s;
                elseif isstruct(varargin{1})
                    if nargin < 2
                        opersinfo = parsed_equation.operinfo;
                    else
                        opersinfo = varargin{2};
                    end
                    fs1 = makeequation({}, varargin{1}, length(varargin{1}), opersinfo);
                    obj.equation = sprintf('%s', fs1{:});
                    obj.fs = fs1;
                    obj.structure = varargin{1};
                elseif iscell(s)&&~isempty(s)
                    obj(1,length(s)) = parsed_equation('');
                    fss=parseeq(s);
                    for i = 1:length(s)
                        obj(i).equation=sprintf('%s',fss{i}{:});
                        obj(i).types=[];
                        obj(i).fs = fss{i};
                    end
                end
            end
        end
        
        function value = get.structure(obj)
            %structure is only calculated if needed
            if isempty(obj.structure)
                if ~isempty(strfind(obj.types,[1 1]))||~isempty(strfind(obj.types,[200 1]))||~isempty(strfind(obj.types,[200 200]))||~isempty(strfind(obj.types,[1 200]))
                    error('parsed_equation:twonumbers','Parsing error in equation "%s" - pair of numbers/variables without an operator',obj.equation);
                end
                if ~isempty(strfind(obj.types,[18 18]))
                    error('parsed_equation:twoopers','Parsing error in equation "%s" - two successive operators without variables',obj.equation);
                end
                sumbrackets=sum(strcmp(obj.fs,'('))-sum(strcmp(obj.fs,')'));
                if sumbrackets < 0
                    obj.fs = [];
                    error('parsed_equation:bracks','Parsing error in equation "%s" - unbalanced brackets, missing "("',obj.equation);
                elseif sumbrackets > 0
                    obj.fs = [];
                    error('parsed_equation:bracks','Parsing error in equation "%s" - unbalanced brackets, missing ")"',obj.equation);
                end
                obj.structure = expandeq(obj.addbrackets(true), obj.fs, obj.types);
                resolvedvars={};
                allvars=obj.symvar;
                usedvars=false(size(allvars));
                for i=1:length(obj.structure)
                    tmpvar=obj.structure(i).args(strncmp(obj.structure(i).args,'#',1));
                    if any(~ismember(tmpvar,resolvedvars))
                        error('parsed_equation:structure','Parsing error in equation "%s" - error in operator %s, character %d',obj.equation,obj.structure(i).oper,obj.structure(i).charnr);
                    end
                    resolvedvars{end+1}=obj.structure(i).leftvar; %#ok<AGROW>
                    usedvars=usedvars | ismember(allvars,obj.structure(i).args);
                end
                usedvars=usedvars | ismember(allvars,obj.structure(end).leftvar);
                if any(~usedvars)
                    error('parsed_equation:varnotused','Parsing error in equation "%s" - syntax error',obj.equation);
                end         
            end
            value = obj.structure;
        end
        
        function value = get.fs(obj)
            %fs is only calculated if needed
            if isempty(obj.fs)
                [obj.fs] = parseeq(obj.equation);
            end
            value = obj.fs;
        end
        
        
        function value = get.types(obj)
            %types is only calculated if needed
            if isempty(obj.types)
                obj.types = updatetypes(obj.fs,obj.vr);
            end
            value = obj.types;
        end
        function set.equation(obj, s)
            if ~strcmp(s, obj.equation)
                obj.equation = s;
                obj.fs = [];  %#ok<MCSUP>
                obj.types = [];  %#ok<MCSUP>
                obj.structure = [];  %#ok<MCSUP>
            end
        end
        function obj1=subs(obj,old,new)
            %substitute old with new (both can be variable or equation)
            %Note: if old is an equation it cannot be a string but should be a
            %parsed_equation
            %new can both be parsed_equation or string (makes no
            %difference)
            if nargin==2
                if ischar(old)
                    old=parsed_equation(old);
                end
                f=find(strcmp(old.fs,'='),1);
                new=old.fs(f+1:end);
                old=old.fs(1:f-1);
            end
            if ischar(new)
                new={new};
            elseif isa(new,'parsed_equation')
                new=new.fs;
            end
            if length(new)>1&&strcmp(new{end},';')
                new=new(1:end-1);
            end
            if isa(old,'parsed_equation')
                old=old.fs;
                %now we can search for the sub-equation in obj
            end
            if length(obj) > 1
                obj1=obj;
                for i = 1:length(obj)
                    obj1(i) = subs(obj(i), old,new);
                end
                return;
            end
            fsold=obj.fs;
            feq=find(strcmp(fsold,'='),1);
            if isempty(feq)
                feq=0;
            end
            if ischar(old)
                f=find(strcmp(fsold,old));
                if ~isempty(f)
                    fs1=fsold;
                    for i=length(f):-1:1
                        if f(i)>feq
                            fs1=[fs1(1:f(i)-1) {'('} new {')'} fs1(f(i)+1:end)];
                        end
                    end
                    obj1=parsed_equation(sprintf('%s', fs1{:}));
                end
            elseif ~isempty(old)
                f=find(strcmp(fsold,old{1}));
                if ~isempty(f)
                    fs1=fsold;
                    for i=length(f):-1:1
                        if (f(i)>feq)&&f(i)<=length(fsold)-length(old)+1
                            found=true;
                            for j=1:length(old)
                                if ~strcmp(fsold{f(i)+j-1},old{j})
                                    found=false;
                                end
                            end
                            if found
                                fs1=[fs1(1:f(i)-1) {'('} new {')'} fs1(f(i)+length(old):end)];
                            end
                        end
                    end
                    obj1=parsed_equation(sprintf('%s', fs1{:}));
                else
                    obj1=parsed_equation(obj.equation);
                end
            else
                obj1=parsed_equation(obj.equation);
            end
            
        end
        function [obj,result] = mupad_syntax(obj)
            %all small words <4 that are reserved in mupad (the entire list is
            %1500 words long) (esp. E I O D can be a problem)
            mupad_words={'new','N_','GC','lhs','All','Re','Bin','SVD','_or','psi','Rem','arg','int','R_',...
                'E','GT','PDF','rhs','ToY','I','i','New','_in','C_','Ci','is','gcd','RGB','Chi','Dom','Z_','Int',...
                'CDF','has','adt','O','QRD','lcm','op','id','_if','ToX','Mem','Quo','Li','Im',...
                'Ax','hfa','Q_','det','rec','D','Not','val','To','End','Shi','Cat','NC','PF',...
                'fp','Raw','RK4','Low','Pop','Any','use','map','fix','Min','Up','Max','ode','Ei','zip','Top','Ssi',...
                'Gap','ToZ','GL','Lex','beta','gamma','theta','zeta'};
            result.allvars={};
            for i=1:length(obj)
                result.allvars = [result.allvars obj(i).fs];
            end
            result.allvars=unique(result.allvars(~strncmp(result.allvars,'%',1)));
            result.renamewords = intersect(result.allvars, mupad_words);
            result.rename2words = cell(size(result.renamewords));
            if ~isempty(result.renamewords)
                for i = 1:length(result.renamewords)
                    k = 1;
                    while any(strcmp(result.allvars, sprintf('%s%d', result.renamewords{i}, k)))
                        k = k + 1;
                    end
                    result.rename2words{i} = sprintf('%s%d', result.renamewords{i}, k);
                end
            end
            stringfuncs=intersect(result.allvars,{'val','solver'});
            m=1;
            if ~isempty(stringfuncs)
                for i=1:length(obj)
                    for j=1:length(stringfuncs)
                        f=find(strcmp(obj(i).fs,stringfuncs{j}));
                        for k=1:length(f)
                            [~,ndx]=getfunctionpars(obj(i),f(k));
                            if size(ndx,1)==1
                                strng=strncmp(obj(i).fs(ndx(1,1):ndx(1,2)),'''',1);
                                if ~any(strng)
                                    strng=strncmp(obj(i).fs(ndx(1,1):ndx(1,2)),'"',1);
                                end
                                fs1=ndx(1,1)+find(strng)-1;
                                result.renamewords(end+1)=obj(i).fs(fs1);
                                result.rename2words(end+1)={sprintf('g_string%d',m)};
                                m=m+1;
                            end
                        end
                    end
                end
            end
            hasiif=any(strcmp(result.allvars,'iif'));
            if hasiif
                for i=1:length(obj)
                    f=find(strcmp(obj(i).fs,'iif'));
                    if ~isempty(f)
                        for j=length(f):-1:1
                            [~,ndx]=getfunctionpars(obj(i),f(j));
                            if size(ndx,1)==3
                                obj(i).fs{f(j)}='piecewise';
                                obj(i).fs{ndx(1,1)-1}=[ obj(i).fs{ndx(1,1)-1} '['];
                                notcond= ['_not(',sprintf('%s',obj(i).fs{ndx(1,1):ndx(1,2)}),')'];
                                %notcond='Otherwise';
                                obj(i).fs{ndx(2,2)+1}=['],[' notcond obj(i).fs{ndx(2,2)+1}];
                                obj(i).fs{ndx(3,2)+1}=[']' obj(i).fs{ndx(3,2)+1} ];
                                obj(i)=parsed_equation(sprintf('%s',obj(i).fs{:}));
                            end
                        end
                    end
                end
            end
            hasminmax=any(strcmp(result.allvars,'min'))||any(strcmp(result.allvars,'max'));
            if hasminmax
                for i=1:length(obj)
                    fminmax=find(strcmp(obj(i).fs,'min')|strcmp(obj(i).fs,'max'));
                    if ~isempty(fminmax)
                        for j=length(fminmax):-1:1
                            [~,ndx]=getfunctionpars(obj(i),fminmax(j));
                            if size(ndx,1)==2
                                fs1=obj(i).fs;
                                if strcmp(obj(i).fs{fminmax(j)},'min')
                                    arg1=sprintf('%s',obj(i).fs{ndx(1,1):ndx(1,2)});
                                    arg2=sprintf('%s',obj(i).fs{ndx(2,1):ndx(2,2)});
                                else
                                    arg2=sprintf('%s',obj(i).fs{ndx(1,1):ndx(1,2)});
                                    arg1=sprintf('%s',obj(i).fs{ndx(2,1):ndx(2,2)});
                                end
                                fs1{fminmax(j)}=sprintf('piecewise([%s<%s',arg1,arg2);
                                fs1{fminmax(j)+1}=',';
                                fs1{ndx(1,2)+1}=['],[Otherwise' obj(i).fs{ndx(1,2)+1}];
                                fs1{ndx(2,2)+1}=[']' obj(i).fs{ndx(2,2)+1} ];
                                obj(i)=parsed_equation(sprintf('%s',fs1{:}));
                            end
                        end
                    end
                end
            end
            if any(strcmp(result.allvars,'log10'))
                for i=1:length(obj)
                    f=find(strcmp(obj(i).fs,'log10'));
                    if ~isempty(f)
                        for j=length(f):-1:1
                            obj(i).fs{f(j)}='(1/log(10))*log';
                        end
                        obj(i)=parsed_equation(sprintf('%s',obj(i).fs{:}));
                    end
                end
            end
            if any(strcmp(result.allvars,'log2'))
                for i=1:length(obj)
                    f=find(strcmp(obj(i).fs,'log2'));
                    if ~isempty(f)
                        for j=length(f):-1:1
                            obj(i).fs{f(j)}='(1/log(2))*log';
                        end
                        obj(i)=parsed_equation(sprintf('%s',obj(i).fs{:}));
                    end
                end
            end
            %these are also changed back
            renameopers={'diff','g_diff1';...
                '.*','*';...
                './','/';...
                '.^','^'};
            %these are changed back by mupad
            renameopers2={'mod','_mod';...
                '=',':=';...
                'asin','arcsin';...
                'acos','arccos';...
                'atan','arctan';...
                'real','Re';...
                'imag','Im';...
                'log','ln';...
                '==','=';...
                '~',' not ';...
                '&&',' and ';...
                '||',' or ';...
                '&',' and ';...
                '|',' or '};
            %  [~,~,IB] = intersect(result.allvars,renameopers(:,1)');
            result.renamewords=[result.renamewords,transpose(renameopers(:,1))];
            result.rename2words=[result.rename2words,transpose(renameopers(:,2))];
            if hasiif||hasminmax
                result.renamewords(end+1)={'iif'};
                result.rename2words(end+1)={'piecewise'};
            end
            for j=length(obj):-1:1
                fs1=obj(j).fs;
                types1=obj(j).types;
                fsnumb=find(types1==obj(j).vr.number&~cellfun(@isempty,regexp(fs1,'[.]$','once')));
                for i=1:length(fsnumb)
                    fs1{fsnumb(i)}=[fs1{fsnumb(i)} '0'];
                end
                fs1=fs1(types1~=obj(j).vr.comment);
                for i=1:size(renameopers2,1)
                    fs1(strcmp(renameopers2{i,1},fs1))=renameopers2(i,2);
                end
                for i=1:length(result.renamewords)
                    fs1(strcmp(result.renamewords{i},fs1))=result.rename2words(i);
                end
                if isempty(fs1)
                    obj(j)=[];
                else
                    fsemicolon=find(strcmp(fs1,';'));
                    if ~isempty(fsemicolon)&&(fsemicolon(end)==length(fs1))
                        fs1{fsemicolon(end)}='';
                        fs1(fsemicolon(1:end-1))={','};
                    end
                    obj(j).fs=fs1;
                    obj(j).equation=sprintf('%s',obj(j).fs{:});
                    obj(j).types=[];
                    obj(j).structure=[];
                end
            end
        end
        function obj1 = runmupad(obj,comm)
            %RUNMUPAD - runs the equation in mupad (you can use the muPAD commands that return equations)
            %obj1 = runmupad(obj,comm) - comm should include % as placeholder for the equation, for instance 'diff(%,x)'
            %
            %example:
            %runmupad(parsed_equation('(Z*e*g)/(A+h)-A*Z*e*g*1/(A+h)^2'),'Simplify(%,Steps=1000)')
            %
            if ~i_hastoolbox('symbolic')
                error('parsed_equation:mupad','No symbolic toolbox');
            end
            [obj2,res] = mupadsyntax(obj);
            comm=sprintf('generate::MATLAB(%s)',strrep(comm,'%','%s'));
            eq=char(obj2);
            isvectorized= strcontains(eq,'.*')|| strcontains(eq,'.^')|| strcontains(eq,'./');
            s=char(evalin(symengine,sprintf(comm,eq)));
            s1=mupad2matlab(s,res);
            s1=regexprep(s1,'[*]1.0\>',''); %remove *1.0
            s1=regexprep(s1,'(?<=\<[1-9])[.][0]\>',''); %replace 1.0 with 1 etc
            
            if ~isvectorized
                s1=strrep(s1,'.*','*');
                s1=strrep(s1,'./','/');
                s1=strrep(s1,'.^','^');
            end
            obj1=parsed_equation(s1);
        end
        
        function [allunits, allvars, iter] = analyseunits(eqs, vars, units, debug)
            if nargin < 4
                debug = false;
            end
            if debug
                disp(eqs);
                line = 1;
                for i = 1:length(eqs)
                    s = structure2str(eqs(i).fs, eqs(i).structure);
                    s = str2cell(s);
                    for k = 1:length(s)
                        fprintf('%d:  %s\n', line, s{k});
                        line = line + 1;
                    end
                end
            end
            if ~isa(units, 'varunit')
                units = varunit(units);
            end
            allvars = transpose(symvar(eqs));
            for i = length(allvars):-1:1
                allunits(i) = varunit;
            end
            for i = 1:length(vars)
                ndx = find(strcmp(vars{i}, allvars));
                allvars(ndx) = vars(i);
                allunits(ndx) = units(i);
            end
            structs = cell(size(eqs));
            tmpvars = cell(size(eqs));
            tmpunits = cell(size(eqs));
            for i = 1:length(eqs)
                s = eqs(i).structure;
                %for efficiency some presets
                for j = 1:length(s)
                    s(j).cons = nan(size(s(j).args));
                    for k = 1:length(s(j).args)
                        s(j).cons(k) = str2double(s(j).args{k});
                    end
                    s(j).evaluated = false;
                end
                structs{i} = s;
                tmpv = cell(size(structs{i}));
                for j = length(tmpv):-1:1
                    if ~strncmp(structs{i}(j).leftvar, '#', 1)
                        tmpv{j} = [];
                    else
                        tmpv{j} = structs{i}(j).leftvar;
                    end
                end
                tmpvars{i} = tmpv;
                clear('tmpu','variable')
                tmpu(length(tmpv)) = varunit;
                for j = 1:length(tmpv) - 1
                    tmpu(j) = varunit;
                end
                tmpunits{i} = tmpu;
            end
            anythingchanged = true;
            iter = 0;
            while anythingchanged&&iter < 1000
                anythingchanged = false;
                iter = iter + 1;
                line = 1;
                for i = 1:length(eqs)
                    s = structs{i};
                    for j = 1:length(s)
                        if ~s(j).evaluated
                            con = [];
                            un = varunit;
                            for k = length(s(j).args):-1:1
                                un(k + 1) = getunit(s(j).args{k}, tmpvars{i}, tmpunits{i}, allvars, allunits);
                                con(k + 1) = s(j).cons(k);
                            end
                            con(1) = NaN;
                            un(1) = getunit(s(j).leftvar, tmpvars{i}, tmpunits{i}, allvars, allunits);
                            %                 con(1) = str2double(s(j).leftvar);
                            con(imag(con) > 0) = NaN;
                            [un2, changed] = solveoper(un, s(j).oper, con, eqs,[i,s(j).charnr]);
                            s(j).evaluated = ~any(un2.isundefined);
                            if any(changed)
                                anythingchanged = true;
                                for k = length(s(j).args):-1:1
                                    [tmpunits,allunits] = setunit(s(j).args{k}, un2(k + 1), i, tmpvars, tmpunits, allvars, allunits);
                                    if debug&&changed(k+1)
                                        fprintf('Iter %d:  line %d changed %s to %s\n', iter, line, s(j).args{k}, un2(k + 1).unit);
                                    end
                                end
                                [tmpunits,allunits] = setunit(s(j).leftvar, un2(1), i, tmpvars, tmpunits, allvars, allunits);
                                if debug&&changed(1)
                                    fprintf('Iter %d:  line %d changed %s to %s\n', iter, line, s(j).leftvar, un2(1).unit);
                                end
                            end
                        end
                        line = line + 1;
                    end
                    structs{i} = s;
                end
            end
            if nargout == 0
                fprintf('%d iters:\n', iter);
                for i = 1:length(allvars)
                    fprintf('%s [%s]\n', allvars{i}, char(allunits(i)));
                end
            end
        end
        
        function [pars,js] =  getfunctionpars(obj, j)
            %j=index of function name
            %brackets
            obj.types(strcmp(obj.fs, '(')) = obj.vr.brack1;
            obj.types(strcmp(obj.fs, ')')) = obj.vr.brack2;
            obj.types(strcmp(obj.fs,',')) = obj.vr.comma;
            obj.types(strcmp(obj.fs, '; ')) = obj.vr.semicolon;
            while j>0&&(obj.types(j)~=obj.vr.funcname)
                j=j-1;
            end
            j = j + 2;
            if obj.types(j)==obj.vr.brack2
                %no arguments
                pars={};
                js=[j j-1];
                return;
            end
            poplevel = 1;
            pars = cell(5, 1);
            js = zeros(5, 2);
            ip = 1;
            while (j < length(obj.fs))&&poplevel > 0
                k = 0;
                while (j+k < length(obj.types)) && poplevel~=0 && ~(poplevel==1&&any(obj.types(j+k)== obj.vr.comma))
                    if obj.types(j + k) == obj.vr.brack1
                        poplevel = poplevel + 1;
                    end
                    k = k + 1;
                    if any(obj.types(j + k) == [obj.vr.semicolon, obj.vr.brack2])
                        poplevel = poplevel - 1;
                    end
                end
                js(ip, :) = [j, j + k - 1];
                pars{ip} = sprintf('%s', obj.fs{js(ip, 1):js(ip, 2)});
                ip = ip + 1;
                j = j + k + 1;
            end
            pars = pars(1:ip - 1);
            js = js(1:ip - 1, :);
        end
        
        function res = eval(obj)
            %eval for parsed_equation
            res = evalin('caller',char(obj));
        end
        
        function res = evalin(caller, obj)
            res = evalin(caller, char(obj));
        end
        
        function str = cell(obj)
            str = cell(length(obj), 1);
            for i = 1:length(obj)
                str{i} = obj(i).equation;
            end
        end
        function s = sym(obj)
            s=sym(char(obj));
        end
        function str = char(obj)
            if length(obj) > 1
                c = cell(obj);
                str = sprintf('%s\n', c{:});
            else
                str = obj.equation;
            end
        end
        
        function str = latex(obj)
            str = str2latex(obj);
        end
        
        function viewlatex(obj)
            str2latex(obj);
        end
        
        function s1 = symvar(obj)
            s1 = {};
            for i = 1:length(obj)
                s1=[s1 obj(i).fs(obj(i).types == parsed_equation.vr.variable)];
            end
            s1 = unique(s1);
        end
        
        function str = msword(obj)
            [~, str] = str2latex(obj);
        end
        
        %         function [str, fs]  = div2power(obj, opersinfo)
        %             if nargin < 2
        %                 opersinfo = parsed_equation.operinfo;
        %             end
        %             if length(obj) > 1
        %                 str = cell(length(obj), 1);
        %                 fs = cell(length(obj), 1);
        %                 for i = 1:length(obj)
        %                     [str{i}, fs{i}] = div2power(obj(i), opersinfo);
        %                 end
        %                 return;
        %             end
        %             s = obj.structure;
        %             f = find(strcmp({s.oper(:)}, '/'));
        %             nr = length(s) + 1;
        %             for j = length(f):-1:1
        %                 arg2 = s.args{f(j)}{2};
        %                 s.args{f(j)}{2} = sprintf('#%03d', nr);
        %                 s.opers=[s.opers(1:f(j)-1) {'^', '*'} s.opers(f(j)+1:end)];
        %                 s.operstype = [s.operstype(1:f(j) - 1) obj.vr.oper s.operstype(f(j):end)];
        %                 s.leftvar = [s.leftvar(1:f(j) - 1) {sprintf('#%03d', nr)} s.leftvar(f(j):end)];
        %                 a = {arg2, '-1'};
        %                 s.args = [s.args(1:f(j) - 1) {a} s.args(f(j):end)];
        %                 nr = nr + 1;
        %             end
        %             disp(structure2str(obj.fs, s))
        %             fs = makeequation({}, s, length(s.opers), opersinfo);
        %             str = sprintf('%s', fs{:});
        %         end
        function [obj] = changevar(obj, oldvar,newvar)
            if length(obj) > 1
                for i = 1:length(obj)
                    obj(i) = changevar(obj(i), oldvar,newvar);
                end
                return;
            end
            ndx = strcmp(oldvar, obj.fs);
            if any(ndx)
                fs1 = obj.fs;
                fs1(ndx) = {newvar};
                obj.equation = sprintf('%s', fs1{:});
                obj.fs = fs1;
            end
        end
        function [obj] = maxbrackets(obj, opersinfo)
            %maxbrackets - add all unnecessary brackets according to the
            %precedence rules.
            %you can define your own precedence rules by adapting the
            %structure that is received from parsed_equation.operinfo
            %
            %use:
            %obj.maxbrackets(operinfo)
            if nargin < 2
                opersinfo = parsed_equation.operinfo;
            end
            if length(obj) > 1
                for i = 1:length(obj)
                    obj(i) = maxbrackets(obj(i), opersinfo);
                end
                return;
            end
            s = obj.structure;
            if ~isempty(s)
                opersinfo.prec = [];
                obj = parsed_equation(s, opersinfo);
            end
        end
        function [obj] = minbrackets(obj, opersinfo)
            %minbrackets - remove unnecessary brackets according to the
            %precedence rules.
            %you can define your own precedence rules by adapting the
            %structure that is received from parsed_equation.operinfo
            %
            %use:
            %obj.minbrackets(operinfo)
            if nargin < 2
                opersinfo = parsed_equation.operinfo;
            end
            if length(obj) > 1
                for i = 1:length(obj)
                    obj(i) = minbrackets(obj(i), opersinfo);
                end
                return;
            end
            s = obj.structure;
            if ~isempty(s)
                obj = parsed_equation(s, opersinfo);
            end
        end
        
        function str = expanded(obj)
            lvar = 'tmp';
            while any(strncmp(obj.fs, lvar, length(lvar)))
                lvar = [lvar '_'];
            end
            s = obj.structure;
            str = structure2str(obj.fs, s);
            str = strrep(str, '#', lvar);
        end
        
        function disp(obj)
            disp(char(obj));
        end
        
        function [haserror,msg]=test(obj)
            fs1 = obj.fs;
            haserror=false;
            f=find(strcmp(fs1, '='));
            if ~isempty(f)
                %remove leftside before equal sign
                fs1 = fs1(f + 1:end);
                obj.equation=sprintf('%s',fs1{:});
            end
            fs1=obj.addbrackets(false);
            s = obj.equation;
            disp('testing addbrackets');
            s1 = sprintf('%s', fs1{:});
            disp(s1); % s1 is used later, no longer a warning in R2018a
            haserror = dotest(s, s1);
            if ~haserror
                disp('testing min/maxbrackets');
                if ~isempty(f)
                    s=strrep(s,' ',''); %else it can be seen as a function command line
                    obj1 = parsed_equation(s);
                    s1 = obj1.minbrackets.equation;
                    s2 = obj1.maxbrackets.equation;
                else
                    s1 = obj.minbrackets.equation;
                    s2 =  obj.maxbrackets.equation;
                end
                disp(s1);
                haserror = dotest(s, s1);
                disp(s2);
                haserror = haserror||dotest(s, s2);
                if ~haserror
                    disp('testing expanded');
                    s1 = obj.expanded;
                    disp(s1);
                    dotest(s, s1, obj.structure(end).leftvar);
                end
            end
        end
        function s = vectorize(obj)
            for i = 1:length(obj)
                fs1 = obj(i).fs;
                ndx = strcmp(fs1, '*');
                fs1(ndx) = {'.*'};
                ndx = strcmp(fs1, '^');
                fs1(ndx) = {'.^'};
                ndx = strcmp(fs1, '/');
                fs1(ndx) = {'./'};
                ndx = strcmp(fs1, '&&');
                fs1(ndx) = {'&'};
                ndx = strcmp(fs1, '||');
                fs1(ndx) = {'|'};
                obj(i).fs = fs1;
                obj(i).equation = sprintf('%s', obj(i).fs{:});
                obj(i).structure = [];
            end
            s = char(obj);
        end
        function s = errorat(obj, pos)
            if length(pos) == 1
                pos = [1 pos];
            end
            s=sprintf('%s\n%s^\n',char(obj(pos(1))),repmat(' ',1,pos(2)-1));
        end
        
        function [fs, types] = addbrackets(obj, usebrackcode, precedence)
            %fs is not updated
            if (nargin < 2)
                usebrackcode = true;
            end
            if usebrackcode
                qbracks={'«','»'};
            else
                qbracks={'(',')'};
            end
            %precedence{1}={'(',')','«','»'}; brackets have precedence anyway
            if nargin<3
                precedence = parsed_equation.operinfo;
            end
            fsunary =  addunarycode(obj.fs, qbracks);
            fs = fsunary;
            fbr = find(strcmp(fs, '['));
            if ~isempty(fbr)
                fis=find(strcmp(fs, '='), 1);
                if ~isempty(fis)
                    fbr = fbr(fbr > fis);
                end
                for j = length(fbr):-1:1
                    fs = [fs(1:fbr(j) - 1) qbracks(1) fs(fbr(j):end)];
                end
                fbr = find(strcmp(fs, ']'));
                if ~isempty(fis)
                    fbr = fbr(fbr > fis);
                end
                for j = length(fbr):-1:1
                    fs = [fs(1:fbr(j)) qbracks(2) fs(fbr(j) + 1:end)];
                end
            end
            isvar = false(size(fs));
            for j = 1:length(fs)
                if ~isempty(fs{j})
                    isvar(j)=isletter(fs{j}(1))||(fs{j}(1)=='_');
                end
            end
            
            isfun = isvar;
            %isvarfun = isvar;
            for i = length(fs):-1:1
                if isfun(i)
                    if (i==length(fs))||~strcmp(fs{i+1}, '(')
                        isfun(i) = false;
                    else
                        isvar(i) = false;
                        j = popbrack(fs, i + 1, 1);
                        fs = [fs(1:i - 1) qbracks(1) fs(i:j) qbracks(2) fs(j + 1:end)];
                    end
                end
            end
            for i = 1:10
                fs=addbrackoper(fs, precedence.opers(precedence.prec==i), i == 2, qbracks, i==5);
            end
            fs=cellrep(fs,'+$','+');
            fs=cellrep(fs,'-$','-');
            if nargout > 1
                types = updatetypes(fs,obj.vr);
            end
        end
    end
    methods (Static)
        function operinf = operinfo
            %     MATLAB precedence rules, used by parsed_equation (from highest to lowest)
            %      1:   ' .' .^ ^
            %      2:   unary+ unary- ~
            %      3:   .* ./ .\ * / \
            %      4:   + -
            %      5:   :
            %      6:   < <= > >= == ~=
            %      7:   &
            %      8:   |
            %      9:   &&
            %      10:   |
            oper.opers={'''','.''','.^','^',...
                '+$','-$','~',...
                '.*','./','.\','*','/','\',...
                '+','-',...
                ':', ...
                '<','<=','>','>=','==','~=',...
                '&','|','&&','||'};
            oper.prec = [1, 1, 1, 1, ...
                2, 2, 2, ...
                3, 3, 3, 3, 3, 3, ...
                4, 4, ...
                5, ...
                6, 6, 6, 6, 6, 6, ...
                7, 8, 9, 10];
            oper.type = zeros(size(oper.prec)) + parsed_equation.vr.oper;
            oper.type(oper.prec == 2)=parsed_equation.vr.unary_oper_pre;
            oper.type(1) = parsed_equation.vr.unary_oper_post;
            oper.type(2) = parsed_equation.vr.unary_oper_post;
            if nargout == 0
                opers = oper.opers;
                ndx = strcmp('+$', opers);
                opers{ndx} = 'unary+';
                ndx = strcmp('-$', opers);
                opers{ndx} = 'unary-';
                disp('MATLAB precedence rules, used by parsed_equation (from highest to lowest)')
                for i = 1:max(oper.prec)
                    s=sprintf('%s ', opers{oper.prec == i});
                    fprintf('   %d:   %s\n', i, s);
                end
                disp('see also <a href="matlab:doc precedence">doc precedence</a>');
            else
                operinf = oper;
            end
        end
    end
end
%END OF CLASS DEFINITION

function str = structure2str(fs, s)
vr = parsed_equation.vr;
str = '';
for i = 1:length(s)
    if ~isempty(s(i).oper)
        sargs = sprintf('%s,',s(i).args{:});
        if ~isempty(sargs)
            sargs = sargs(1:end - 1);
        end
        switch s(i).opertype
            case vr.unary_oper_pre
                str=sprintf('%s%s = %s%s;\n', str, s(i).leftvar, s(i).oper, sargs);
            case vr.unary_oper_post
                str=sprintf('%s%s = %s%s;\n', str, s(i).leftvar, sargs, s(i).oper);
            case vr.oper
                str=sprintf('%s%s = %s%s%s;\n', str, s(i).leftvar, s(i).args{1}, s(i).oper, s(i).args{2});
            case vr.funcname
                str=sprintf('%s%s = %s(%s);\n', str, s(i).leftvar, s(i).oper, sargs);
        end
    end
end
if (length(s)==1)&&isempty(s.oper)
    f=find(strcmp(fs, '='));
    if isempty(f)
        f = 0;
    end
    str=sprintf('%s%s = %s;\n', str, s(1).leftvar, sprintf('%s', fs{f+1:end}));
end
% str = strrep(str, '#', lvar);
end

function fs = makeequation(fs, s, i, opers)

    function fs = makeeq(fs, s, i, opers)
        %
        mypreced = getprecedence(s(i).oper, s(i).opertype, opers);
        arg = cell(size(s(i).args));
        for j = 1:length(s(i).args)
            arg{j} = s(i).args(j);
            if strncmp(arg{j}, '#', 1)
                k = find(strcmp(arg{j}, {s.leftvar}));
                if ~isempty(k)
                    arg{j} = makeeq({}, s, k, opers);
                    if s(i).opertype ~= parsed_equation.vr.funcname
                        if isempty(opers.prec)
                            arg{j}=[{'('},arg{j},{')'}];%sprintf('(%s)',arg{j});
                        else
                            preced1 = getprecedence(s(k).oper, s(k).opertype, opers);
                            if preced1 > mypreced||(preced1==mypreced&&j>1) %larger value is lower precedence
                                %same precedence is processed from left to right unless the user uses brackets
                                %This order remains unchanged here even for
                                %commutative operators (+*)
                                arg{j}=[{'('},arg{j},{')'}];
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(s(i).opertype)
            switch s(i).opertype
                case parsed_equation.vr.unary_oper_pre
                    fs = [fs, s(i).oper, arg{1}];
                case parsed_equation.vr.unary_oper_post
                    fs = [fs,arg(1), s(i).oper];
                    % str=sprintf('%s%s%s', str,  arg{1}, s.opers{i});
                case parsed_equation.vr.oper
                    fs = [fs, arg{1}, s(i).oper, arg{2}];
                    %    str=sprintf('%s%s%s%s', str, arg{1}, s.opers{i}, arg{2});
                case parsed_equation.vr.funcname
                    n = 0;
                    for m = 1:length(arg)
                        n = n + 1 + length(arg{m});
                    end
                    sargs1 = cell(1, n - 1);
                    n = 1;
                    for m = 1:length(arg)
                        sargs1(n:n + length(arg{m}) - 1) = arg{m};
                        n = n + length(arg{m}) + 1;
                        if m < length(arg)
                            sargs1(n - 1) = {','};
                        end
                    end
                    fs=[fs, {s(i).oper,'('} sargs1, {')'}];
                    %     sargs1 = sprintf('%s,',arg{:});
                    %         if ~isempty(sargs1)
                    %             sargs1 = sargs1(1:end - 1);
                    %         end
                    %         str=sprintf('%s%s(%s)', str,  s.opers{i}, sargs1);
            end
        else
            if ~isempty(s(i).args)
                fs = [fs s(i).args];
            end
        end
    end
%check for circular arguments (else MATLAB may hang in makeeq)
for i = 1:length(s)
    for j1 = 1:length(s(i).args)
        arg1 = s(i).args{j1};
        if ~isempty(arg1)&&arg1(1)=='#'
            k1 = find(strcmp(arg1, {s.leftvar}));
            if k1 >= i
                error('parsed_equation:syntaxerror','Syntax error in the equation, double operators?, oper=%s',s.opers{i});
            end
        end
    end
end
fs = makeeq(fs, s, i, opers);
if ~strcmp(s(end).leftvar, 'ans')
    fs=[s(end).leftvar {'='} fs];
end
end

function preced = getprecedence(oper, opertype, operinfo)
if isempty(operinfo.prec)
    preced = [];
else
    if ~isempty(oper)&& (opertype==parsed_equation.vr.unary_oper_pre)&&any(strcmp(oper,{'+','-'}))
        oper = [oper '$'];
    end
    ndx = strcmp(oper, operinfo.opers);
    if ~any(ndx)
        preced = 0;
    else
        preced = operinfo.prec(ndx);
    end
end
end



function types = updatetypes(fs,vr)

fsunary =  addunarycode(fs);
fis=find(strcmp(fsunary, '='));
if isempty(fis)
    fis = 0;
end
opers={'.*','./','.\','*','/','\','.^','^','+','-',':','<','<=','>','>=','==','~=','&','|','&&','||'};
keywords={'pi','inf','Inf','nan','NaN','eps','if','case','else','elseif','function','while','for','switch','otherwise'};
kwtype = [vr.constant, vr.constant, vr.constant, vr.constant, vr.constant, vr.constant, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw];
operletter='.*/\^+-:<>=~&|';
unary_opers_pre={'+$','-$','~'};
unary_opers_post={'''','.'''};
types =  zeros(size(fsunary));
for i = 1:length(unary_opers_pre)
    types(strcmp(fsunary, unary_opers_pre{i})) = vr.unary_oper_pre;
end
for i = 1:length(unary_opers_post)
    types(strcmp(fsunary, unary_opers_post{i})) = vr.unary_oper_post;
end
isend = strcmp(fs, 'end');
if any(isend)
    types(isend) = vr.kw;
    fbr1 = find(strcmp(fs, '('), 1);
    if ~isempty(fbr1)
        fbr2=find(strcmp(fs,')'),1,'last');
        f = find(isend);
        ndx = f > fbr1 & f < fbr2;
        types(f(ndx)) = vr.variable;
    end
end
%       ndx = strcmp(fs, '(') | strcmp(fs, qbracks{1});
%       types(ndx) = vr.brack1;
%       ndx = strcmp(fs, ')') | strcmp(fs, qbracks{2});
%       types(ndx) = vr.brack2;
%       feq=strcmp(fs, '=');
%       types(feq) = vr.assign;
if any(strcmp(' $', fsunary)|strcmp(' ', fsunary))&&~any(strcmp(fsunary{1},keywords))
    types(1) = vr.funcname;
else
    for j = 1:length(fsunary)
        if ~isempty(fsunary{j})&&~any(types(j)==[vr.kw vr.variable])
            ch = fsunary{j}(1);
            if isletter(ch)||(ch=='_')
                if (j < fis(1))||(j==length(fsunary))||~strcmp(fsunary{j+1}, '(')
                    s = strcmp(fsunary{j}, keywords);
                    if any(s)
                        types(j) = kwtype(s);
                    else
                        types(j) = vr.variable;
                    end
                else
                    types(j) = vr.funcname;
                end
            elseif strcontains('1234567890', ch)
                types(j) = vr.number;
            elseif ch == '%'
                types(j) = vr.comment;
            elseif strcontains(operletter, ch)&&any(strcmp(fsunary{j},opers))
                types(j) = vr.oper;
            elseif ch == '.'
                if isnan(str2double(fsunary{j}))
                    types(j) = vr.field;
                else
                    types(j) =  vr.number;
                end
            end
        end
    end
end
end
function un = getunit(cunit, tmpvars, tmpunits, allvars, allunits)
if strncmp(cunit, '#', 1)
    f1 = find(strcmp(tmpvars, cunit));
    f2 = [];
else
    f1 = [];
    f2 = find(strcmp(allvars, cunit));
end
if isempty(f1)&&isempty(f2)
    un = varunit;
elseif isempty(f2)
    un = tmpunits(f1(1));
else
    un = allunits(f2(1));
end
end
function [tmpunits,allunits] = setunit(cunit, un, i, tmpvars, tmpunits, allvars, allunits)
%[tmpunits,allunits]=setunit(s(j).leftvar, un(1), i, tmpvars, tmpunits, allvars, allunits);
if strncmp(cunit, '#', 1)
    f1 = find(strcmp(tmpvars{i}, cunit));
    if ~isempty(f1)
        tmpunits{i}(f1(1)) = un;
    end
else
    f2 = find(strcmp(allvars, cunit));
    if ~isempty(f2)
        allunits(f2(1)) = un;
    end
end
end

function fsnew  = parseeq(fs1)

%disp('parseeq intern')

% %remove ... and parse equation
% if ischar(fs1)
%     fs1=regexprep(fs1,'(?<!%.*)[.][.][.]$',' ');
%     fs1 = doparseeq1(fs1);
% %    fs1 = doparseeq(fs1);
% elseif iscell(fs1)
%     %only if last in the line the dots should be replaced
%     %but not if it is in a comment
%     fddd=regexp(fs1,'(?<!%.*)[.][.][.]$');
%     for i=length(fddd):-1:1
%         if ~isempty(fddd{i})
%             fs1{i}=strrep([fs1{i} ' ' fs1{i+1}],'...','');
%             fs1(i+1)={'%'};
%         end
%     end
%      fs1=doparseeq1(fs1);
% %    for i=1:length(fs1)
% %        fs1{i} = doparseeq(fs1{i});
% %    end
% end
fs1=doparseeq1(fs1);
fsnew=fs1;
if ~isempty(fsnew)&&iscell(fsnew{1})
    for i=1:length(fsnew)
        fsnew{i}=afterparse(fsnew{i});
    end
else
    fsnew=afterparse(fsnew);
end
end
function fs=afterparse(fs)
f = strcmp(fs, ' ');

if any(f)
    f2 = [(f(1:end - 1) + f(2:end)) > 1 false];
    fs = fs(~f2);
    f = strcmp(fs, ' ') ;
    
    isvar = false(size(fs));
    isnum = isvar;
    isfun = isvar;
    for j = 1:length(fs)
        if ~isempty(fs{j})
            ch = fs{j}(1);
            if isletter(ch)||(ch=='_')
                if (j==length(fs))||~strcmp(fs{j+1}, '(')&&~any(strcmp(fs{j},{'if','case','else','elseif','end','function'}))
                    isvar(j) = true;
                else
                    isfun(j) = true;
                end
            else
                %           isnum(j) = ~isempty(strfind('1234567890', ch));
                isnum(j) = any('1234567890'==ch);
            end
        end
    end
    isvar1 = isvar | isfun | isnum;
    f2=strcmp(fs,'(')|strcmp(fs,'[')|strcmp(fs,'{');
    f3=strcmp(fs,')')|strcmp(fs,']')|strcmp(fs,'}');
    for j = length(fs)-1:-1:2
        if f(j)
            if (isvar1(j - 1)||(f3(j - 1)))&&isvar1(j+1)||(f3(j - 1)&&f2(j+1))
                f(j) = false;
            end
        end
    end
    fs = fs(~f);
end

if nargout == 0
    fs=cellrep(fs,' $',' ');
    fprintf('%s\n',strrep(sprintf('%s',fs{:}),' ',char(183)));
else
    fs = cellrep(fs,' $',' ');
end

end


% ks = find(strcmp(sfrom, fs));
% for i = 1:length(ks)
%    fs{ks(i)} = sto;
% end
function fsnew = addunarycode(fs, ~)
opers={'=',',','.^','^','+$','-$','~','.*','./','.\','*','/','\','+','-','~',':','<','<=','>','>=','==','~=','&','|','&&','||'};
unaries={'-','+'};
fsnew = fs;
for i = 1:2
    f = find(strcmp(fs, unaries{i}));
    for j = 1:length(f)
        if f(j) == 1
            fsnew{f(j)} = sprintf('%s$', unaries{i});
        elseif any(strcmp(fs{f(j) - 1}, [opers, {'(', '«'}]))
            fsnew{f(j)} = sprintf('%s$', unaries{i});
        end
    end
end
end

function f = popbrack(s,istart, istep)
%s and ivar should be created with parseq
%use istart=find(fpos>=istart,1) to create an istart from a position in the
%original string

if istep > 0 %with powers unaries have precedence if they are after ^
    while istart<length(s)&&any(strcmp(s{istart},{'-$','+$'}))
        istart = istart + 1;
    end
end
isvar = false(size(s));
isfield = false(size(s));
for i = 1:length(s)
    if ~isempty(s{i})
        isvar(i)=isletter(s{i}(1))||(s{i}(1)=='_');
        if (s{i}(1)=='.')&&(length(s{i}) > 1)
            isfield(i) = isletter(s{i}(2));
        end
    end
end
if length(s) < 2
    f = istart;
    return;
end
f = istart;
found = 0;
brack = 0;
if istep > 0
    brack1={'(','«'};
    brack2={')','»'};
else
    brack2={'(','«'};
    brack1={')','»'};
end
while (f > 0)&&(f<=length(s))&&strcmp(s{f}, ' ')
    f = f + istep;
end
isbrack1 = strcmp(s, brack1{1}) | strcmp(s, brack1{2});
isbrack2 = strcmp(s, brack2{1}) | strcmp(s, brack2{2});

while ~found
    while ((f > 0)&&(f<=length(s)))&&(isvar(f)||isfield(f)||((s{f}(1)=='.')||(s{f}(1)>='0')&&(s{f}(1)<='9')))
        f = f + istep;
    end
    if (f<=0)||(f > length(s))
        found = true;
    else
        %     s1 = s(f);
        prec = false;
        if isbrack1(f) %any(strcmp(s1, brack1))
            brack = brack + 1;
        elseif isbrack2(f) %any(strcmp(s1, brack2))
            brack = brack - 1;
            if brack == 0 &&~prec
                prec = f < length(s); %&&any(strcmp(precedence, s{f+1}));
                if ~prec
                    f = f + istep;
                end
            end
        end
        if brack <= 0 &&~prec
            found = true;
            f = f - istep;
        else
            f = f + istep;
        end
    end
end
if f == 0
    f = 1;
end
if f > length(s)
    f = length(s);
end
%if (istep==-1)&&(f > 1)&&(s(f)=='(')&&isvarletter(s(f-1)) %exp(x)^exp(x) left side solved
if (istep==-1)&&(f > 1)&&(strcmp(s{f}, '(')&&isvar(f-1)) %exp(x)^exp(x) left side solved
    f = f - 1;
    while (f > 0)&&isvar(f)
        f = f - 1;
    end
    f = f + 1;
end
end

function fs = addbrackoper(fs, opers, isunary, qbracks, isthreeoper)
f = false(size(fs));
for i = 1:length(opers)
    f = f | strcmp(fs, opers{i});
end
i = find(f, 1);
while ~isempty(i)
    isunary2=any(strcmp(fs{i},{'''','.'''}));
    f(i) = false;
    if ~isunary&&(i > 1)
        j = popbrack(fs,i - 1, - 1);
    else
        j = i;
    end
    if isunary %this is to handle double opers: --a
        while (i+1 < length(fs))&&any(strcmp(fs{i+1}, opers))
            i = i + 1;
        end
    end
    if ~isunary2&&i < length(fs)
        j2 = popbrack(fs, i + 1, 1) + 1;
    else
        j2 = i + 1;
    end
    if isthreeoper&&j2 < length(fs)&&any(strcmp(opers, fs{j2}))
        f(j2) = false;
        j2 = popbrack(fs, j2 + 1, 1) + 1;
    end
    %   if isunary2 %this is to handle double opers:''''
    %      while (i-1 >0)&&any(strcmp(fs{i+1}, {'''','.'''}))
    %         i = i - 1;
    %      end
    %   end
    if j < 1
        j = 1;
    end
    if j2 > length(fs) + 1
        j2 = length(fs) + 1;
    end
    if j==1||~(strcmp(fs{j-1},'(')&&(strcmp(fs{j2},')')))&&(isunary||isunary2||j2>j+2)
        %    fs=[fs(1:j-1) {sprintf('«(%d)',i)} fs(j:j2-1) {sprintf('»(%d)',i)} fs(j2:end)];
        fs = [fs(1:j - 1) qbracks(1) fs(j:j2 - 1) qbracks(2) fs(j2:end)];
        f = [f(1:j - 1) false f(j:j2 - 1) false f(j2:end)];
    end
    i = find(f, 1);
end
end
function fss = doparseeq1(eqs,pattern)
if nargin<2
    pattern=['([%].*)|',... comment
        '(["][^"]*["])|',... %double quote string
        '([.][.][.]|[.][.])|'... %dots or .. or ...
        '([.$][\*/\\^''])|'...%dot opers
        '([.]?[0-9]+[.]?[0-9]*([eE][-+]?)?[0-9]*)|'... %numbers
        '[@][A-Za-z_]+[A-Za-z0-9_]*|',...%function handles
        '([.][A-Za-z_]+[A-Za-z0-9_]*)|',...%fields
        '([A-Za-z_][A-Za-z0-9_#]*)|'...%symbols
        '([&|=~<>][&|=<>])|'...%longer opers
        '([ .~><\(\)^\*\\/\+='':;,&|%!{}\[\]#@"])|[-\n]'];%any other char
end
if ischar(eqs)
    themodel=strtrim(regexprep({eqs},'[ \t]*',' '));
else
    themodel=strtrim(regexprep(eqs,'[ \t]*',' '));
end
%split roughtly
%it seems not possible to let 2.*pi translate to "2" ".*" "pi"
% replace first dot to $ (not very efficient though)
themodel1=regexprep(themodel,'[\.](?=[\*/\\^''])','$');
%magical function!
[ss1,fss,splits]=regexp(themodel1,pattern ,'start','match','split');
for i=1:length(fss)
    fss{i}=regexprep(fss{i},'[$](?=[\*/\\^''])','.'); %repair
end
fsymbols=regexp(themodel, '([A-Za-z_][A-Za-z0-9_#]*)','start');%symbols
%fdotop=regexp(themodel,'([%].*)|[.][*/\\^'']','start');%dot opers to solve the problem of 2.* ->'2','.*'
faccent=strfind(themodel,'''');
feq=strfind(themodel,'=');
naccent=cellfun('length',faccent);
%all operators with length 2
%
for i=1:length(themodel)
    fs1=fss{i};
    if any(~cellfun('isempty',splits{i}))
        fprintf('|%s',fss{i}{:});
        fprintf('\n\n');
        disp(themodel{i});
        fprintf('\n\n');
        s=sprintf('%s',splits{i}{:});
        error('parsed_equation:parse','unknown chararacter %s',s);
    end
    %     if ~isempty(fdotop{i}) %inefficient and probably not necessary
    %         %dot opers to solve the problem of 2.* ->'2','.*'
    %         wrongdot=setdiff(fdotop{i},ss1{i});
    %         if ~isempty(wrongdot)
    %             for k=length(wrongdot):-1:1
    %                 fndx=find(ss1{i}<wrongdot(k),1,'last');
    %                 s1=fs1{fndx};
    %                 s2=fs1{fndx+1};
    %                 fs1{fndx}=s1(1:end-1);
    %                 fs1{fndx+1}=[s1(end) s2];
    %                 ss1{i}(fndx+1)=ss1{i}(fndx+1)-1;
    %             end
    %         end
    %     end
    if (naccent(i)>0)
        if ~isempty(feq{i})&&feq{i}(1)>faccent{i}(1)
            naccent(i)=naccent(i)-1;
            hasis=true;
        else
            hasis=false;
        end
        if (naccent(i)>0)&&mod(naccent(i),2)==0
            s=themodel{i};
            if hasis
                s(faccent{i}(1))='*';
            end
            [ss1{i},fs1]=regexp(s,['\<[%].*|[''][^'']*['']|' pattern ],'start','match');
            if hasis
                fs1{find(strcmp(fs1,'*'),1)}='''';
            end
        end
    end
    if length(fsymbols{i})>=2&&fsymbols{i}(1)==ss1{i}(1)&&fsymbols{i}(2)==ss1{i}(3)&&strcmp(fs1{2},' ')
        %command style
        ss=intersect(fs1,{'if','elseif','for','while','case','switch'});
        if isempty(ss)
            %'([%].*)|[ ]|((?<=[''])[^'']*(?=[\'']))|([\[][^\]]*[\]])|([^ ''\[\]]*)','start','match')
            
            
            pattern1=['([%].*)'... %comment
                '|[ ]',... %space
                '|((?<=[''])[^'']*(?=[\'']))'... %string ' '
                '|([\[][^\]]*[\]])',... %string [  ]
                '|([^ ''\[;\]]*)']; %only
            [ss1{i},fs1]=regexp(themodel{i},pattern1,'start','match');
            if strncmp(fs1{end},'%',1)
                fs1=fs1(1:end-1);
                if ~isempty(fs1)&&strncmp(fs1{end},' ',1)
                    fs1=fs1(1:end-1);
                end
                %          else
                %              fs1=fs1;
            end
            ndx=strcmp(fs1,' ');
            fs1(ndx)={' $'};
        end
    end
    fss{i}=fs1;
end
if ischar(eqs)
    fss=fss{1};
end
end
function fs = doparseeq(s)
%breaks an equation into words (operands and variables etc)
if ~isempty(s)&&s(1)=='%' % do not analyse comment lines
    fs = {s};
    return;
end
s = strtrim(s);
fs = cell(1, length(s));
sletter = isletter(s);
k = 0;
i = 1;
while i <= length(s)
    j = i;
    breakch = false;
    isval = false;
    prevch = ' ';
    while (j<=length(s))&&~breakch
        ch = s(j);
        isletterch = sletter(j);
        if isletterch
            while j < length(s)&&sletter(j+1)
                j = j + 1;
                ch = s(j);
            end
        end
        if (i==j)&&((ch>='0')&&(ch<='9'))||(j<length(s)&&(s(j)=='.')&&((s(j+1)>='0')&&(s(j+1)<='9')))
            isval = true;
            hasdot = false;
        end
        breakch = strcontains(' .~><()^*\/+-='':;,&|%!{}[]', ch);
        if ch=='.' && (j+1<=length(s)) && sletter(j+1)
            breakch = false;
        end
        if (isletterch|| ch=='_' || (~isval && (ch>='0')&&(ch<='9'))) && (j+1<=length(s)) && (s(j+1)=='.')
            breakch = true;
            j = j + 1;
        end
        if k==1&&ch==' '&& ~any(strcmp(fs{k},{'if','elseif','for','while','case','switch'}))
            %analyse command syntax
            if j+1 > length(s)||s(j+1)~='='
                prevop = fs{k};
                if isletter(prevop(1))||(prevop(1)=='_') %first must be a function
                    l = length(s);
                    s1 = s;
                    while 1
                        s1=strrep(s1,'  ',' ');
                        l1 = length(s1);
                        if l1 == l
                            break;
                        end
                        l = l1;
                    end
                    s1=strrep(s1,' ','\n');
                    pars = str2cell(sprintf(s1));
                    if length(pars) > 1
                        if ~strcontains('=(',pars{2}(1))&&~any(strcmp(pars{2},{'*','-','.*','^','.^','''','/','./'}))
                            for ii = 1:length(pars)
                                if ii == 1
                                    fs = pars(1);
                                elseif ~isempty(pars{ii})
                                    if ii==length(pars)&&pars{ii}(end)==';'
                                        pars{ii} = pars{ii}(1:end - 1);
                                    end
                                    fs = [fs {' $', pars{ii}}];
                                end
                            end
                            return;
                        end
                    end
                end
            end
        end
        %command syntax
        if isval&&(((~hasdot)&&(ch=='.'))||(strcontains('-+',ch)&&(prevch=='E')))
            if ~((j < length(s)&& strcontains('*\/^''', s(j+1))))
                breakch = false;
            end
        end
        if ch == '.'
            hasdot = true;
        end
        if i==j && ch == '%' % start of a comment, rest is ignored
            j = length(s);
        elseif ~breakch
            prevch = upper(ch);
            j = j + 1;
        elseif i ~= j
            j = j - 1;
        end
    end
    if j > length(s)
        j = length(s);
    end
    op = s(i:j);
    if (k > 0)&&(length(op) == 1) && (length(fs{k}) <= 2)
        op3 = [fs{k}, op];
        if any(strcmp(op3,{'...','..','.*','.\','./','.^','.''','==','~=','!-','>=','<=','&&','||'}))
            k = k - 1;
            op = op3;
        end
    end
    k = k + 1;
    fs{k} = op;
    i = j + 1;
end
fs = fs(1:k);

faccent=find(strcmp('''',fs));
%check whether there are strings in the equation, a single ' is interpreted as transpose,
%a string cannot occur before =
if ~isempty(faccent)
    feq=find(strcmp('=', fs), 1);
    if ~isempty(feq)&&(faccent(1) < feq)
        faccent(1) = [];
    end
    nfaccent = length(faccent);
    if nfaccent > 1 && (mod(nfaccent, 2)==0)%iseven
        for i = 1:2:nfaccent
            fs{faccent(i)} = sprintf('%s', fs{faccent(i):faccent(i + 1)});
            for k = faccent(i + 1):-1:faccent(i) + 1
                fs{k} = '';
            end
        end
    end
    fs = fs(~strcmp(fs, ''));
end
end
function fs =  cellrep(fs, sfrom, sto)
ks = strcmp(sfrom, fs);
fs(ks) = {sto};
% ks = find(strcmp(sfrom, fs));
% for i = 1:length(ks)
%    fs{ks(i)} = sto;
% end
end
%MATLAB precedence).
function res = expandeq(fs, fsorig,typesorig)
% if nargin == 1
%    [fs,types]=parseeq(equation,true,{'(',')'});
% elseif nargin == 2
%    if ischar(equation)
%       dotest = types;
%       equation=strrep(equation,char(183),' ');
%       [fs,types]=parseeq(equation,true,{'(',')'});
%    else
%       fs = equation;
%       equation=sprintf('%s',fs{:});
%    end
% end
vr = parsed_equation.vr;

%replace : with colon (no brackets supported only values/variables a:b
%or a:b:c)
f = strcmp(fs, ':');
if any(f)
    fis=find(strcmp(fs, '='), 1);
    if ~isempty(fis)
        f(1:fis) = false;
    end
    isnumvar = zeros(size(fs));
    isnumvar(f) = 2;
    for j = 1:length(fs)
        if isletter(fs{j}(1))||(fs{j}(1)=='_')||~isnan(str2double(fs(j)))
            isnumvar(j) =  1;
        end
    end
    f1 = strfind(isnumvar, [1 2 1 2 1]);
    for j = length(f1):-1:1
        fs=[fs(1:f1(j)-1) {'colon','(',fs{f1(j)},',',fs{f1(j)+2},',',fs{f1(j)+4},')'} fs(f1(j)+5:end) ];
        isnumvar = [isnumvar(1:f1(j) - 1) [0, 0, 0, 0, 0, 0, 0, 0] isnumvar(f1(j) + 5:end) ];
    end
    f1 = strfind(isnumvar, [1 2 1]);
    for j = length(f1):-1:1
        fs=[fs(1:f1(j)-1) {'colon','(',fs{f1(j)},',',fs{f1(j)+2},')'} fs(f1(j)+3:end) ];
        isnumvar = [isnumvar(1:f1(j) - 1) [0, 0, 0, 0, 0, 0] isnumvar(f1(j) + 3:end) ];
    end
    f1 = strfind(isnumvar, [0 2 0]);
    for j = length(f1):-1:1
        fs=[fs(1:f1(j)) {'colon','(','1',',','end',')'} fs(f1(j)+2:end) ];
        isnumvar = [isnumvar(1:f1(j)) [0, 0, 0, 0, 0, 0] isnumvar(f1(j) + 2:end) ];
    end
end
% replace [ with horzcat and vertcat
f = strcmp(fs, '[');
if any(f)
    f = find(f);
    f2 = find(strcmp(fs, ']'));
    if length(f)==1&&length(f2)==1
        fsnew = makecat(fs(f(1):f2(1)));
        fs = [fs(1:f(1) - 1) fsnew fs(f2(1) + 1:end)];
    end
end

types = updatetypes(fs,vr);
isoper1=(types==vr.oper) | (types==vr.funcname) | (types==vr.unary_oper_pre) | (types==vr.unary_oper_post);
isoper=(typesorig==vr.oper) | (typesorig==vr.funcname) | (typesorig==vr.unary_oper_pre) | (typesorig==vr.unary_oper_post);
if sum(isoper1) == sum(isoper)
    counts1 = ones(size(fsorig));
    for i = 2:length(fsorig)
        counts1(i) = counts1(i - 1) + length(fsorig{i - 1});
    end
    counts = zeros(size(fs));
    counts(isoper1) = counts1(isoper);
else
    counts = ones(size(fs));
end


l_startfunc = 9991;
l_endfunc = 9992;
%l_semicolon=13;
s = [];
if ~isempty(types)&&types(1)==vr.funcname
    %command syntax of functions - replace with function syntax
    fspace = find(strcmp(fs, ' '));
    if ~isempty(fspace)&&fspace(1)==2
        fs{2} = '(';
        if strcmp(fs{end}, '; ')
            fs{end} = ')';
        else
            fs{end + 1} = ')';
        end
        fs=cellrep(fs,' ',',');
        for j = 2:length(fs)
            if ~strncmp(fs{j},'''',1)&&~any(strcmp(fs{j},{',','(',')'}))
                fs{j}=sprintf('''%s''',fs{j});
            end
        end
    end
end
types(strncmp(fs, '@',1)) = vr.function_handle;
types(strncmp(fs, '''',1)&~strcmp(fs, '''')) = vr.string;
types(strcmp(fs, '(')|strcmp(fs, '«')) = vr.brack1;
types(strcmp(fs, ')')|strcmp(fs, '»')) = vr.brack2;
types(strcmp(fs, '=')) = vr.assign;
types(strcmp(fs,',')) = vr.comma;
f=find(strcmp(fs,';')|strncmp(fs,'%',1));
if ~isempty(f)
    fs = fs(1:f(1) - 1);
    types = types(1:f(1) - 1);
    counts = counts(1:f(1) - 1);
end
leftside = 'ans';
if any(types == vr.assign)
    f=find(types == vr.assign);
    f1=find(types(1:f) == vr.variable, 1);
    if ~isempty(f1)
        leftside = sprintf('%s', fs{f1});
        if fs{f1+1}==''''
            %keep the differentiate sign
            leftside=[leftside ''''];
        end
        k = 1;
        if f1+k<=length(types)&&types(f1+k)==vr.field%
            leftside  = sprintf('%s%s', leftside, fs{f1 + k});
            %   k = k + 1;
        end
    end
    fs = fs(f + 1:end);
    types = types(f + 1:end);
    counts =  counts(f + 1:end);
end
%remove extra brackets around variables
fv = 1;
while ~isempty(fv)
    fv = [strfind(types, [vr.brack1 vr.variable vr.brack2]) strfind(types, [vr.brack1 vr.number vr.brack2]) ...
        strfind(types, [vr.brack1 vr.string vr.brack2]) strfind(types, [vr.brack1 vr.function_handle vr.brack2])];
    fv = sort(fv);
    for i = length(fv):-1:1
        if fv(i)==1||(types(fv(i) - 1)~=vr.funcname)
            ndx = true(size(types));
            ndx(fv(i)) = false;
            ndx(fv(i) + 2) = false;
            types = types(ndx);
            counts =  counts(ndx);
            fs = fs(ndx);
        else
            fv(i) = [];
        end
    end
end
%remove extra brackets (add brackets adds brackets round a function
%(sin(x)) remove the inside brackets
f=find(types == vr.funcname);
jstart = f + 1;
jend = f;
for i = 1:length(f)
    if f(i) < length(types)
        %   jstart(i) = f(i) + 1;
        pop = 1;
        if types(jstart(i)) == vr.brack1
            jend(i) = jstart(i) + 1;
            while (jend(i)<=length(types))&&(pop > 1||types(jend(i))~=vr.brack2)
                if types(jend(i)) == vr.brack1
                    pop = pop + 1;
                elseif types(jend(i)) == vr.brack2
                    pop = pop - 1;
                end
                jend(i) = jend(i) + 1;
            end
            types(jstart(i)) = l_startfunc;
            types(jend(i)) = l_endfunc;
        else
            types(f) = vr.variable;
        end
    end
end
%add brackets in function arguments (is not done in parseeq)
f1=find(types == vr.comma);
for i = length(f1):-1:1
    fs=[fs(1:f1(i)-1) {')',',','('} fs(f1(i)+1:end)];
    types = [types(1:f1(i) - 1) [vr.brack2 vr.comma vr.brack1] types(f1(i) + 1:end)];
    counts = [counts(1:f1(i) - 1) [0 0 0] counts(f1(i) + 1:end)];
end
f1=find(types == l_startfunc);
for i = length(f1):-1:1
    fs=[fs(1:f1(i)-1) {'(','('} fs(f1(i)+1:end)];
    types = [types(1:f1(i) - 1) [vr.comma vr.brack1] types(f1(i) + 1:end)];
    counts = [counts(1:f1(i) - 1) [0 0 ] counts(f1(i) + 1:end)];
end
f1=find(types == l_endfunc);
for i = length(f1):-1:1
    fs=[fs(1:f1(i)-1) {')',')'} fs(f1(i)+1:end)];
    types = [types(1:f1(i) - 1) [vr.brack2 vr.comma] types(f1(i) + 1:end)];
    counts = [counts(1:f1(i) - 1) [0 0] counts(f1(i) + 1:end)];
end

br2=cumsum(types == vr.brack2);
poplevel=cumsum(types == vr.brack1) - [0 br2(1:end - 1)];

eqnrs = zeros(size(poplevel));
diffpop = [diff(poplevel) 0 - 1];
br2 = find(diffpop < 0);
isoper=(types==vr.oper) | (types==vr.funcname) | (types==vr.unary_oper_pre) | (types==vr.unary_oper_post);
eqnr = 1;
if isempty(poplevel)
    res = [];
    return;
end
for i = 1:length(br2)
    for j = diffpop(br2(i)):-1
        pop = poplevel(br2(i)) + 1 + j;
        ndx = 1:br2(i);
        ndx=poplevel(ndx)==pop & eqnrs(ndx)==0;
        if any(isoper(ndx))
            eqnrs(ndx) = eqnr;
            eqnr = eqnr + 1;
        else
            f = find(ndx);
            if ~isempty(f)
                poplevel(f(1):f(end)) = poplevel(f(1):f(end)) - 1;
            end
        end
    end
end
[sub.eqnr, ndx] = unique(eqnrs);
sub.level = poplevel(ndx);
[~,ndx] = sortrows([sub.level',sub.eqnr'],[ - 1 2]);
sub.level = sub.level(ndx);
sub.eqnr = sub.eqnr(ndx);


isvar=(types==vr.number) | (types==vr.string)| (types==vr.function_handle)| (types==vr.variable)| (types==vr.constant);

sub.operndx = zeros(size(sub.eqnr));
sub.args = cell(size(sub.eqnr));
sub.leftvar = cell(size(sub.eqnr));
n = 0;
for i = 1:length(sub.eqnr)
    f=find(eqnrs==sub.eqnr(i) & isoper);
    if ~isempty(f)
        n = n + 1;
        sub.charnr(n) = counts(f);
        sub.operndx(n) = f;
        sub.leftvar{n} = sprintf('#%03d', i);
        switch types(f)
            case vr.unary_oper_pre
                sub.args{n} = {findarg(f, 1, types, eqnrs,sub.eqnr, fs, isvar, isoper, poplevel)};
            case vr.unary_oper_post
                sub.args{n} = {findarg(f,  - 1,types,eqnrs,sub.eqnr,fs,isvar, isoper, poplevel)};
            case vr.oper
                sub.args{n} = {findarg(f,  - 1,types,eqnrs,sub.eqnr,fs,isvar, isoper, poplevel),findarg(f,1,types,eqnrs,sub.eqnr,fs,isvar, isoper, poplevel)};
            case vr.funcname
                sub.args{n} = funargs(f, types, eqnrs, sub.eqnr,fs, isvar, isoper, poplevel);
        end
    end
end
if n == 0
    s.args = fs;
    s.leftvar = leftside;
    s.oper = '';
    s.opertype = [];
    s.charnr = 1;
else
    sub.leftvar{n} = leftside;
    s=struct('leftvar','','args',{},'oper','','opertype',[],'charnr',0);
    for i = n:-1:1
        if sub.operndx(i) > 0
            s1.leftvar = sub.leftvar{i};
            s1.args = sub.args{i};
            s1.oper = fs{sub.operndx(i)};
            s1.opertype = types(sub.operndx(i));
            s1.charnr = sub.charnr(i);
            s(i) = s1;
        end
    end
end
%s.fs=fs;
%s.types=types;
%disp(equation);
%disp(sexpanded);
if nargout > 0
    res = s;
end


end

function error1 = dotest(orig_eq, new_eq, avar)
if nargin < 3
    avar = [];
end
error1 = true;
try
    strm = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(strm);
    res1 = evalin('base', orig_eq);
catch
    disp('Cannot test because I cannot run the original equation');
    return;
end
try
    strm = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(strm);
    if isempty(avar)
        res2 = evalin('base', new_eq);
    else
        evalin('base', new_eq);
        res2  = evalin('base', avar);
    end
catch err
    disp('error: test failed cannot run resulting equation');
    rethrow(err)
end
if (all(isnan(res1(:)))&&all(isnan(res2(:))))||(all(isinf(res1(:)))&&all(isinf(res2(:))))
    diff1 = 0;
elseif ~(size(res2,1)==size(res1,1)&&size(res2,2)==size(res1,2))
    diff1=inf;
else
    diff1 = sum(sum((res2 - res1).^2));
end
if diff1 < 1E-10
    disp('Success: the resulting equation gives the same result as the original');
    disp(diff1)
    error1 = false;
else
    disp(res1)
    disp(res2);
    disp(new_eq);
    error('grind:expandeq','Error in the new expanded equation with brackets: the result is different');
end
end
function fsnew = makecat(fs)
fsnew = fs;

if strcmp(fs{1},'[')&&strcmp(fs(end),']')
    iscomma=strcmp(fs,' ')|strcmp(fs,',');
    issemi = strcmp(fs, '; ');
    isbrack1=strcmp(fs,'[')|strcmp(fs,'(')|strcmp(fs,'«');
    isbrack2=strcmp(fs,']')|strcmp(fs,'}')|strcmp(fs,'»');%'«','»'
    poplevel = 0;
    for i = 1:length(fs)
        if isbrack1(i)
            poplevel = poplevel + 1;
        elseif isbrack2(i)
            poplevel = poplevel - 1;
        end
        if iscomma(i)&&poplevel > 1
            iscomma(i) = false;
        end
        if issemi(i)&&poplevel > 1
            iscomma(i) = false;
        end
        nrow = sum(issemi) + 1;
        ncol = round((sum(iscomma | issemi)+1) / nrow);
    end
    fc = cell(1, nrow * ncol);
    iscomma(1) = true;
    iscomma(end) = true;
    fcom = find(iscomma | issemi);
    for j = 2:length(fcom)
        fc(j - 1) = {[fs(max(1,fcom(j - 1) + 1):min(length(fs),fcom(j) - 1)),{','}]};
    end
    fc = transpose(reshape(fc, ncol, nrow));
    for j = 1:nrow
        fc(j,1)={[{'«','horzcat','('},fc{j,1}]};
        fend = fc{j, ncol};
        fend{end} = ')';
        fend{end + 1} = '»';
        if nrow > 1&&j < nrow
            fend{end + 1} = ',';
        end
        fc(j, ncol) = {fend};
    end
    fc = transpose(fc);
    fsnew = horzcat(fc{:});
    fsnew = fsnew(1:end);
    if nrow > 1
        fsnew=[{'«','vertcat','('} fsnew,{')','»'}];
    end
end
fspace = strcmp(fsnew, ' ');
fsnew(fspace) = {','};
end


%get function arguments (first and last brackets are also commas)
function args = funargs(f, types, eqnrs,eqnames,fs, isvar, isoper, poplevel)
nr = eqnrs(f);
iscomma=types == parsed_equation.vr.comma;
commandx=find(iscomma & eqnrs==nr);
args = cell(length(commandx) - 1, 1);
for j = 1:length(commandx) - 1
    eqnrs1 = eqnrs(commandx(j) + 1:commandx(j + 1) - 1);
    oper1 = isoper(commandx(j) + 1:commandx(j + 1) - 1);
    if ~any(oper1)
        var1 = isvar(commandx(j) + 1:commandx(j + 1) - 1);
        f1 = find(var1, 1);
        if ~isempty(f1)
            args{j} = fs{commandx(j) + f1};
            k = 1;
            while (f1+commandx(j)+k<=length(fs))&&types(commandx(j)+f1+k)==parsed_equation.vr.field
                args{j} = sprintf('%s%s', args{j}, fs{commandx(j) + f1+k});
                k = k + 1;
            end
        end
    else
        oper1 = isoper(commandx(j) + 1:commandx(j + 1) - 1);
        pop1 = poplevel(commandx(j) + 1:commandx(j + 1) - 1);
        pop2 = pop1(oper1);
        eq2 = eqnrs1(oper1);
        [~, imin] = min(pop2);
        args{j} = sprintf('#%03d', find(eqnames == eq2(imin)));
    end
end
end
function s = findarg(f, direc, types, eqnrs, eqnames, fs, isvar, ~, ~)
nr = eqnrs(f);
ndx = false(size(types));
j = 1;
while (f+direc * j < length(ndx))&&(f+direc * j > 0)&&((eqnrs(f+direc * j)==0))
    j = j + 1;
end

if direc == 1
    ndx(f + j:end) = true;
else
    ndx(1:f - j) = true;
end
ndx=ndx & eqnrs==nr;
f1 = find(isvar & ndx);
if ~isempty(f1)
    s = fs{f1};
    k = 1;
    while f1(1)+k<=length(fs)&&types(f1(1)+k)==parsed_equation.vr.field
        s = sprintf('%s%s', s, fs{f1 + k});
        k = k + 1;
    end
else
    while (f +direc *  j > 1)&&f +direc *  j < length(eqnrs)&& ...
            eqnrs(f +direc *  j)==nr&&any(types(f +direc *  j)==[parsed_equation.vr.brack1, parsed_equation.vr.brack2])
        j = j + 1;
    end
    s = sprintf('#%03d', find(eqnames == eqnrs(f  + direc *  j)));
end

end

function result = mupad2matlab(mupadstring, res)
    function equation = renamevars(equation, oldvars, newvars)
        fs = parsed_equation(equation).fs;
        for k = 1:length(oldvars)
            fs(strcmp(fs, strtrim(oldvars{k}))) = newvars(k);
        end
        
        equation = sprintf('%s', fs{:});
    end
if iscell(mupadstring)
    result=cell(size(mupadstring));
    for i=1:length(mupadstring)
        result{i}=mupad2matlab(mupadstring{i},res);
    end
    return;
end
s=strtrim(char(mupadstring));
s=regexprep(s,'[*]1.0\>','');
if strncmp(s, 't0 = ', 5)
    s = s(6:end);
end
if any(strcmp(res.renamewords, 'iif'))
    cc = strtrim(str2cell(sprintf(s)));
    xbrack = 0;
    for m = 1:length(cc)
        s1 = cc{m};
        if strncmp(s1, 't0 =', 4)
            s1 = s1(5:end);
        elseif strncmp(s1, 'if ', 3)
            s1 = ['iif(' s1(4:end)];
            if s1(end) ~= ','
                s1 = [s1 ','];
            end
            
        elseif strncmp(s1, 'elseif ', 7)
            xbrack = xbrack + 1;
            s1 = [',iif(' s1(8:end)];
            if s1(end) ~= ','
                s1 = [s1 ','];
            end
            
        elseif strncmp(s1, 'else', 4)
            s1 = [',' s1(5:end)];
        elseif strncmp(s1, 'end', 3)
            s1 = ')';
        end
        
        cc{m} = s1;
    end
    
    while xbrack>0
        cc{end} = [cc{end} ')'];
        xbrack=xbrack-1;
    end
    
    s = sprintf('%s', cc{:});
    s(s==';')='';
end

s=strrep(s,'\n','');
s=strrep(s,'&&','&');
s=strrep(s,'||','|');
if ~isempty(res.rename2words)
    s = renamevars(s, res.rename2words, res.renamewords);
end

s(s=='#')='''';
if s(end)==';'
    s=s(1:end-1);
end

s=s(s ~= ' ');
result= s;

end
