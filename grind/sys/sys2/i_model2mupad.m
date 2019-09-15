function syms1=i_model2mupad(dosensitiv,syms,maxtime,runwhich)
%this function translates the model to mupad and generates all analytical
%results that are needed in various functions
%function handles are created (including the list of parameters)
global g_grind;
%Help from symbolic toolbox
if nargin==0
    dosensitiv=true;
end
if nargin<2
    if ~isfield(g_grind,'syms')
        syms=[];
    else
        syms=g_grind.syms;
    end
end
if nargin<3
    maxtime=30;
end
if nargin<4
    runwhich=ones(7,1);
end
if numel(g_grind.pars)==0
    dosensitiv=false;
    dojacp=false;
else
    dojacp=true;
end
if (~isempty(syms.Jacobian)&&~isempty(syms.Hessian)&&(~dojacp||~isempty(syms.Jacobianp)))&&(~dosensitiv||~isempty(syms.Sensitivp))
    if nargout==1
        syms1=syms;
    end
    return;
    %symbolics already done
end

i_parcheck;
if isempty(syms)
    %    syms.fulleq={};
    syms = struct('Jacobian',{{}},'Jacobianp',{{}},'Sensitivp',{{}},'Hessian',{{}},'Hessianp',{{}},'der3',{{}},'der4', {{}},'der5',{{}});
end
syms.errormsg = {};
if any(strcmp(g_grind.funcnames.names,'pi'))
    error('grind:symbolic','Assignments to pi not supported');
end

if any(strcmp(g_grind.model,'%external odefile'))
    g_grind.syms.errormsg={'cannot derive symbolic differentials'};
    syms1=g_grind.syms;
    return;
end

model=regexprep(g_grind.model,'[%].*','');
model=regexprep(model,'\((t|t-1|t\+1)\)',''); %remove (t-1) (t) (t+1)


[allfunct,types]=i_localfunctions(model);
ndx=strcmp(allfunct,'ln');
types(ndx)={'built-in'};
nondiffs=[allfunct(strcmp(types,'local')) intersect(allfunct,{'if','for','Dx','Dxx','Dy','Dyy','floor','ceil','lag','neighborcells','sum'})];
if ~isempty(nondiffs)
    s = sprintf('%s, ',nondiffs{:});
    syms.errormsg{end + 1} = sprintf('Cannot differentiate the following functions in the model:\n%s', s(1:end - 2));
    if nargout==1
        syms1=syms;
    end
    return;
end

res1 = openmupad(true);  %translates the model to mupad
warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
warning( 'off', 'symbolic:mupadmex:MuPADTextWarning'); %R2011

dim = g_grind.statevars.dim;
statevars = cell(parsed_equation(g_grind.statevars.names).mupad_syntax);
pars = cell(parsed_equation(g_grind.pars).mupad_syntax);

if isempty(syms.Jacobian)
    if ~isfield(syms,'fulleq')||isempty(syms.fulleq)
        %expanded equations (can be used in mupad)
        eqns = getmupadstrings('eqns#', false);
        syms.fulleq = fulljac(eqns.eq);
        %     eqs = sprintf('%s#,',statevars{:});
        %     eqs=regexp(eqs,',','split');
        %     eqs=eqs(1:end-1);
        %     eq= translateback(eqs, size(eqs), res1);
        %     syms.fulleq=eq.eq;
    end
    %Analytical Jacobian for state variables
    if g_grind.solver.haslags
        syms.errormsg{end + 1} = 'Jacobian for delay differential equation not yet supported';
    else
        if ~g_grind.solver.isimplicit
            tic
            %calculation of Jacobians (getJacobi# see grind.mu) is usually not the bottleneck
            %translating the higher derivatives back to MATLAB can be very slow
            %iif statements can lead to extremely long equations
            srw=sprintf('%d,',runwhich);
            evalin(symengine,sprintf('getJacobi#(%g,[%s])',maxtime*1000/2,srw(1:end-1)));
            jac = getmupadstrings('Jacobian#',  false);
            syms.Jacobian = fulljac(jac.eq);
            jac = getmupadstrings('Jacobianp#',  false);
            syms.Jacobianp = fulljac(jac.eq);
            if toc<maxtime
                hess = getmupadstrings('Hessian#', true);
                syms.Hessian=hess.eq;
            end
            if toc<maxtime
                hess = getmupadstrings('Hessianp#',false);
                syms.Hessianp=hess.eq;
            end
            if toc<maxtime
                der3 = getmupadstrings('der3#', true);
                syms.der3=der3.eq;
            end
            if toc<maxtime
                der4 = getmupadstrings('der4#', true);
                syms.der4=der4.eq;
            end
            if toc<maxtime
                der5 = getmupadstrings('der5#', true);
                syms.der5=der5.eq;
            end
            
        else
            s1 = sprintf('%s#,',statevars{:});
            evalin(symengine,sprintf('implicitJacs#([%s])',s1(1:end - 1)));
            jac = getmupadstrings('Jacobian#', false);
            syms.Jacobian = fulljac(jac.eq);
            jac = getmupadstrings('Jacobian_y#', false);
            syms.Jacobian_y = fulljac(jac.eq);
            jac = getmupadstrings('Jacobian_yp#', false);
            syms.Jacobian_yp = fulljac(jac.eq);
            %we can use the general formula vfor implicit differentiation of a
            %function R(x,y)=0
            %
            %dy/dx= -(dR/dx)/(dR/dy)
            %comes from the generalized chain rule
            %0 = R(x,y)
            %d 0/dx = d(R(x,y)/dx
            % 0= d(R(x,y)/dx
            %
            %Another method is to calculate the implicit differntial and solve unknowns:
            %R:=subs(R,y=y(x));
            %dRdx :=diff(R,x);
            %dydt:=solve(dRdx=0,diff(y(x),y))
            %same result as above
        end
        
    end
    [syms,hasdiff]=translateallback(syms,res1,dim);
    if hasdiff
        syms.errormsg='Could not find derivatives of all functions';
    end
end



%parameter sensitivities are a bit slow

if isempty(syms.Sensitivp)&&dosensitiv
    evalin(symengine,'delete(g_res)');
    origpars=[g_grind.pars(:);g_grind.statevars.names(:)];
    pars=[pars;statevars];
    symvars = cell(length(statevars), length(pars));
    for m = 1:length(pars)
        apar = pars{m};
        ipar = i_getno(origpars{m});
        if ipar.isvar
            apar = [apar '##'];
        end
        
        for k = 1:length(statevars)
            avar = statevars{k};
            s = '';
            for i = 1:dim
                s=sprintf('%s%s=%s(%s),',s,statevars{i},statevars{i},apar);
            end
            
            evalin(symengine,sprintf('g_res1:=subs(eqns#[%d],%s);',k,s(1:end - 1)));
            evalin(symengine,sprintf('g_res1:=diff(g_res1,%s)',apar));
            s = '';
            % sensname=sprintf('g_d%s_%s',avar,apar);
            % statevarname=statevars{i};
            for i = 1:g_grind.statevars.dim
                statevarname = sprintf('g_X1[%d]', i);
                sensname = sprintf('g_X1[%d]', g_grind.statevars.dim * m + i);
                s=sprintf('%sdiff(%s(%s),%s)=%s,',s,statevars{i},apar,apar,sensname);
                s=sprintf('%sD(%s)(%s)=%s,',s,statevars{i},apar,sensname);
                s=sprintf('%s%s(%s)=%s,',s,statevars{i},apar,statevarname);
            end
            
            symvars{k, m} = sprintf('d%s/d%s', avar, apar);
            evalin(symengine,sprintf('g_res[%d,%d]:=subs(g_res1,%s);',k,m,s(1:end - 1)));
        end
        
    end
    
    res =  translateback('g_res', size(symvars), res1);
    res.symvars = symvars;
    for i=1: g_grind.statevars.dim*(length(pars)+1)
        res.eq=regexprep(res.eq,sprintf('(?<![a-zA-Z_0-9])g_X1[(]%d[)]',i),sprintf('g_X1(%d,:)',i));
    end
    
    %reset(symengine);
    syms.Sensitivp = res.eq;
    if res.hasdiff
        s = sprintf('%s, ',res.diff_funcs{:});
        g_sens.symbolic = false;
        syms.errormsg{end + 1} = sprintf('Cannot differentiate the following functions in the model:\n%s', s(1:end - 2));
    end
    
end

warning('on','symbolic:generate:FunctionNotVerifiedToBeValid');
warning( 'on', 'symbolic:mupadmex:MuPADTextWarning'); %old version
if nargout==0
    g_grind.syms=syms;
else
    syms1=syms;
end

function jac = replacestatevars(jac, dim, isimplicit)
if isimplicit
 for i = 1:dim
    %replace statevar symbols with g_X1(i)
    jac=regexprep(jac,sprintf('\\<%s#\\>',i_statevars_names(i)),sprintf('g_X2(%d,:)', i));
 end
end
for i = 1:dim
    %replace statevar symbols with g_X1(i)
    jac=regexprep(jac,sprintf('\\<%s\\>',i_statevars_names(i)),sprintf('g_X1(%d,:)', i));
end

%jac=regexprep(jac,'\<([0-9])[.][0]\>','$1');

function jac=fulljac(sjac)
if isempty(sjac)
    jac={};
    return;
end
jac=cell(sjac.size);
jac(:)={'0'};
if isfield(sjac,'unique')
    for i=1:length(sjac.unique)
        jac(sjac.unique(i).indices)={sjac.unique(i).equation};
    end
end

%TRANSLATE THE MODEL TO MUPAD
function result = openmupad(doevaluate)
global g_grind;
%start of translation to mupad
if ~i_hastoolbox('symbolic')
    error('grind:mupad','Symbolic toolbox is required for this option');
end

if g_grind.statevars.vector
    error('grind:mupad','Vector notation not supported for symbolic Jacobians');
end

if doevaluate
    reset(symengine);
    %  evalin(symengine, 'g_running:=1'); %flag to be able to check if there is a current session running
end

themodel = g_grind.model(g_grind.modelndx);
themodel=regexprep(themodel,'[%].*',''); %remove comments
if g_grind.solver.isdiffer
    %remove brackets after = sign
    themodel=regexprep(themodel,'(?<=[=].*)\((t|t-1|t\+1)\)',''); %may be slow
end
themodel=strrep(themodel,';;',';');



if g_grind.solver.isimplicit
    k=1;
    for i = 1:length(themodel)
        s = strtrim(themodel{i});
        s(s=='''')='#';
        if ~isempty(s)&&s(1) == '0'
            s1= regexp(s,'(?<=[\[])[^\]]*','match');
            if isempty(s1)
                s(s==';')=[];
                s={s(4:end)};
            else
                s=regexp(s1{1},'[;]','split');
            end
            for j=1:length(s)
                themodel{k}=sprintf('%s'' = %s;',g_grind.statevars.names{j},s{j});
                k=k+1;
            end
        else
            themodel{k} = s;
            k=k+1;
        end
    end
    
end

% if g_grind.solver.isdiffer
%    themodel = removebracks(themodel);
% end
statevars = cell(parsed_equation(g_grind.statevars.names).mupad_syntax);
pars = cell(parsed_equation(g_grind.pars).mupad_syntax);

[obj, result] = parsed_equation(themodel).mupad_syntax;
evalin(symengine,sprintf('eqns#:=array(1..%d);',length(statevars)));
s=sprintf('%s,',statevars{:});
evalin(symengine,sprintf('vars#:=[%s];',s(1:end-1)));
s=sprintf('%s,',pars{:});
evalin(symengine,sprintf('pars#:=[%s];',s(1:end-1)));
result.modeleq = obj;
if doevaluate
    k = 0;
    for i = 1:length(result.modeleq)
        feq=find(strcmp(result.modeleq(i).fs, '='));
        if ~isempty(feq)&&feq(1)==3
            seq = char(result.modeleq(i));
            %  disp(seq)
            evalin(symengine, seq);
            %assign aux vars
        elseif ~isempty(feq)
            seq = sprintf('%s', result.modeleq(i).fs{feq(1) + 1:end});
            k = k + 1;
            seq=sprintf('eqns#[%d]:=%s; ',k,seq);
            %    disp(seq)
            evalin(symengine, seq);
        end
        
    end
    oldcd=cd;
    cd(fullfile(grindpath,'sys2'))
    evalin(symengine,'read("grind.mu")');
    cd(oldcd);
end

function result = getmupadstrings(symbolics,doperm)
%other ways can be EXTREMELY slow
symbs=char(evalin(symengine,sprintf('expr2text(%s)',symbolics)));
if strcmp(symbs,'NIL')
    result.eq=[];
    return;
end
symbs=symbs(1:end-1);

eqns=regexp(symbs,'[,][ ][\(]([0-9\, ])*[\)][ ][=][ ]','split');
sindices=regexp(symbs,'(?<=[ ][\(])([0-9\, ])*(?=[\)][ ])','match');
if isempty(sindices)
    eqns=regexp(symbs,'[,][ ]([0-9])*[ ][=][ ]','split');
    sindices=regexp(symbs,'(?<=[,][ ])([0-9])*(?=[ ][=])','match');
end
%eval(['[' sprintf('%s; ',sindices{:}) ']']);
s=regexp(eqns{1},'(?<=[\.][\.])[0-9]*' ,'match');
eqns=eqns(2:end);
result.eq.size=str2num(sprintf('%s ',s{:})); %#ok<ST2NM>
subss=sscanf(sprintf('%s, ',sindices{:}),'%d,',[numel(result.eq.size),numel(sindices)])';
if length(result.eq.size)==1
    result.eq.size=[1 result.eq.size];
end
if ~doperm
    c=cell(1,size(subss,2));
    for i=1:size(subss,2)
        c{i}=subss(:,i);
    end
    indices=sub2ind(result.eq.size,c{:});
end
uni=unique(eqns);
k=1;
for i=1:length(uni)
    if ~any(strcmp({'','0'},uni{i}))
        if doperm
            ndx=strcmp(uni{i},eqns);
            subs1=subss(ndx,:);
            psubs=[];
            for j=1:size(subs1,1)
                us=unique(perms(subs1(j,2:end)),'rows');
                us=[zeros(size(us,1),1)+subs1(j,1) us]; %#ok<AGROW>
                psubs=[psubs; us]; %#ok<AGROW>
            end
            c=cell(1,size(psubs,2));
            for j=1:size(psubs,2)
                c{j}=psubs(:,j);
            end
            indices=sub2ind(result.eq.size,c{:});
            result.eq.unique(k)=struct('equation',uni{i},'indices',indices);
        else
            result.eq.unique(k)=struct('equation',uni{i},'indices',indices(strcmp(uni{i},eqns)));
        end
        k=k+1;
    end
end

function [syms,hasdiff] = translateallback(syms,res,dim)
%TRANSLATE The mupad result back (fast)

%extract all strings (eqns) for efficiency
fsyms=fieldnames(syms);
eqns={};
startndxs=zeros(size(fsyms));
k=1;
for i=1:length(fsyms)
    startndxs(i)=k;
    if iscell(syms.(fsyms{i}))
        n=numel(syms.(fsyms{i}));
        eqns=[eqns;syms.(fsyms{i})(:)];
    elseif isstruct(syms.(fsyms{i}))&&isfield(syms.(fsyms{i}), 'unique')
        n=numel(syms.(fsyms{i}).unique);
        eqns=[eqns; {syms.(fsyms{i}).unique(:).equation}'];
    else
        n=0;
    end
    k=k+n;
end
isimplicit=any(strcmp(fsyms,'Jacobian_yp'));
%replacements in vectorized regexp 

repstr={'(?<=[^:<>])[=]','==';...
    '<>','~=';...
    '[&][&]','&';...
    '[|][|]','|';...
    '\<not\>','~';...
    '\<or\>','|';...
    '\<and\>','&';...
    '\<Re\>','real';...
    '\<Im\>','imag';...
    '\<PI\>','pi';...
    '\<arcsin\>','asin';...
    '\<arccos\>','acos';...
    '\<arctan\>','atan'};
eqns=regexprep(eqns,repstr(:,1),repstr(:,2));
ndx= regexp(eqns,'\<piecewise\>', 'once');
ndx=find(~cellfun('isempty',ndx));

for j=1:length(ndx)
    %replace piecewise with iif (slow)
    s=eqns{ndx(j)};
    % if ~isempty(regexp(s,'\<piecewise\>', 'once'))
    eq=parsed_equation(s);
    f=find(strcmp(eq.fs,'piecewise'),1);
    while any(f)
        %get the arguments of piecewise (nested possible)
        [args,js]=eq.getfunctionpars(f);
        args=regexprep(args,'[\[\]]','');
        %uneven arguments are conditions, even are values (except
        %the last)
        conds=false(size(args));
        conds(1:2:end-1)=true;
        %if there are 4 arguments often the second condition is the
        %opposite of the first, for efficiency can be removed.
        if numel(args)==4&&~any(strcmp(eq.fs,'&'))&&~any(strcmp(eq.fs,'|'))
            args2=regexprep(args([1,3]),'~=','==');
            if (strcmp(args{1},['~' args{3}])||strcmp(args{3},['~' args{1}])||strcmp(args2{1},args2{2}))
                args(3)=[];
                conds(3)=[];
            end
        end
        if conds(end-1)
            args{end+1}='NaN';
            conds(end+1)=false;
        end
        for i=1:length(args)
            if conds(i)
                if strcontains(args{i},'~')
                    %AU! the precedence rules are different in
                    %mupad: not has lower precedence
                    mupadprec=parsed_equation.operinfo;
                    mupadprec.prec(mupadprec.prec>6)=mupadprec.prec(mupadprec.prec>6)+1;
                    mupadprec.prec(strcmp(mupadprec.opers,'~'))=7;
                    arg=parsed_equation(args{i});
                    fs=arg.addbrackets(false,mupadprec);
                    args{i}=sprintf('%s',fs{:});
                end
                args{i}=['iif(' args{i}];
            end
        end
        thefunct=sprintf('%s,',args{:});
        thefunct=sprintf('%s%s',thefunct(1:end-1),repmat(')',1,sum(conds)));
        s=sprintf('%s%s%s',sprintf('%s',eq.fs{1:js(1,1)-3}),thefunct,sprintf('%s',eq.fs{js(end,2)+2:end}));
        eq=parsed_equation(s);
        f=find(strcmp(eq.fs,'piecewise'),1);
    end
    eqns{ndx(j)}=char(eq);
end

if ~isempty(res.rename2words)
    eqns = renamevars(eqns, res.rename2words, res.renamewords);
end

eqns=regexprep(eqns,';$','');
%eqns= strrep(eqns,'#','''');
eqns= strrep(eqns,' ','');
difndx=regexp(eqns,'[^a-zA-Z_0-9#]diff(','once');
hasdiff=any(~cellfun('isempty',difndx));
if ~hasdiff
    difndx=regexp(eqns,'(?<![A-Za-z0-9])D\((\[[0-9]*\], )?','once');
    hasdiff=any(~cellfun('isempty',difndx));
end
%result.hasdiff =  result.hasdiff||~isempty(||~isempty(regexp(s,'(?<![A-Za-z0-9])D\((\[[0-9]*\], )?','once'));
eqns = replacestatevars(eqns, dim, isimplicit);

%put the strings back where they belong
for i=1:length(fsyms)
    if iscell(syms.(fsyms{i}))
        n=numel(syms.(fsyms{i}));
        syms.(fsyms{i})(:)=eqns(startndxs(i):startndxs(i)+n-1);
    elseif isstruct(syms.(fsyms{i}))&&isfield(syms.(fsyms{i}), 'unique')
        n=numel(syms.(fsyms{i}).unique);
        for k=1:n
            syms.(fsyms{i}).unique(k).equation=eqns{startndxs(i)+k-1};
        end
    end
end
 


function result = translateback(symbolics,siz, res, vfull)
%slow but good
%siz = size(symbolics);
if nargin<4
    vfull=true(siz);
end
result.eq = cell(siz);
result.hasdiff = false;
siz3=1;
siz4=1;
siz5=1;
siz6=1;
if length(siz)>=3
    siz3=siz(3);
end
if length(siz)>=4
    siz4=siz(4);
end
if length(siz)>=5
    siz5=siz(5);
end
if length(siz)>=6
    siz6=siz(6);
end
for i = 1:siz(1)
    for j = 1:siz(2)
        for k=1:siz3
            for l=1:siz4
                for m=1:siz5
                    for n=1:siz6
                        if ischar(symbolics)
                            if ~vfull(i,j,k,l,m,n)
                                s='0';
                            else
                                if length(siz)==2
                                    s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(%s[%d,%d])',symbolics, i, j))));
                                elseif length(siz)==3
                                    s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(%s[%d,%d,%d])',symbolics, i, j, k))));
                                elseif length(siz)==4
                                    s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(%s[%d,%d,%d,%d])',symbolics, i, j, k, l))));
                                elseif length(siz)==5
                                    s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(%s[%d,%d,%d,%d,%d])',symbolics, i, j, k, l,m))));
                                elseif length(siz)==6
                                    s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(%s[%d,%d,%d,%d,%d,%d])',symbolics, i, j, k, l,m, n))));
                                end
                            end
                        else
                            s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(%s)', char(symbolics(i, j, k))))));
                        end
                        
                        %      s = strtrim(char(evalin(symengine, sprintf('generate::MATLAB(g_res[%d,%d])', i, j))));
                        if vfull(i,j,k,l,m,n)
                            if any(strncmp(s, {'t0 = ','t1 = '}, 5))
                                s = s(6:end);
                            end
                            s=regexprep(s,'[*]1.0\>','');
                            if any(strcmp(res.renamewords, 'iif'))
                                cc = strtrim(str2cell(sprintf(s)));
                                xbrack = 0;
                                for m1 = 1:length(cc)
                                    s1 = cc{m1};
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
                                    
                                    cc{m1} = s1;
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
                            
                            %       [f1,f2]=regexp(s,'[(][0-9[.]eE]*[)]','start','end');
                            %       for k = length(f1):-1:1
                            %          d = str2double(s(f1(k) + 1:f2(k) - 1));
                            %          s = [s(1:f1(k)) sprintf('%d', d) s(f2(k):end)];
                            %       end
                            
                            s(s=='#')='''';
                            if s(end)==';'
                                s=s(1:end-1);
                            end
                            
                            s=s(s ~= ' ');
                            result.hasdiff =  result.hasdiff||~isempty(regexp(s,'[^a-zA-Z_0-9#]diff(','once'))||~isempty(regexp(s,'(?<![A-Za-z0-9])D\((\[[0-9]*\], )?','once'));
                        end
                        result.eq{i, j, k, l, m, n} = s;
                    end
                end
            end
        end
    end
    
end

if result.hasdiff
    s = sprintf('%s, ',result.eq{:});
    fdiff = regexp(s, '[^a-zA-Z_0-9#]diff(');
    fD = regexp(s, '[^a-zA-Z_0-9#]D(');
    result.diff_funcs = cell(length(fdiff) + length(fD), 1);
    fcomma=s == ',';
    for i = 1:length(fdiff)
        f1 = find(fcomma(fdiff(i):end), 1) + fdiff(i) - 1;
        result.diff_funcs{i} = s(fdiff(i) + 5:f1 - 1);
    end
    fbrack=s == ')';
    for i = 1:length(fD)
        f1 = find(fbrack(fD(i):end), 1) + fD(i) - 1;
        result.diff_funcs{i} = s(fD(i) + 3:f1 - 1);
        f=find(result.diff_funcs{i} == ',',1);
        if ~isempty(f)
            result.diff_funcs{i} = result.diff_funcs{i}(f + 1:end);
        end
        
    end
    
    result.diff_funcs = unique(result.diff_funcs);
end


% function res=makesparse(c)
% res.size=size(c);
% uni=unique(c(:));
% k=1;
% for i=1:length(uni)
%     if ~any(strcmp({'','0'},uni{i}))
%         res.unique(k)=struct('equation',uni{i},'indices',find(strcmp(uni{i},c(:))));
%         k=k+1;
%     end
% end

function equation = renamevars(equation, oldvars, newvars)

%fs = parsed_equation(equation).fs;
for k = 1:length(oldvars)
    equation= regexprep(equation,sprintf('\\<%s\\>',oldvars{k}),newvars{k});
    % fs(strcmp(fs, strtrim(oldvars{k}))) = newvars(k);
end

