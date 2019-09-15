function [ggrind, parsedmodel, parsedspec] = i_analysemodel(amodel)
%analyse the model as typed by the user to extract all parameters
%state variables, auxiliaries etc.
%The auxiliary variables are sorted based on dependency if necessary
%The odefile is constucted and the g_grind stucture filled.
%first parse is rather slow, the subsequent analysis is very fast and efficient.

%constants
vr = parsed_equation.vr;
keywords={'t', 'pi',  'inf','Inf', 'nan','NaN','true','false','elseif','else','if','end','while','for',...
    'switch','otherwise','case','global','function','return','definepars','defextern','defpermanent',...
    'lag', 'externlag','dwiener', 'djump',   'rednoise', 'definespace', 'rand',     'randn',      'randlogn',  'randi',   'randperm',...
    'boxcartrain','setevent','implicitdisperse','implicitvars',':', '(', ')',',', ' ',  ';','{','}','setodefile'};
typkw = [vr.t, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.else, vr.if, vr.block2, vr.block1, vr.block1, ...
    vr.block1,   vr.kw,    vr.kw,vr.kw, vr.kw,     vr.kw,vr.definepars,vr.defextern,vr.defpermanent, ...
    vr.lag,vr.externlag,vr.dwiener,vr.djump, vr.rednoise,vr.definespace, vr.stochfun,vr.stochfun,vr.stochfun,vr.stochfun, vr.stochfun, ...
    vr.boxcartrain, vr.setevent, vr.implicitdisperse, vr.implicitvars,vr.colon, vr.brack1, vr.brack2, vr.comma, vr.space, vr.semicolon,vr.curlbrack1,vr.curlbrack2, vr.setodefile];
%remove dots... in g_grind.model
%maybe not elegant, but it is better not to support them (as it makes
%things complicated)
f1=strcontains(amodel,'%');
fddd=find(~cellfun('isempty',regexp(amodel,'[.][.][.]$'))&~f1);
%fddd=find(~cellfun(@isempty,regexp(amodel,'(?<!%.*)[.][.][.]$'))); VERY
%slow
for i=length(fddd):-1:1
    if fddd(i)<length(amodel)
        s=amodel{fddd(i)};
        s=[s(1:end-3) amodel{fddd(i)+1}];
        amodel{fddd(i)}=s;
        amodel(fddd(i)+1)=[];
    end
end
ggrind=i_init_g_grind;
ggrind.model = transpose(amodel);
ggrind.locfunc = {};
ggrind.locfuncname = {};
ggrind.theodefile = {};
aparsedmodel = cell(size(amodel));
aparsedspec = cell(size(amodel));
maxlen = 0;
infunction = 0;
ggrind.errors = {};
k = 1;
k1 = 0;

%parse and remove local functions
objs = parsed_equation(amodel);

for i = 1:length(amodel)
    %parsing, slowest step:
    obj = objs(i);
    aparsedmodel{i} = obj.fs;
    aparsedspec{i} = obj.types;
    %   [aparsedmodel{i}, aparsedspec{i}] = parseeq(strtrim(amodel{i}));
    if ~isempty(aparsedmodel{i})&&strcmp(aparsedmodel{i}{1}, 'function')
        infunction = true;
        f=find(strcmp(aparsedmodel{i}, '='));
        k1 = k1 + 1;
        if ~isempty(f)
            ggrind.locfuncname{k1, 1} = aparsedmodel{i}{f + 1};
        else
            ggrind.locfuncname{k1, 1} = '';
        end
        
    end
    
    if ~infunction
        if length(aparsedspec{i}) > maxlen
            maxlen = length(aparsedspec{i});
        end
        
    else
        ggrind.locfunc{k, 1} = amodel{i};
        k = k + 1;
        if strcmp(aparsedmodel{i}{1}, 'return')
            infunction = false;
        end
        
        aparsedmodel{i} = {''};
        aparsedspec{i} = vr.empty;
    end
    
end

%merge '...' and remove empty lines
% for i = length(aparsedmodel):-1:1
%    if (length(aparsedmodel{i})==1 && isempty(aparsedmodel{i}{1})) ||isempty(aparsedmodel{i})
%       aparsedmodel(i) = [];
%       aparsedspec(i) = [];
%    elseif i < length(aparsedmodel) && strcmp('...', aparsedmodel{i}{end})
%       aparsedmodel{i} = [aparsedmodel{i}(1:end - 1) aparsedmodel{i + 1}];
%       aparsedspec{i} = [aparsedspec{i}(1:end - 1) aparsedspec{i + 1}];
%       if length(aparsedspec{i}) > maxlen
%          maxlen = length(aparsedspec{i});
%       end

%       aparsedmodel(i + 1) = [];
%       aparsedspec(i + 1) = [];
%    end

% end

%make a single matrix of the model (for efficiency)
parsedmodel = cell(length(aparsedmodel), maxlen);
parsedmodel(:) = {''};
parsedspec = zeros(size(parsedmodel)) + vr.empty;
for i = 1:length(aparsedmodel)
    if ~isempty(aparsedmodel{i})
        parsedmodel(i, 1:length(aparsedmodel{i})) = aparsedmodel{i};
        parsedspec(i, 1:length(aparsedmodel{i})) = aparsedspec{i};
    end
    
end

[index] = makeindex(parsedmodel, parsedspec, vr);
[parsedspec,parsedmodel,index,beforeeq]=updatespecs(parsedspec,parsedmodel,index,vr,keywords,typkw);

%if definepars is used in command style make the arguments variables
[fi, fj] = findstrcmpndx(parsedmodel, 'definepars', index);
[fi2, fj2] = findstrcmpndx(parsedmodel, 'global', index);
if ~isempty(fi2)
    warning('grind:model','Do not use global in your GRIND model, it is ignored');
    parsedmodel(fi2,fj2)={'%global'};
end

[fi3, fj3] = findstrcmpndx(parsedmodel, vr.implicitvars, index);

fi = [fi; fi2; fi3];
fj = [fj; fj2; fj3];

for i = 1:length(fi)
    for j = fj(i) + 1:size(parsedmodel, 2)
        if ~isempty(strtrim(parsedmodel{fi(i),j}))&&~any(strcmp(parsedmodel{fi(i),j},{'(',',',')'}))
            parsedspec(fi(i), j) = vr.variable;
            v = parsedmodel{fi(i), j};
            f = strfind(v,',');
            if ~isempty(f)
                f = [0 f length(v) + 1];
                for ii = 2:length(f)
                    parsedmodel{fi(i), j + ii - 2} = v(f(ii - 1) + 1:f(ii) - 1);
                    parsedspec(fi(i), j + ii - 2) = vr.variable;
                end
                
                index = makeindex(parsedmodel, parsedspec, vr);
                index.fbeforeeq = beforeeq(index.fndx);
                break;
            end
            
        end
        
    end
    
end

%replace 'dN','/','dt' by  '','N','''
fdt = findstrcmpndx(parsedmodel, '/', index, index.fbeforeeq);
if ~isempty(fdt)
    [is, js] = ind2sub(size(parsedspec), fdt);
    % if is>1 && strcmp(parsedspec(is-1,js),'/')
    changed = false;
    for i = 1:length(is)
        if js(i) > 1&&js(i) < size(parsedspec,2)&&strcmp(parsedmodel(is(i), js(i) + 1), 'dt')
            var = parsedmodel{is(i), js(i) - 1};
            var = var(2:end);
            parsedmodel(is(i),js(i)-1:js(i)+1)={'',var,''''};
            parsedspec(is(i), js(i) - 1:js(i) + 1) = [vr.empty, vr.variable, 0];
            changed = true;
        end
        
    end
    
    if changed
        index = makeindex(parsedmodel, parsedspec, vr);
        index.fbeforeeq = beforeeq(index.fndx);
    end
    
end


fisacc1=findstrcmpndx(parsedmodel,'''',index,index.fbeforeeq);
fisacc=findstrcmpndx(parsedmodel,'''',index);

%isacc=strcmp(parsedmodel,'''');
parsedspec(fisacc) = vr.accent;
parsedspec(fisacc1) = vr.diff;

isdiff = ~isempty(fisacc1);

ft = findstrcmpndx(parsedmodel, vr.t,index);
ggrind.solver.nonautonomous = false;
if ~isempty(ft)
    %parsedspec(ft) = vr.t;
    [fdiff] = findstrcmpndx(parsedmodel, vr.t, index, index.fbeforeeq);
    ggrind.solver.isdiffer = ~isempty(fdiff);
    if isdiff&&ggrind.solver.isdiffer
        ggrind.errors{end + 1} = 'Cannot combine difference equation with differential equation';
    end
    
    if isdiff
        ggrind.solver.nonautonomous = true;
    else
        [lineswt,~] = findstrcmpndx(parsedmodel, vr.t, index);
        lineswt = unique(lineswt);
        for i = 1:length(lineswt)
            %remove (t+1) or (t-1), except  *(t+1) or other oper. In case of ambiguity use double brackets: ((t+1))
            f1 =  strfind(parsedspec(lineswt(i), :), [vr.brack1 vr.t vr.oper vr.number vr.brack2]);
            for j = 1:length(f1)
                if ~(f1(j)>1&&any(parsedspec(lineswt(i), f1(j)-1)==[vr.oper vr.brack1])) %exception if you have *(t+1) or ((t+1))
                    parsedspec(lineswt(i), f1(j):f1(j) + 4) = vr.empty;
                end
                
            end
            
            %remove (t), except *(t) or other oper. In case of ambiguity use double brackets: ((t))
            f1 =  strfind(parsedspec(lineswt(i), :), [vr.brack1 vr.t  vr.brack2]);
            for j = 1:length(f1)
                if ~(f1(j)>1&&any(parsedspec(lineswt(i), f1(j)-1)==[vr.oper vr.brack1])) %exception if you have *(t) or ((t))
                    parsedspec(lineswt(i), f1(j):f1(j) + 2) = vr.empty;
                end
                
            end
            
        end
        parsedspec(fdiff) = vr.diff;
    end
    
else
    ggrind.solver.isdiffer = false;
end

if any(index.utypes == vr.setodefile)
    [fi, fj]=find(parsedspec == vr.setodefile);
    pars = getfunctionpars(parsedmodel,parsedspec, fi(1),fj(1),vr);
    ggrind=setodefile(pars{:});
    ggrind.errors={};
    if ~strcmp(amodel{1},'%external odefile')
        if size(amodel,1)>1
            amodel=amodel';
        end
        amodel=[{'%external odefile'} amodel];
    end
    ggrind.model=amodel(:)';
    return;
end
%replace double quote with single (use double quotes in implicit models to


if ~isdiff&&~ggrind.solver.isdiffer
    %is implicit?
    [fimplic, fimplic2] = findstrcmpndx(parsedmodel, vr.implicitvars,index);
    if ~isempty(fimplic)
        flenacc= cellfun(@length,strfind(parsedmodel,''''));
        if any(flenacc(:)==2)
            %an even number of y' in an equation is interpreted as string
            %in implicit model (DAE) this should not be done
            %here this is repaired
            fstring=flenacc==2;
            [c,~]=find(fstring); %find columns
            c=unique(c);
            repairstr=cell(length(c));
            maxlen=size(parsedmodel,2);
            oldlen=maxlen;
            for i=1:length(c)
                line=sprintf('%s',parsedmodel{c(i),:});
                line=[line ''''];
                repairstr{i}=parsed_equation(line);
                if length(repairstr{i}.fs)-1>maxlen
                    maxlen=length(repairstr{i}.fs)-1;
                end
                
            end
            
            if maxlen>oldlen
                for i=1:size(parsedmodel,1)
                    parsedmodel(i,oldlen+1:maxlen)={''};
                    parsedspec(i,oldlen+1:maxlen)= vr.empty;
                end
                
            end
            
            for i=1:length(c)
                ffs=repairstr{i}.fs;
                ttypes=repairstr{i}.types;
                n=length(ffs)-1;
                parsedmodel(c(i),1:n)=ffs(1:n);
                parsedspec(c(i),1:n)=ttypes(1:n);
            end
            
            [index] = makeindex(parsedmodel, parsedspec, vr);
            [parsedspec,parsedmodel,index,beforeeq]=updatespecs(parsedspec,parsedmodel,index,vr,keywords,typkw);
            %   fimplicit=findstrcmpndx(parsedmodel, vr.implicitvars, index);
            %    parsedspec(fimplicit)=vr.implicitvars;
            ft = findstrcmpndx(parsedmodel, 't',index);
            parsedspec(ft)=vr.t;
        else
            parsedspec(fimplic,fimplic2)=vr.implicitvars;
        end
        fisacc=findstrcmpndx(parsedmodel,'''',index);
        parsedspec(fisacc) = vr.diff;
        for i=1:length(fimplic)
            beforeeq(fimplic(i),parsedspec(fimplic(i),:)==vr.variable) = true;
        end
        
        index.fbeforeq=beforeeq(index.fndx);
        ggrind.solver.nonautonomous= any(any(parsedspec==vr.t));
        ggrind.solver.isimplicit=true;
    else
        ggrind.errors{end + 1} = 'Not any model equations entered';
    end
    
end


%avoid confusion)
fstring=find(strncmpndx(parsedmodel, '"', 1, index));
if ~isempty(fstring)
    for i=1:length(fstring)
        s=parsedmodel{fstring(i)};
        s(s=='"')='''';
        parsedmodel{fstring(i)}=s;
        parsedspec(fstring(i))=vr.string;
    end
    
end


lineswdiff=repmat(any(parsedspec == vr.diff, 2), 1, size(parsedspec, 2));
if ggrind.solver.isimplicit
    lineswdiff=repmat(any(parsedspec == vr.implicitvars, 2), 1, size(parsedspec, 2));
end

beforeeqdiff = beforeeq & lineswdiff;

index.utypes = parsedspec(index.fndx(index.undxs1));
fbeforeeqdiff = beforeeqdiff(index.fndx);
fbeforeeq = beforeeq(index.fndx);
varnrs=find(index.utypes >= vr.funcname);
ggrind.explvector = false;
%Main loop to assign the types of variables.
haswrongorder = false(size(index.uelems));

index.uvarfirstass = zeros(size(varnrs)) - 1;
for k = 1:length(varnrs)
    i = varnrs(k);
    %  uvar = index.uelems(i); %convenient with debugging
    andx = transpose(index.undxs1(i):index.undxs2(i));
    avarndx = index.fndx(andx);
    assigned = fbeforeeq(andx);
    if any(assigned)
        index.utypes(i) = vr.auxil;
        iv = index.fi(andx(assigned));
        abeforeeqdiff = fbeforeeqdiff(andx(assigned));
        ms = abeforeeqdiff;
        if any(abeforeeqdiff)
            if any(~abeforeeqdiff)
                disp('Warning: it is not a good practice to assign to state variables, it is better to use auxiliary variables');
            end
            
            dd.dim1 = 0;
            dd.dim2 = 0;
            for m = 1:length(ms)
                if ms(m)
                    index.utypes(i) = vr.statevar;
                    elems = parsedspec(iv(m), :);
                    elems = elems(beforeeq(iv(m), :));
                    selems = parsedmodel(iv(m), :);
                    selems = selems(beforeeq(iv(m), :));
                    f = strfind(elems, [vr.brack1 vr.number vr.brack2]);
                    colons=find(elems == vr.colon);
                    isvector =  ~isempty(colons);
                    if isvector
                        d.dim1 = str2double(selems{colons(1) + 1});
                        f1=find(elems == vr.brack1);
                        f2=find(elems == vr.brack2);
                        if length(f1) == 1
                            parsedspec(iv(m), f1:f2) = vr.empty;
                        end
                        
                        if length(colons) > 1
                            d.dim2 = str2double(selems{colons(2) + 1});
                        else
                            d.dim2 = 1;
                        end
                        
                        if d.dim1 * d.dim2 == 1
                            isvector = false;
                        end
                        
                    elseif ~isempty(f)
                        ggrind.explvector = true;
                        d.dim1 = str2double(selems{f(1) + 1});
                        d.dim2 = 1;
                    else
                        d.dim1 = 1;
                        d.dim2 = 1;
                    end
                    
                    if dd.dim1 < d.dim1
                        dd.dim1 = d.dim1;
                    end
                    
                    if dd.dim2 < d.dim2
                        dd.dim2 = d.dim2;
                    end
                    
                end
                
            end
            
            ggrind.statevars.dims{end + 1} = d;
            ggrind.statevars.vector = ggrind.statevars.vector||isvector||ggrind.explvector;
            index.uvarfirstass(i) = iv(1);
        else
            iva = index.fi(andx(~assigned));
            firstuse = min(iva);
            index.uvarfirstass(i) = min(iv);
            parsedspec(avarndx) = vr.auxil;
            if ~isempty(firstuse)&&(firstuse <= index.uvarfirstass(i))
                haswrongorder(i) = true;
            end
            
        end
        
    else
        %         kw = strcmp(index.uelems{i}, keywords);
        %         if any(kw)
        %             f = find(kw, 1);
        %             index.utypes(i) = typkw(f);
        %             parsedspec(avarndx) = typkw(f);
        %         else
        if all(parsedspec(avarndx) == vr.funcname)
            if strcmp(index.uelems{i}, ggrind.locfuncname)
                index.utypes(i) = vr.locfun;
            elseif ~isempty(which(index.uelems{i}))
                index.utypes(i) = vr.funcname;
            else
                index.utypes(i) = vr.parameter;
                parsedspec(avarndx) = vr.parameter;
                
            end
            
        elseif index.utypes(i)==vr.variable
            index.utypes(i) = vr.parameter;
            parsedspec(avarndx) = vr.parameter;
        end
        
    end
    
    %     end
    
end

firststateass=index.uvarfirstass(index.utypes == vr.statevar);
ggrind.statevars.names=transpose(index.uelems(index.utypes == vr.statevar));
[~, ndx1] = sort(firststateass);
if ggrind.statevars.vector
    ggrind.statevars.vectnames = ggrind.statevars.names(ndx1);
    ggrind.statevars.names = {};
    ggrind.statevars.dims = ggrind.statevars.dims(ndx1);
    offset = 0;
    for i = 1:length(ggrind.statevars.dims)
        d = ggrind.statevars.dims{i};
        d.from = offset + 1;
        d.to = offset + d.dim1 * d.dim2;
        offset = d.to;
        ggrind.statevars.dims{i} = d;
    end
    ggrind.statevars.dim = offset;
else
    ggrind.statevars.names = ggrind.statevars.names(ndx1);
    ggrind.statevars.vectnames = {};
    ggrind.statevars.dims = {};
    ggrind.statevars.dim = length(ggrind.statevars.names);
end

if ggrind.solver.isimplicit
    for j=1:size(parsedspec,1)
        fstatep=strfind(parsedspec(j,:),[vr.variable vr.diff]);
        for k=1:length(fstatep)
            parsedspec(j,fstatep(k))=vr.statevar_p;
            parsedspec(j,fstatep(k)+1)=vr.empty;
            parsedmodel{j,fstatep(k)+1}='';
        end
        
    end
    
end


if any(index.utypes == vr.defextern)
    [fi, fj]=find(parsedspec == vr.defextern);
    ggrind.externvars = cell(length(fi), 1);
    for i = 1:length(fi)
        pars = getfunctionpars(parsedmodel,parsedspec, fi(i),fj(i),vr);
        f =  strcmp(index.uelems, pars{1});
        if ~isempty(f)
            index.utypes(f) = vr.extern;
        end
        
        e.name = pars{1};
        if length(pars) > 1
            e.default = pars{2};
            if length(pars) > 2
                e.options = pars(3:end);
            end
            
        else
            e.default = '0';
        end
        ggrind.externvars{i} = e;
    end
    
end


if any(index.utypes == vr.defpermanent)
    [fi, fj]=find(parsedspec == vr.defpermanent);
    for i = 1:length(fi)
        pars =  getfunctionpars(parsedmodel,parsedspec, fi(i),fj(i),vr);
        f =  strcmp(index.uelems, pars{1});
        indexcode = index.uelemnr(f);
        parsedspec(index.elems == indexcode)=vr.permanent;
        if ~isempty(f)
            index.utypes(f) = vr.permanent;
            perm.name = pars{1};
            if length(pars) > 1
                perm.currentval = pars{2};
            end
            
            ggrind.permanent{i} = perm;
            haswrongorder(f) = false;
        end
        
    end
    
end

ggrind.definespace = {};
if any(index.utypes == vr.definespace)
    [fi, ~]=find(parsedspec == vr.definespace);
    ggrind.definespace = cell(length(fi), 1);
    for i = 1:length(fi)
        ggrind.definespace{i} = sprintf('%s', parsedmodel{fi(i), :});
        parsedspec(fi(i), :) = vr.empty;
    end
    
end

%make an unique number for each block of code (only level 1)
startblock=sum(parsedspec == vr.block1 | parsedspec==vr.if, 2);
endblock=sum(parsedspec == vr.block2, 2);
endblock = [0; endblock(1:end - 1)];
blocklevel = cumsum(startblock - endblock);
blocknr = zeros(size(startblock));
bnr = 1;
for i = 1:length(startblock)
    if startblock(i) > 0||blocklevel(i)==0
        bnr = bnr + 1;
    end
    
    blocknr(i) = bnr;
end

diffblocks=sum(parsedspec == vr.diff, 2);
if ggrind.solver.isimplicit
    %for implicit model remove lines with 0 in the first position
    %and lines with vr.implicitvars
    diffblocks=parsedspec(:,1) == vr.number|sum(parsedspec==vr.implicitvars,2)>0;
end

diffblocks = unique(blocknr(diffblocks > 0));

for i = 1:length(diffblocks)
    blocknr(diffblocks(i) == blocknr)=diffblocks(i) + vr.mindiff;
end

deleterows=sum(parsedspec==vr.definepars | parsedspec==vr.defextern | parsedspec==vr.defpermanent, 2) > 0;
blocknr(deleterows) = 0;
ggrind.pars=index.uelems(index.utypes == vr.parameter);
if size(ggrind.pars, 1) > 1
    ggrind.pars = transpose(ggrind.pars);
end

p = lower(ggrind.pars);
[~, ndx1] = sort(p);
ggrind.pars = ggrind.pars(ndx1);
ggrind.funcnames.names=index.uelems(index.utypes == vr.auxil);
if ~isempty(ggrind.funcnames.names)
    firstuse=index.uvarfirstass(index.utypes == vr.auxil);
    [~, ndx] = sort(firstuse);
    ggrind.funcnames.names = transpose(ggrind.funcnames.names(ndx));
end

%if ~isempty(ggrind.funcnames.names) %for testing
if any(haswrongorder)
    %sort auxvars
    disp('sorting auxiliary variables');
    selaux=parsedspec == vr.auxil;
    %    auxblocks=auxblocks(auxblocks<vr.mindiff);
    assignedaux = index.elems.*(selaux & beforeeq);
    assignedaux=assignedaux(:, sum(assignedaux, 1) ~= 0);
    doubleaux=~([true; diff(assignedaux(:, 1))] | assignedaux(:, 1)==0);
    auxblocks = unique(blocknr(any(selaux, 2) & blocknr < vr.mindiff));
    usedaux = index.elems.*(selaux & ~beforeeq);
    usedaux = sort(usedaux, 2);
    usedaux=usedaux(:, sum(usedaux, 1) ~= 0); %auxilaries that are used per block
    assignedset = 0;
    block = zeros(size(auxblocks));
    k = 0;
    usedblocks = false(size(auxblocks));
    %loop to sort the blocks, loop till nothing is changed
    ischanged = true;
    while ischanged
        ischanged = false; %this is also safe if not all blocks can be assigned
        i = 1;
        while i  <= length(auxblocks)
            usedau=usedaux(auxblocks(i) == blocknr, :);
            if size(usedau, 1) > 1
                %more than one lines in the block, it can be that it is
                %self-assigning (don't check details)
                assignedau=assignedaux(auxblocks(i) == blocknr, :);
                usedau = usedau(~ismember(usedau(:), assignedau(:)));
            end
            
            usedau = usedau(usedau > 0);
            % a block of code can be added if all used auxiliaries are already
            % defined
            if ~usedblocks(i)&&all(ismember(usedau, assignedset))
                while 1  % if there are two the same assigned auxiliary variables both lines are added
                    usedblocks(i) = true;
                    ischanged = true;
                    k = k + 1;
                    block(k) = auxblocks(i);
                    assignedset=union(assignedset, assignedaux(auxblocks(i) == blocknr, :));
                    if (i==length(auxblocks))||usedblocks(i+1)||~any(doubleaux(auxblocks(i+1)==blocknr))
                        break;
                    end
                    
                    i = i + 1;
                end
                
            end
            
            i = i + 1;
        end
        
    end
    
    if all(usedblocks)
        %only if all blocks are assigned we can trust the result
        [~, ndx] = sort(block);
        auxblocks1 = auxblocks(ndx);
        blockndx = any(selaux, 2) & blocknr < vr.mindiff;
        if sum(blockndx) == length(auxblocks) % no blcoks  > 1
            blocknr(blockndx) = auxblocks1;
        else
            blocknr1 = blocknr;
            for i = 1:length(auxblocks)
                blocknr1(blocknr == auxblocks(i)) = auxblocks1(i);
            end
            
            blocknr = blocknr1;
        end
        
    else
        %write the partly sorted list of auxiliaries to help finding the
        %circular dependence
        fprintf(2, 'Cannot sort all auxiliary variables:\nThe following auxiliary variables are sorted in sequential dependence:\n');
        f = find(selaux);
        for i = 1:length(f)
            parsedmodel{f(i)}=sprintf('<a href="matlab:%%" style="font-weight:bold">%s</a>', parsedmodel{f(i)});
        end
        
        usedb=block(block ~= 0);
        if isempty(usedb)
            disp('none');
        else
            for i = 1:length(usedb)
                ndx=find(blocknr == usedb(i));
                for j = 1:length(ndx)
                    disp(writeline(ndx(j), parsedspec, parsedmodel, vr));
                end
                
            end
            
        end
        
        fprintf(2,'\nThese cannot be sorted, probably because of circular dependence:\n');
        unusedb = auxblocks(~usedblocks);
        for i = 1:length(unusedb)
            ndx=find(blocknr == unusedb(i));
            for j = 1:length(ndx)
                disp(writeline(ndx(j), parsedspec, parsedmodel, vr));
            end
            
        end
        
        ggrind.errors{end + 1} = 'Cannot sort auxuliary variables, possibly due to circular dependence or double assignments or two statements on one row';
    end
    
end


% if nargin > 1
%    if g_grind.statevars.vector
%       statenames1 = g_grind.statevars.vectnames;
%       statenames2 = ggrind.statevars.vectnames;
%    else
%       statenames1 = g_grind.statevars.names;
%       statenames2 = ggrind.statevars.names;
%    end

%    haserror=(length(g_grind.pars)~=length(ggrind.pars))|| any(~strcmp(g_grind.pars, ggrind.pars));
%    haserror= haserror || (length(statenames1)~=length(statenames2))||any(~strcmp(statenames1, statenames2));
%    if 0 % the current grind has too many errors to check this
%       uoldfuncs = unique(g_grind.funcnames.names); %old names are not unique
%       unewfuncs = sort(ggrind.funcnames.names);
%       if ~isempty(uoldfuncs)&&~isempty(unewfuncs)
%          haserror= haserror || (length(uoldfuncs)~=length(unewfuncs))||any(~strcmp(uoldfuncs, unewfuncs));
%       end

%    end
%    if haserror
%       error('test:analysemodel','params not correct');
%    end

% end


% ndx1 = ndx(b < vr.mindiff);
% funcs = cell(size(ndx1));
% for i = 1:length(ndx1)
%   funcs{i} = printline(parsedmodel,parsedspec,ndx1(i),vr);
% end

% ggrind.funcs = sprintf('%s\n', funcs{:});
parsedspec(parsedspec == vr.diff)=vr.empty;
%ggrind.parsedmodel = parsedmodel(ndx, :);
[b, fndx] = sort(blocknr);
ggrind.modelndx=fndx(b~=0);
fndx=fndx(b ~= 0 & b < vr.mindiff);
ggrind.funcs = cell(length(fndx), 1);
kf = 1;
for i = 1:size(fndx, 1)
    s1 = writeline(fndx(i), parsedspec, parsedmodel, vr);
    if ~isempty(s1)
        ggrind.funcs{kf} = s1;
        kf = kf + 1;
    end
    
end

ggrind.funcs = sprintf('%s\n', ggrind.funcs{1:kf - 1});
%Adapt the code of parsedmodel such that we can write an odefile:
ggrind = makeodefile(ggrind, parsedmodel, parsedspec, blocknr, index, vr, beforeeqdiff);


function ggrind = makeodefile(ggrind, parsedmodel, parsedspec, blocknr, index, vr, beforeeqdiff)
fbeforeeqdiff = beforeeqdiff(index.fndx);
ggrind.solver.name = 'ode45';
odefile = cell(length(parsedmodel) + 20, 1);
fparsedspec = parsedspec(index.fndx);
felemnr = index.elems(index.fndx);
k = 1;
% if ~ggrind.newodefile
%    odefile{k} = '%function created by GRIND'; k = k + 1;
%    odefile{k}='function g_X2=curr_ode0(t,g_X1)';kfun=k; k=k + 1;
%    odefile{k} = ''; kglob = k; k = k + 1;
%    sglobals = ggrind.pars;
% else
odefile{k} = '%function created by GRIND'; k = k + 1;
if isempty(ggrind.pars)
    odefile{k}=sprintf('function g_X2=curr_ode0(t,g_X1)');  kfun=k; k=k + 1;
    ppar = '';
else
    ppar = sprintf(',%s',ggrind.pars{:});
    odefile{k}=sprintf('function g_X2=curr_ode0(t,g_X1%s)',ppar);  kfun=k; k=k + 1;
end

odefile{k} = ''; kglob = k; k = k + 1;
sglobals = {};
%end

if ~isempty(ggrind.externvars) ||~isempty(ggrind.permanent)||ggrind.statevars.vector
    sglobals = [sglobals {'g_grind'}];
    if ~isempty(ggrind.permanent)
        odefile{k} = 'i_updatepermanent(t);';k=k + 1;
        for i = 1:length(ggrind.permanent)
            p = ggrind.permanent{i};
            iuelem = strcmp(p.name, index.uelems);
            elnr = index.uelemnr(iuelem);
            parsedmodel(index.fndx(felemnr == elnr)) ={sprintf('g_grind.permanent{%d}.currvalue', i)};
            %        permndx = strcmpndx(parsedmodel, p.name, index);
            %        parsedmodel(permndx) = {sprintf('g_grind.permanent{%d}.currvalue', i)};
        end
        
    end
    for i = 1:length(ggrind.externvars)
        e = ggrind.externvars{i};
        odefile{k} = sprintf('%s=externvar(%d,%s,t);',e.name,i,e.default);k = k + 1;
    end
    
end

%analyse lags
flagndx=fparsedspec == vr.lag;
ggrind.solver.haslags = any(flagndx);
if ggrind.solver.haslags
    if ggrind.statevars.vector
        ggrind.errors{end + 1} = 'Lags are not yet supported in vector models';
    end
    
    if ggrind.solver.isdiffer
        ggrind.errors{end + 1} = 'Lags are not supported in difference equations';
    end
    
    %[lagi, lagj] = findstrcmpndx(parsedmodel, 'lag', index);
    lagi = index.fi(flagndx);
    lagj = index.fj(flagndx);
    ggrind.dde.lags = {};
    ggrind.dde.isvariable=false;
    ggrind.solver.history = [];
    changed=false;
    for i = 1:length(lagi)
        [pars,js] =  getfunctionpars(parsedmodel,parsedspec, lagi(i),lagj(i),vr);
        avar = pars{1};
        isstrng = false;
        if  strcontains(avar,'''')
            avar = avar(2:end - 1);
            isstrng = true;
        end
        
        ivar =  find(strcmp(avar, ggrind.statevars.names));
        parsedmodel{lagi(i), js(1,1)} = int2str(ivar);
        if ~isstrng
            parsedspec(lagi(i), js(1,1))  = 0;
            index.elems(lagi(i), js(1, 1)) =  max(max(index.elems)) + 1;
            felemnr = index.elems(index.fndx);
        end
        
        lagstr = sprintf('%s', pars{2});
        lagno = find(strcmp(lagstr, ggrind.dde.lags), 1);
        if isempty(lagno)
            lagno = length(ggrind.dde.lags) + 1;
            ggrind.dde.lags{lagno} = lagstr;
        end
        
        %the lag can only be a combination of parameters in dde23 for ddesd
        %more is possible (isvariable=true)
        typesoflag=parsedspec(lagi(i),js(2,1):js(2,2));
        if any(typesoflag==vr.t)||any(typesoflag==vr.variable)
            ggrind.dde.isvariable=true;
        end
        if any(typesoflag==vr.auxil)||any(typesoflag==vr.permanent)||any(typesoflag==vr.stochfun)
            error('grind:dde:typeoflag','Only parameters, time (t), or state variables can be used in the time delays')
        end
        
        parsedmodel{lagi(i), js(2,1)} = int2str(lagno);
        for j=js(2,1)+1:js(2,2)
            parsedmodel{lagi(i), j}='';
            parsedspec(lagi(i), j)  = 0;
            changed=true;
        end
        parsedmodel{lagi(i), lagj(i)} = 'g_lags';
    end
    if ggrind.dde.isvariable
        ggrind.solver.name='ddesolsd';
    else
        ggrind.solver.name='ddesol';
    end
    %    if ~any(strcmp('g_grind', sglobals))
    %       sglobals = [sglobals {'g_grind'}];
    %    end
    if changed
        index = makeindex(parsedmodel, parsedspec, vr);
        fbeforeeqdiff=beforeeqdiff(index.fndx);
        fparsedspec = parsedspec(index.fndx);
        felemnr = index.elems(index.fndx);
    end
    odefile{kfun}=sprintf('function g_X2=curr_ode0(t,g_X1,g_lags%s)',ppar);
    odefile{k}=sprintf('if isempty(g_lags)\n  g_lags=repmat(g_X1,1,%d);\nend',length(ggrind.dde.lags));k=k + 1;
end

%boxcartrain
%example:
% juv_=boxcartrain(juv,devjuv,sd);
% [juv_,juv]=boxcartrain(1,juv,devjuv,sd);
fbox=fparsedspec == vr.boxcartrain;
if any(fbox)
    is = index.fi(fbox);
    js = index.fj(fbox);
    ggrind.boxcar.names = cell(1, length(is));
    ggrind.boxcar.gcycl = zeros(1, length(is));
    for i = 1:length(is)
        pars =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
        j = 1;
        while j < size(parsedspec, 1)&&parsedspec(is(i), j)~=vr.auxil
            j = j + 1;
        end
        
        parsedmodel{is(i),j} = sprintf('[%s,%s]',parsedmodel{is(i),j},pars{1});
        parsedmodel{is(i),js(i) + 1} = sprintf('(%d,',i);
        ggrind.boxcar.names{i} = pars{1};
    end
    ggrind.solver.eulerneeded=true; 
end

%setevent
if any(fparsedspec == vr.setevent)
    ggrind.errors{end + 1} = '"setevent" should be in the lower parameters panel';
end

%implicitdisperse
if ggrind.statevars.vector
    fimplicit=fparsedspec == vr.implicitdisperse;
    if any(fimplicit)
        is = index.fi(fimplicit);
        js = index.fj(fimplicit);
        ggrind.implicdisp = cell(1, length(js));
        for i = 1:length(is)
            [pars,js1] =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
            parsedmodel{is(i), js1(1, 1)} = int2str(i);
            parsedmodel(is(i), js1(1, 1) + 1:js1(1, 2)) = {''};
            parsedmodel{is(i), js1(2, 1)}  = pars{1};
            parsedmodel(is(i), js1(2, 1) + 1:js1(2, 2)) = {''};
            p.Name = pars{1};
            p.D = pars{2};
            ggrind.implicdisp{i} = p;
        end
        
        ggrind.solver.name = 'euler';
        ggrind.solver.eulerneeded=1;
    end
    
end



% globalpar = true;
% if globalpar
%    sglobals = [ggrind.pars sglobals];
% else
%    for i = 1:length(ggrind.pars)
%       fndx = findstrcmpndx(parsedmodel, ggrind.pars{i}, index);
%       parsedmodel(fndx) = {sprintf('g_grind.parvalues{%d}', i)};
%    end
%    %     ppars=index.fndx(fparsedspec == vr.parameter);
%    %     for i = 1:length(ppars)
%    %        parsedmodel{ppars(i)} = sprintf('g_grind.parvalues.%s', parsedmodel{ppars(i)});
%    %     end

%    if ~any(strcmp('g_grind', sglobals))
%       sglobals = [sglobals {'g_grind'}];
%    end
% end

%externlag
fexternlag=fparsedspec == vr.externlag;
if any(fexternlag)
    is = index.fi(fexternlag);
    js = index.fj(fexternlag);
    for i = 1:length(is)
        [pars,js1] =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
        par1=regexp(pars{1},'[^'']*','match','once');
        ivar=1;
        while (ivar<=length(ggrind.externvars))&&(~strcmp(par1,ggrind.externvars{ivar}.name))
            ivar=ivar+1;
        end
        doclear=false(1,size(parsedmodel,2));
        if ivar>length(ggrind.externvars)
            warning('grind:externlag','The function "externlag" ignored, it works only for external variables (see <a href="matlab:help defextern">defextern</a>)');
            doclear(js1(1,1)-2:js1(1,1)-2)=true;
            doclear(js1(2,1)-1:js1(end,2))=true;
        else
            parsedmodel{is(i),js1(1,1)}=sprintf('%d',ivar);
            parsedmodel{is(i),js(i)}='externvar';
            parsedmodel{is(i),js1(end,2)}=sprintf('%s,t-(%s)',ggrind.externvars{ivar}.default,pars{2});
            doclear(js1(end,1):js1(end,2)-1)=true;
        end
        parsedmodel(is(i),doclear)={''};
        parsedspec(is(i),doclear)=-99;
%         if ~strncmp(pars{1},'''',1)
%             parsedmodel{is(i),js1(1,1)}=sprintf('''%s''',pars{1});
%         end

        % parsedmodel(is(i),js1(1,1)+1:js1(end))={''};
    end   
end

if ~isempty(intersect(fparsedspec,[vr.dwiener,vr.djump,vr.rednoise, vr.stochfun]))
    ggrind.solver.isstochastic = true;
    %dwiener
    fdwiener=fparsedspec == vr.dwiener;
    if any(fdwiener)
        if ggrind.solver.isdiffer
            ggrind.errors{end + 1} = '"dwiener" cannot be used in difference equations';
        end
        is = index.fi(fdwiener);
        js = index.fj(fdwiener);
        [~,ndx]=sortrows([is,js],[1,2]);
        is=is(ndx);
        js=js(ndx);
        oldis=is(1);
        numj=0;
        numi=1;
        for i = 1:length(is)
            if is(i)==oldis
                numj=numj+1;
            else
                oldis=is(i);
                numj=1;
                numi=numi+1;
            end
            [pars,js1] =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
            doclear=false(1,size(parsedmodel,2));
            doclear(js1(1):js1(end))=true;
            pars=regexprep(pars,'^''|''$','');
            dwien=dwiener(pars{:});
            ndx1=find(strcmp(pars,dwien.fun),1);
            doclear(js1(ndx1,1):js1(ndx1,2))=false;
            parsedmodel{is(i),js1(1,1)-1}='(t,';
            parsedmodel{is(i),js1(end,2)+1}=sprintf(',%d%s',i,parsedmodel{is(i),js1(end,2)+1});
            if isfield(dwien,'dfundx')&&~isempty(dwien.dfundx)
                ndx1=find(strcmp(pars,dwien.dfundx),1);
                doclear(js1(ndx1,1)-1:js1(ndx1,2))=false;
            end
            parsedmodel(is(i),doclear)={''};
            parsedspec(is(i),doclear)=-99;
            if ~isfield(dwien,'corr')
                dwien.corr=[];
            end
            ggrind.solver.dwiener.args(i)=struct('corr',dwien.corr,'L',[],'var',numi,'num',numj, 'alpha',dwien.alpha,'beta',dwien.beta);
            % parsedmodel(is(i),js1(1,1)+1:js1(end))={''};
        end
        ggrind.solver.eulerneeded=1;
    end
    
    if any(fparsedspec == vr.dwiener)
        ggrind.solver.eulerneeded=1;
    end
    
    %djump
    fdjump=fparsedspec == vr.djump;
    if any(fdjump)
        is = index.fi(fdjump);
        js = index.fj(fdjump);
        [~,ndx]=sortrows([is,js],[1,2]);
        is=is(ndx);
        js=js(ndx);
        oldis=is(1);
        numj=0;
        for i = 1:length(is)
            if is(i)==oldis
                numj=numj+1;
            else
                oldis=is(i);
                numj=1;
            end
            [pars,js1] =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
            doclear=false(1,size(parsedmodel,2));
            doclear(js1(1):js1(end))=true;
            pars=regexprep(pars,'^''|''$','');
            djum=djump(pars{:});
            %             ndx1=find(strcmp(pars,djum.fun),1);
            %             doclear(js1(ndx1,1):js1(ndx1,2))=false;
            parsedmodel{is(i),js(i)}='i_djump';
            jumppars=[djum.timing_pars(:); djum.size_pars(:)];
            if ggrind.statevars.vector
            else
                for k1=1:length(ggrind.statevars.names)
                    f1=strcmp(jumppars,ggrind.statevars.names{k1});
                    if any(f1)
                        jumppars(f1)={sprintf('g_X1(%d,:)',k1)};
                    end
                end
            end
            parslist=sprintf('%s,',jumppars{:});
            if ~isempty(parslist)
                parslist=parslist(1:end-1);
                parsedmodel{is(i),js1(1,1)-1}=sprintf('(%d,t,%s',i,parslist);
            else
                parsedmodel{is(i),js1(1,1)-1}=sprintf('(%d,t',i);
            end
            %parsedmodel{is(i),js1(end,2)+1}=sprintf('%d%s',i,parsedmodel{is(i),js1(end,2)+1});
            parsedmodel(is(i),doclear)={''};
            parsedspec(is(i),doclear)=-99;
            if ~isfield(djum,'corr')
                djum.corr=[];
            end
            djum.var=is(i);
            djum.num=numj;
            djum.nextt=[];
            ggrind.solver.djump.args(i)=orderfields(djum);
            % parsedmodel(is(i),js1(1,1)+1:js1(end))={''};
        end
        ggrind.solver.eulerneeded=1;
    end
    
    if any(fparsedspec == vr.djump)
        ggrind.solver.eulerneeded=1;
    end
    
end

if ggrind.solver.isdiffer
    ggrind.solver.name = 'i_differ';
end

if ggrind.solver.eulerneeded
    ggrind.solver.name='euler';
end

if ggrind.statevars.vector
    odefile{k}='g_X2=zeros(g_grind.statevars.dim,1);'; k=k + 1;
    for i = 1:length(ggrind.statevars.vectnames)
        if  ggrind.statevars.dims{i}.dim2 == 1
            odefile{k}= sprintf('%s = g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to);',ggrind.statevars.vectnames{i}, i, i);
        else
            odefile{k}= sprintf('%s = reshape(g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to),g_grind.statevars.dims{%d}.dim1,g_grind.statevars.dims{%d}.dim2);', ...
                ggrind.statevars.vectnames{i}, i, i, i ,i);
        end
        
        k = k + 1;
        iuelem = strcmp(ggrind.statevars.vectnames{i}, index.uelems);
        elnr = index.uelemnr(iuelem);
        if ggrind.explvector
            andx=felemnr==elnr & fbeforeeqdiff;
            iv = index.fi(andx);
            jv = index.fj(andx);
            for j = 1:length(iv)
                f = strfind(parsedspec(iv(j),:), [vr.brack1 vr.number vr.brack2]);
                if ~isempty(f)
                    parsedspec(iv(j), f(1):f(1) + 2) = vr.empty;
                    aval = str2double(parsedmodel(iv(j), f(1) + 1));
                    parsedmodel{iv(j),jv(j)} = sprintf('g_X2(g_grind.statevars.dims{%d}.from+%d,1)',i,aval - 1);
                end
                
            end
            
        else
            parsedmodel(index.fndx(felemnr==elnr & fbeforeeqdiff)) ={sprintf('g_X2(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to)', i, i)};
        end
        
    end
    
else
    for i = 1:length(ggrind.statevars.names)
        iuelem = strcmp(ggrind.statevars.names{i}, index.uelems);
        elnr = index.uelemnr(iuelem);
        %      parsedmodel(index.fndx(felemnr==elnr & fbeforeeqdiff)) ={sprintf('g_X2(%d,1)',i)};
        %      parsedmodel(index.fndx(felemnr==elnr & ~fbeforeeqdiff)) ={sprintf('g_X1(%d)', i)};
        parsedmodel(index.fndx(felemnr==elnr & fbeforeeqdiff)) ={sprintf('g_X2(%d,:)',i)};
        parsedmodel(index.fndx(felemnr==elnr & ~fbeforeeqdiff)) ={sprintf('g_X1(%d,:)', i)};
        opers={'*','^','/'};
        for j=1:length(opers)
            indx = strcmp(opers{j}, index.uelems);
            elnr=index.uelemnr(indx);
            if ~isempty(elnr)
                fndx=index.fndx(felemnr==elnr);
                ndx=strcmp(parsedmodel(fndx),opers{j});  %lags can have emptied opers that are still in de index
                parsedmodel(fndx(ndx))={['.' opers{j}]};
            end
        end
        if isempty(intersect({'Dxx','Dx','Dyy','Dy','iif','if','neighborcells','dwiener','switch','sum'}, index.uelems))
            ggrind.solver.opt.Vectorized='on';
        else
            ggrind.solver.opt.Vectorized='off';
        end
        
    end
    
end

if ggrind.solver.isimplicit
    ggrind.solver.opt.Vectorized='off';
    [f1,f2]=find(parsedspec==vr.statevar_p);
    for i = 1:length(f1)
        s=parsedmodel{f1(i),f2(i)};
        if strcontains(s,'g_X1')
            s(4)='3';
            parsedmodel{f1(i),f2(i)}=s;
        end
        
    end
    
    dellines=index.fi(fparsedspec==vr.implicitvars);
    for i=1:length(dellines)
        parsedmodel(dellines(i),:)={''};
        parsedspec(dellines(i),:)=vr.empty;
    end
    
    odefile{kfun}=sprintf('function g_res=curr_ode0(t,g_X1,g_X3%s)',ppar);
    parsedmodel(index.fndx(fparsedspec==vr.number & index.fbeforeeq))={'g_res'};
    ggrind.solver.name='i_ode15sol';
end

[b, ndx] = sort(blocknr);
ndx=ndx(b ~= 0);
for i = 1:size(ndx, 1)
    s1 = writeline(ndx(i), parsedspec, parsedmodel, vr);
    if ~isempty(s1)
        odefile{k} = s1;
        k = k + 1;
    end
    
end

if ~isempty(sglobals)
    odefile{kglob} = sprintf('global %s;', strtrim(sprintf('%s ',sglobals{:})));
else
    odefile(kglob) = [];
    k = k - 1;
end

if ~isempty(ggrind.permanent)
    odefile{k} = 'i_updatepermanent;';k=k + 1;
end


ggrind.theodefile = odefile(1:k - 1);

function s = writeline(i, parsedspec, parsedmodel, vr)
s2 = parsedspec(i, :);
ndxs2=(s2~=vr.empty & s2~=vr.comment);
f = find(ndxs2, 1, 'last');
if s2(f) == vr.semicolon
    ndxs2(f) = false;
end

s1 = parsedmodel(i, ndxs2);
if ~isempty(s1)
    s2 = s2(ndxs2);
    %  if ~isempty(intersect([vr.else,vr.if],s2))]
    if any(s2==vr.else | s2==vr.if)
        s=sprintf('%s',sprintf('%s',s1{:}));
    else
        s=sprintf('%s;',sprintf('%s',s1{:}));
    end
    
else
    s = '';
end

function index = makeindex(parsedmodel, parsedspec, vr)
[b, ndx] = sort(parsedmodel(:));
if size(b, 1) == 1
    b = transpose(b);
    ndx = transpose(ndx);
end

nempty=parsedspec(ndx) ~= vr.empty & parsedspec(ndx) ~= vr.comment;
b = b(nempty);
index.fndx = ndx(nempty);
[index.fi, index.fj] = ind2sub(size(parsedspec), index.fndx);
ddiff = ~strcmp(b(1:end - 1), b(2:end));
sameelems = cumsum([1; ddiff(1:end)]);
index.uelems = b([ddiff; true]);
index.undxs1 = find([true; ddiff]);
index.undxs2 = find([ddiff; true]);
index.uelemnr = transpose(1:length(index.undxs1));
index.elems = zeros(size(parsedmodel));
index.elems(index.fndx) = sameelems;
index.utypes = parsedspec(index.fndx(index.undxs1));

function [ndx1, ndx2] = findstrcmpndx(~, elem, index, afilter)
if ischar(elem)
    iuelem = find(strcmp(elem, index.uelems));
else
    iuelem = find(index.utypes==elem);
end
if ~isempty(iuelem)
    if nargin == 4
        endx = false(size(index.fndx));
        endx(index.undxs1(iuelem):index.undxs2(iuelem)) = true;
        endx(~afilter) = false;
        if nargout == 2
            ndx1 = index.fi(endx);
            ndx2 = index.fj(endx);
        else
            ndx1 = index.fndx(endx);
        end
        
    else
        if nargout == 2
            ndx1 = index.fi(index.undxs1(iuelem):index.undxs2(iuelem));
            ndx2 = index.fj(index.undxs1(iuelem):index.undxs2(iuelem));
        else
            ndx1 = index.fndx(index.undxs1(iuelem):index.undxs2(iuelem));
        end
        
    end
    
else
    ndx1 = [];
    ndx2 = [];
end



function [ndx] = strcmpndx(parsedmodel, elem, index)
iuelem = strcmp(elem, index.uelems);
elemnr = index.uelemnr(iuelem);
if ~isempty(elemnr)
    ndx=index.elems == index.uelemnr(iuelem);
else
    ndx = false(size(parsedmodel));
end

function [ndx] = strncmpndx(parsedmodel, elem, len, index)
iuelem = strncmp(elem, index.uelems,len);
elemnr = index.uelemnr(iuelem);
if ~isempty(elemnr)
    ndx=index.elems == index.uelemnr(iuelem);
else
    ndx = false(size(parsedmodel));
end

function [pars,js] =  getfunctionpars(parsedmodel,parsedspec, i,j,vr)
if parsedspec(i, j + 1) == vr.space
    %command style function
    specs = parsedspec(i, j + 1:end);
    js=find(specs == vr.space);
    jend=find(specs == vr.semicolon);
    if isempty(jend)
        jend = length(specs);
        while jend > 0&&specs(jend)==vr.empty
            jend = jend - 1;
        end
        
        jend = jend + 1;
    end
    
    jend = jend + j;
    js = js + j;
    js = [js'+1,[js(2:end)'; jend] - 1];
    pars = cell(size(js, 1), 1);
    for k = 1:size(js, 1)
        pars{k} = sprintf('%s', parsedmodel{i, js(k, 1):js(k, 2)});
        f= strfind(pars{k},'''');
        if length(f) == 2
            pars{k} = pars{k}(f(1) + 1:f(2) - 1);
        end
        
    end
    
else
    %brackets
    j = j + 2;
    poplevel = 1;
    pars = cell(5, 1);
    js = zeros(5, 2);
    ip = 1;
    while (j < size(parsedspec, 2))&&poplevel > 0
        k = 0;
        while (j+k < size(parsedspec, 2)) && poplevel~=0 && ~(poplevel==1&&any(parsedspec(i, j+k)== vr.comma))
            if any(parsedspec(i, j + k) == [vr.brack1,vr.curlbrack1])
                poplevel = poplevel + 1;
            end
            
            k = k + 1;
            if any(parsedspec(i, j + k) == [vr.semicolon, vr.brack2,vr.curlbrack2])
                poplevel = poplevel - 1;
            end
            
        end
        
        js(ip, :) = [j, j + k - 1];
        pars{ip} = sprintf('%s', parsedmodel{i, js(ip, 1):js(ip, 2)});
        ip = ip + 1;
        j = j + k + 1;
    end
    
    pars = pars(1:ip - 1);
    js = js(1:ip - 1, :);
end

function [parsedspec,parsedmodel,index,beforeeq]=updatespecs(parsedspec,parsedmodel,index,vr,keywords,typkw)


%parsedspec(findstrcmpndx(parsedmodel, 'for', index)) = vr.block1;

for i=1:length(keywords)
    kw = find(strcmp(index.uelems, keywords{i}),1);
    for j=1:length(kw)
        index.utypes(kw(j)) = typkw(i);
        parsedspec(index.fndx(index.undxs1(kw(j)):index.undxs2(kw(j)))) = typkw(i);
    end
end
%parsedspec(findstrcmpndx(parsedmodel, 'for', index)) = vr.block1;
%find elements before equal sign
iseq=strcmpndx(parsedmodel, '=',index);
parsedspec(iseq) = vr.assign;
ieq = sum(iseq, 2);
beforeeq = false(size(parsedspec));
for i = 1:length(ieq)
    if ieq(i) > 0
        fj = find(iseq(i, :));
        beforeeq(i, 1:fj(1)) = true;
        if length(fj) > 1
            fsemicolon=find(parsedspec(i, :) == vr.semicolon);
            for j = 2:length(fj)
                beforeeq(i, fsemicolon(j - 1):fj(j)) = true;
            end
            
        end
        
    end
    
    %f1 = find(iseq(f(i), :), 1);
end

index.fbeforeeq = beforeeq(index.fndx);



