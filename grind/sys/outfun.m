%OUTFUN   Flexible way to get output from the model
%   Make a matrix from the last run with any auxiliary/state variables from the model.
%   Note that there is no rerun of the model if parameters have been changed.
%
%
%   Usage:
%   A=OUTFUN('FUN') - makes a matrix A with the output FUN.
%   A=OUTFUN({'FUN1','FUN2'}) - makes a matrix A with the output FUN1 and FUN2.
%   OUTFUN('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'catfun' [struct or string or number] - Get categorial statistics of the state variable
%     'datastruc' [struct] - You can use a data structure like g_paranal.run to extract the data
%     'fun' [equation] - The equation(s) to show
%     'runfield' [string] - The name of the field in the datastruc that contains the run information (e.g. 'run' or 'prevrun')
%     'times' [number] - The times for outputs
%   OUTFUN('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-#ndxofparanal' - used internally to get only on value and the result in a correct shape
%     '-b' - use the last <a href="matlab:help contbif">contbif</a> result
%     '-c' - use the last <a href="matlab:help conteq">conteq</a> result
%     '-m' - use the last <a href="matlab:help mcarlo">mcarlo</a> results.
%     '-n' - use the last normal run
%     '-p' - use the last <a href="matlab:help paranal">paranal</a> or  <a href="matlab:help paranal2d">paranal2d</a> result
%     '-t' FUN TIMES - get a single output FUN at time TIMES. If the output is a matrix, it has a correct shape.
%
%   Examples:
%   A=OUTFUN('x*y/par') x, y = state variable par = parameter
%   A=OUTFUN('cons') cons = permanent variable
%   A=OUTFUN({'FUN1','FUN2'},'-p',{[],'Minima+Maxima'}) - makes a matrix A with the output FUN1 and
%      FUN2 using the last paranal results. Show only the Minima+Maxima of FUN2.
%
%
%   See also out, time, paranal, conteq
%
%   Reference page in Help browser:
%      <a href="matlab:commands('outfun')">commands outfun</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [F,ndxs] = outfun(varargin)
%(s, opt, acatfun)
global g_grind g_paranal g_permanent g_cont g_mcarlo g_Y g_t;
fieldnams={'fun', 'U1', 'The equation(s) to show';...
   'datastruc', 'r', 'You can use a data structure like g_paranal.run to extract the data';...
   'catfun', 'r#s#n', 'Get categorial statistics of the state variable';...
   'times', 'n', 'The times for outputs';...
   'runfield', 's', 'The name of the field in the datastruc that contains the run information (e.g. ''run'' or ''prevrun'')'}';
args=i_parseargs(fieldnams,'fun,catfun','-#ndxofparanal,-n,-p,-b,-c,-m,-t',varargin,false,{@i_is_outfun_equation});
if isfield(args,'catfun')&&isstruct(args.catfun)
    args.datastruc=args.catfun;
    args=rmfield(args,'catfun');
end
if ~isfield(args,'runfield')
    args.runfield='run';
end
if nargin<=1||numel(args.opts)>0
    i_parcheck;
end
ndxs = [];
if nargin == 0
    if ~isfield(g_grind, 'outfun')
        g_grind.outfun.fun = '';
        g_grind.outfun.opt = 1;
    end
    res = i_outfundlg( g_grind.outfun);
    if isempty(res)
        return;
    end
    if ischar(res.fun)&&~isempty(strfind(res.fun,' '))
       %if spaces split in different variables/functions
       res.fun=regexp(res.fun,'([''][^'']*[''])|([^ ]*)','match');
       res.fun=regexp(res.fun,'[^'']*','match','once'); 
    end
    args.fun = res.fun;
    if ~isempty(res.optstr)
        opt = ['-', res.optstr(1)];
        args.opts={opt};
    end
    g_grind.outfun.fun = res.fun;
    g_grind.outfun.opt = res.opt;
end

if ~isfield(args,'catfun')
    args.catfun = [];
end
if ~isfield(args,'fun')||isempty(args.fun)
    F = [];
    return;
end
if iscell(args.fun)&&numel(args.fun)==1
    args.fun=args.fun{1};
end
if iscell(args.fun) &&~isfield(args,'datastruc')%Multiple inputs
    Fs = cell(size(args.fun));
    Ndxs = cell(size(args.fun));
    maxn = 0;
    ny = 0;
    args1=args;
    for i = 1:length(args.fun)
        args1.fun=args.fun{i};
        if iscell(args.catfun)
            args1.catfun=args.catfun{i};
        end
        [Fs{i}, Ndxs{i}] = outfun(args1);
        if max(Ndxs{i}) > maxn
            maxn = max(Ndxs{i});
        end
        ny = ny + size(Fs{i}, 2);
    end
    F = zeros(maxn, ny) + NaN;
    k1 = 1;
    k2 = 0;
    if isempty(args.catfun) %significantly faster
        for i = 1:length(args.fun)
            if ~isempty(Fs{i})
                k2 = k2 + size(Fs{i}, 2);
                F(:, k1:k2) = Fs{i};
                k1 = k2 + 1;
            end
        end
        ndxs = Ndxs{1};
    else
        for i = 1:length(args.fun)
            k2 = k2 + size(Fs{i}, 2);
            F(Ndxs{i}, k1:k2) = Fs{i};
            k1 = k2 + 1;
        end
        ndxs = transpose((1:size(F, 1)));
    end
    return;
end
%ps=[];
if isfield(args,'datastruc')
    F=[];
    if iscell(args.fun)
        s=args;
        for i=1:length(args.fun)
           s.fun=args.fun{i};
           F = [F,outfun(s)];
        end
        return;
    else
        F = getoutfun(args.fun, args.datastruc.t, args.datastruc.Y, args.datastruc.perm, args.datastruc.parvalues, args.datastruc.pars);
    end
    if ~isempty(args.datastruc.parvalues)
        p  = transpose(repmat(args.datastruc.parvalues(:, 1), 1,size(args.datastruc.Y,1)));
    else
        p = [];
    end
else
    if isempty(args.opts)
        if isfield(args,'times')
            args.opts={'-t'};
        else
            args.opts={'-n'};
        end
    end
    anatype=find(strcmp(args.opts,{'-n','-p','-c','-b','-m','-t','-#ndxofparanal'}));
    if anatype == 1
        F = getoutfun(args.fun, g_t, g_Y, [], [], {});
        p = [];
    else
        oldY = g_Y;
        oldt = g_t;
        oldperm = g_permanent;
        switch anatype
            case 2 %-p
                if ~isfield(g_paranal, args.runfield)||isempty(g_paranal.(args.runfield))
                    error('GRIND:outfun:paranal','Error: for outfun(''-paranal'') it is needed that <a href="matlab:paranal">paranal</a> has been run');
                end
                if isfield(g_paranal.(args.runfield), 'perm')
                    F = getoutfun(args.fun, g_paranal.(args.runfield).t, g_paranal.(args.runfield).Y, g_paranal.(args.runfield).perm, g_paranal.(args.runfield).parvalues, g_paranal.(args.runfield).pars);
                else
                    F = getoutfun(args.fun, g_paranal.(args.runfield).t, g_paranal.(args.runfield).Y, [], g_paranal.(args.runfield).parvalues, g_paranal.(args.runfield).pars);
                end
                p  = transpose(repmat(g_paranal.(args.runfield).parvalues(:, 1), 1,size(g_paranal.(args.runfield).Y,1)));
            case 3 %-c
                if isempty(g_cont)||isempty(g_cont.(args.runfield).pars)
                    error('GRIND:outfun:conteq','Error: for outfun(''-conteq'') it is needed that <a href="matlab:conteq">conteq</a> has been run');
                end
                F = getoutfun(args.fun, 0, g_cont.(args.runfield).Y, [], g_cont.(args.runfield).parvalues, g_cont.(args.runfield).pars);
                p = g_cont.(args.runfield).parvalues;
            case 4 %-b  contbif
                if isempty(g_cont)||isempty(g_cont.contbifrun.pars)
                    error('GRIND:outfun:contbif','Error: for outfun(''-b'') it is needed that <a href="matlab:conteq">contbif</a> has been run');
                end
                F = getoutfun(args.fun, 0, g_cont.contbifrun.Y, [], g_cont.contbifrun.parvalues, g_cont.contbifrun.pars);
                p = g_cont.contbifrun.parvalues;
            case 5  %-m
                if isempty(g_mcarlo) ||~isfield(g_mcarlo, args.runfield)||isempty(g_mcarlo.(args.runfield))
                    error('GRIND:outfun:mcarlo','Error: for outfun(''-mcarlo'') it is needed that <a href="matlab:mcarlo">mcarlo</a> has been run');
                end
                F = getoutfun(args.fun, g_mcarlo.(args.runfield).t, g_mcarlo.(args.runfield).Y, g_mcarlo.(args.runfield).perm, g_mcarlo.(args.runfield).parvalues, g_mcarlo.(args.runfield).pars);
                p = [];
            case 6  %-t
                if ~isfield(args,'times')
                    if isfield(args,'catfun') %this is not really nice
                        args.times=str2num(args.catfun);
                        if isempty(args.times)
                            args.times=1;
                        end
                    else
                        args.times=1;
                    end
                end
                %not only one value, but also the correct shape of matrix
                F =  analysesingle(args.fun, args.times(1));
                if numel(args.times)>1 %typically the number of times should be small, otherwise interp1 is more efficient
                    F1=F;
                    F=zeros(size(F,1),size(F,2),length(args.times));
                    F(:,:,1)=F1;
                    for i=2:length(args.times)
                        F(:,:,i)=analysesingle(args.fun, args.times(i));
                    end
                end
                return;
            case 7
                if strcmp('-#ndxofparanal',args.opts) %used internally
                    g_t = g_paranal.(args.runfield).t(:);
                    g_Y = transpose(reshape(permute(g_paranal.(args.runfield).Y, [2, 1, 3]), size(g_paranal.(args.runfield).Y, 2), size(g_paranal.(args.runfield).Y, 1) * size(g_paranal.(args.runfield).Y, 3)));
                    if ~isempty(g_paranal.(args.runfield).perm)
                        g_permanent.Y = transpose(reshape(permute(g_paranal.(args.runfield).perm, [2, 1, 3]), size(g_paranal.(args.runfield).perm, 2), size(g_paranal.(args.runfield).perm, 1) * size(g_paranal.(args.runfield).perm, 3)));
                        g_permanent.t = g_t;
                    end
                    if nargin < 3
                        ndx = 1;
                    else
                        ndx = i_checkstr(args.catfun);
                    end
                    %not only one value, but also the correct shape of matrix
                    F =  analyseparanalsingle(args.fun, ndx, args.runfield);
                    return;
                end
            otherwise
                error('grind:outfun','unknown option');
        end
        g_Y = oldY;
        g_t = oldt;
        g_permanent = oldperm;
    end
end
if ~isempty(args.catfun)
    [F,ndxs] = catfun(args.catfun, F, p(:));
end

if (nargout > 1) && isempty(ndxs)
    ndxs = transpose(1:size(F, 1));
end

function F = getoutfun(s, t, Y, perm, parvalues, pars)
global g_Y g_permanent g_t g_data g_grind;
if isempty(Y)
    F = [];
    return;
end
iX = i_getno(s);
if ~isempty(iX.no)&&~iX.isfun
    F = getvariable(iX, t, Y, perm, parvalues, pars);
elseif strcontains(s, '<param')
    for ipar = 1:length(pars)
        s = strrep(s, sprintf('<param%d>', ipar), pars{ipar});
    end
    iX = i_getno(s);
    if ~isempty(iX.no)
        F = getvariable(iX, t, Y, perm, parvalues, pars);
    else
        F = getoutfun(s, t, Y, perm, parvalues, pars);
    end
elseif strncmpi(s, 'Observed ',9) %obsolete better to use observed('A')
    if ~isempty(g_data)
        ivar2 = strtrim(s(10:end));
        if strcmp(ivar2, 't')
            F = g_data.t;
        else
            indx =  strcmp(g_data.varlist, ivar2);
            F = g_data.obs(:, indx);
        end
    else
        F = [nan; nan];
    end
    return;
else
    s = outf('changeshortcut', s);
    if g_grind.statevars.vector
        s=changevect_varcol(s);
    end
    obj = parsed_equation(s);
    vars = symvar(obj);
    if any(obj.types==obj.vr.funcname)&&numel(Y)~=numel(g_Y) %some functions (for instance outf) are dependent on g_Y
        g_t = t(:);
        g_Y = transpose(reshape(permute(Y, [2, 1, 3]), size(Y, 2), size(Y, 1) * size(Y, 3)));
        if ~isempty(perm)
            g_permanent.Y = transpose(reshape(permute(perm, [2, 1, 3]), size(perm, 2), size(perm, 1) * size(perm, 3)));
            g_permanent.t = g_t;
        end
    end
    for i = length(vars):-1:1
        if i == length(vars)
            iXs = i_getno(vars{i});
            iXs(i) = iXs(1);
        else
            iXs(i) = i_getno(vars{i});
        end
        if iXs(i).isfun
            F = getauxvar(s, vars,t, Y, perm, parvalues, pars);
            return;
        elseif isempty(iXs(i).no)&&~isempty(which(vars{i}))
            vars(i) = [];
            iXs(i) = [];
        elseif isempty(iXs(i).no)
            error('grind:outfun','Unknown variable: %s',vars{i})
        end
    end
    g_l_values123 = cell(length(vars), 1);
    maxsiz=1;
    for i = 1:length(vars)
        g_l_values123{i} = double(getvariable(iXs(i), t, Y, perm, parvalues, pars));
        if size(g_l_values123{i},2)>maxsiz
            maxsiz=size(g_l_values123{i},2);
        end
        obj = obj.changevar(vars{i}, sprintf('g_l_values123{%d}', i));
    end
%     if maxsiz>1
%         for i = 1:length(vars)
%             if size(g_l_values123{i},2)==1
%                 g_l_values123{i}=repmat(g_l_values123{i},1,maxsiz);
%             end
%         end
%     end
        
    % if ~g_grind.statevars.vector %for vector models the equations should already be vectorized
    obj.equation =  vectorize(obj);
    %   end
    %   s1 = sprintf('%s', obj.fs{:});
    try
        F = eval(obj.equation);
    catch err
        fprintf(2, 'Error in equation "%s"\n', s);
        rethrow(err);
    end
    if length(F) == 1
        F = zeros(numel(t), 1) + F;
    end
    % fun=str2func(sprintf('@(g_vars) %s',s1));
    % F=fun(g_vars);
end
function s=changevect_varcol(s)
global g_grind;
%get variables with an index for instance A(1) and replace them by
%g_valcol(A,1)
[f1,f2,vars1]=regexp(s,'(%.*)|(''[^'']*'')|([a-zA-Z_][a-zA-Z0-9_]*\([0-9,]*\))','start','end','match');
for j=length(f1):-1:1
    if ~any(s(f1(j))=='''%')
        iXs = i_getno(vars1{j});
        if iXs.isvar||iXs.ispar
            s1=vars1{j};
            ff1=strfind(s1,'(');
            if ~strcontains(s1,',')
                s=sprintf('%sg_valcol(%s,%s%s',s(1:f1(j)-1),s1(1:ff1-1),s1(ff1+1:end-1),s(f2(j):end));
            else
                ndx=str2num(s1(ff1+1:end-1)); %#ok<ST2NM>
                if iXs.isvar
                    siz=[g_grind.statevars.dims{iXs.vecno}.dim1 g_grind.statevars.dims{iXs.vecno}.dim2];
                else
                    siz=size(evalin('base',s1));
                end
                index1=sub2ind(siz,ndx(1),ndx(2));
                s=sprintf('%sg_valcol(%s,%d%s',s(1:f1(j)-1),s1(1:ff1-1),index1,s(f2(j):end));
            end
        end
    end
end


function v = g_valcol(var,i) 
v=var(:,i);
function F = analysesingle(s, at)
%currently only supported for state variables and permanent vars
global g_Y g_t g_grind;
varno = i_getno(s);
if varno.isvar
    it=find(g_t >= at, 1);
    if isempty(it)
        it = length(g_t);
    end
    if g_grind.statevars.vector&&isempty(varno.ndx)
        dims = g_grind.statevars.dims{varno.vecno};
        F = g_Y(it, dims.from:dims.to);
        F = reshape(F, dims.dim1, dims.dim2);
    elseif g_grind.statevars.vector
        fromno = g_grind.statevars.dims{varno.vecno}.from;
        F = g_Y(it, varno.ndx+fromno-1);
    else
        F = g_Y(it, varno.no);
    end
elseif varno.ispar
    F=evalin('base',s);
elseif varno.istime
    F=at;
elseif varno.isperm
    F = defpermanent('-p', at);
elseif varno.isext
    F = externvar(varno.no, g_grind.externvars{varno.no}.default, at);
else
    if g_grind.statevars.vector
        s=changevect_varcol(s);
        s=strrep(s,'g_valcol(','g_valcol1(');
    end
    %if there are no auxiliary variables it still can be fast
    %get all variables no functions
    allvar=unique(regexp(s,'([A-Za-z_][A-Za-z_0-9]*)(?![\(A-Za-z0-9_])','match'));
    g_l_value123=cell(size(allvar));
    auxvars=false;
    for i=1:length(allvar)
        varno=i_getno(allvar{i});
        if ~varno.isfun
            g_l_value123{i}=analysesingle(allvar{i},at);
        else
            auxvars=true;
            break
        end
    end
    if ~auxvars
        for i=1:length(allvar)
            s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',allvar{i}),sprintf('g_l_value123{%d}', i));
        end
        F=eval(s);
    else
        %*****VERY INEFFICIENT!!****
        it=find(g_t >= at, 1);
        if isempty(it)
            it = length(g_t);
        end
        
        res = outfun(s);
        F = res(it, :);
    end
end
function v = g_valcol1(var,i) 
v=var(i);
function F = analyseparanalsingle(s, ndx,runfield)
%currently only supported for state variables and permanent vars
global g_paranal g_grind;
if ndx > numel(g_paranal.(runfield).t)
    ndx = numel(g_paranal.(runfield).t);
end
if ndx < 1
    ndx = 1;
end
[tndx, ~, stepndx] = ind2sub(size(g_paranal.(runfield).t), ndx);
varno = i_getno(s);
if varno.isvar
    if g_grind.statevars.vector
        dims = g_grind.statevars.dims{varno.vecno};
        F = g_paranal.(runfield).Y(tndx, dims.from:dims.to,stepndx);
        F = reshape(F, dims.dim1, dims.dim2);
    else
        F = g_paranal.(runfield).Y(tndx, varno.no,stepndx);
    end
    return;
elseif varno.isperm
    permvar = g_grind.permanent{varno.no};
    F = g_paranal.(runfield).perm(tndx, permvar.from:permvar.to,stepndx);
    F = reshape(F, permvar.dims(1), permvar.dims(2));
elseif varno.isext
    error('grind:outfun','not yet implemented');
    %F=externvar(varno.no,g_grind.externvars{varno.no}.default,at);
elseif strcmp(s, 't')
    F = g_paranal.(runfield).t(tndx, 1, stepndx);
else
    %not efficient!! and the shape may be wrong
    res = outfun(s, '-p');
    F = res(ndx, :);
end
%get a simple variable
function F = getvariable(iX, t, Y, perm, parvalues, pars)
global g_grind;
ntime = size(Y, 1);
if iX.ispar
    s = g_grind.pars{iX.no};
    if ~isempty(iX.ndx)
        siz=evalin('base',sprintf('size(%s)',s));
        [i,j]=ind2sub(siz,iX.ndx);
        if siz(2)==1
            s=sprintf('%s(%d)',s,i);
        elseif siz(1)==1
            s=sprintf('%s(%d)',s,j);
        else
            s=sprintf('%s(%d,%d)',s,i,j);
        end
    end
    ipar = strcmp(s, pars);
    if any(ipar)
        F = transpose(repmat(parvalues(:, ipar), 1,ntime));
        F = F(:);
    else
        F = evalin('base', s);
        if ~isempty(iX.ndx)&&numel(F)>1
            F=F(iX.ndx);
        end
        F = repmat(transpose(F(:)),numel(t), 1);
    end
    return;
elseif iX.isvar
    if ~isempty(iX.vecno) &&  isempty(iX.ndx)
        F = Y(:, g_grind.statevars.dims{iX.vecno}.from:g_grind.statevars.dims{iX.vecno}.to,:);
        if size(F,3)>1
            F = transpose(reshape(permute(F, [2, 1, 3]), size(F, 2), size(F, 1) * size(F, 3)));
        end
    else
        F = Y(:, iX.no,:);
        F = F(:);
    end
elseif iX.istime
    F = t(:);
elseif iX.isperm
    if isempty(perm)
        F = defpermanent('-g', iX.no);
    else
        p =  g_grind.permanent{iX.no};
        F1 = perm(:, p.from:p.to,:);
        F = transpose(reshape(permute(F1, [2, 1, 3]), size(F1, 2), size(F1, 1) * size(F1, 3)));
    end
elseif iX.isext
    F = externvar(iX.no,str2num(g_grind.externvars{iX.no}.default), t(:));  %#ok<ST2NM>
end
if ~isempty(iX.transform)
    F=iX.transform.fun(F);
end

function g_res1 = getauxvar(g_afun, g_l_symvar_afun, g_t, g_Y, g_perm, g_parvalues, g_pars)
global g_grind;
g_l_ntime = size(g_Y, 1);
g_Y =  transpose(reshape(permute(g_Y, [2, 1, 3]), size(g_Y, 2), g_l_ntime * size(g_Y, 3)));
if ~isempty(g_perm)
    g_perm =  transpose(reshape(permute(g_perm, [2, 1, 3]), size(g_perm, 2), size(g_perm, 1) * size(g_perm, 3)));
end
g_t =  g_t(:);

if isempty(g_parvalues)
    if ~isempty(g_grind.pars)
        eval(i_globalstr(g_grind.pars));
    end
else
    for g_i = 1:length(g_grind.pars)
        g_l_f = strcmp(g_grind.pars{g_i}, g_pars);
        if ~any(g_l_f)
            g_l_p = evalin('base', char(g_grind.pars{g_i}));  %#ok<NASGU>
            eval(sprintf('%s=g_l_p;', g_grind.pars{g_i}));
        else
            g_l_p =  transpose(repmat(g_parvalues(:, g_l_f), 1,g_l_ntime));
            g_l_p = g_l_p(:);  %#ok<NASGU>
            eval(sprintf('%s=g_l_p;', g_grind.pars{g_i}));
        end
    end
end

[g_afun] = checkvecfuncs(g_afun);

if ~isempty(g_Y)
    if g_grind.statevars.vector
        for g_i = 1:length(g_grind.statevars.vectnames)
            eval(sprintf('%s = g_Y(:,%d:%d);',g_grind.statevars.vectnames{g_i}, ...
                g_grind.statevars.dims{g_i}.from, g_grind.statevars.dims{g_i}.to));
        end
    else
        for g_l_i = 1:g_grind.statevars.dim
            eval([char(i_statevars_names(g_l_i)) '= g_Y(:,' num2str(g_l_i) ');']);
        end
    end
end
if ~isempty(g_grind.externvars)
    for g_l_m = 1:length(g_grind.externvars)
        eval(sprintf('%s=externvar(%d, %s,g_t);',g_grind.externvars{g_l_m}.name, ...
            g_l_m, g_grind.externvars{g_l_m}.default));
    end
end
perms = {};
if ~isempty(g_grind.permanent)
    perms = cell(length(g_grind.permanent), 1);
    if isempty(g_perm)
        for g_l_m = 1:length(g_grind.permanent)
            eval(sprintf('%s=defpermanent(''-get'',%d);', g_grind.permanent{g_l_m}.name,g_l_m));
            perms{g_l_m} = g_grind.permanent{g_l_m}.name;
        end
    else
        for g_l_m = 1:length(g_grind.permanent)
            eval(sprintf('%s=g_perm(:,%d);', g_grind.permanent{g_l_m}.name,g_l_m));
            perms{g_l_m} = g_grind.permanent{g_l_m}.name;
        end
    end
end



t = g_t; 
%temporary variable t
%evaluate functions if necessary (it supports matrix notation)
if ~isempty(g_grind) && isfield(g_grind, 'funcnames')&&~isempty(g_grind.funcnames.names) ...
        && isoverlap(g_l_symvar_afun, g_grind.funcnames.names) && ~isoverlap(g_l_symvar_afun,perms)
    if g_grind.statevars.vector
        error('grind:outfun:auxvar','Error: auxilary variables are no longer supported for output \nof vector/matrix models, define them instead \nas permanent variables using "defpermanent"'); 
        %       %  i_evalfuncs; %update g_grind.funcnames
        %       global g_func;
        %       i_update_g_func;
        %       for g_l_i = 1:length(g_grind.funcnames.names)
        %          eval(sprintf('%s = g_func(:,%d:%d);',g_grind.funcnames.names{g_l_i}, g_grind.funcnames.dims{g_l_i}.from,g_grind.funcnames.dims{g_l_i}.to));
        %       end
    else
        funs=regexprep(g_grind.funcs,'\<dwiener(','dwiener(t,');
        i_evalarray(funs);
    end
end

g_l_s1 = g_afun;
g_l_i=strfind(g_l_s1, '''');
if ~isempty(g_l_i) && (length(g_l_i) == 1)
    g_N0 = zeros(size(g_Y, 1), g_grind.statevars.dim);
    for g_l_i = 1:g_grind.statevars.dim
        g_N0(:, g_l_i) = g_Y(:, g_l_i);
    end
    NRes=i_runsinglestep(1,g_N0,true);
    if ~g_grind.statevars.vector
        for g_l_i = 1:g_grind.statevars.dim
            eval([char(g_grind.statevars.names{g_l_i}) '_ddif = NRes(:,' num2str(g_l_i) ');']);
        end
    else
        error('GRIND:outfun:NotImplemented','_ddif not implemented')
    end
    g_l_s1=strrep(g_l_s1,'''','_ddif');
end

%set all operands to array operands
g_res1 = i_evalarray(g_l_s1);
if size(g_res1, 1) == 1
    g_res1 = g_res1 * ones(size(g_Y, 1), 1);
end
function g_res1 = observed(ivar2)  %#ok<DEFNU>
global g_data;
if ~isempty(g_data)
    if strcmp(ivar2, 't')
        g_res1 = g_data.t;
    else
        indx =  strcmp(g_data.varlist, ivar2);
        g_res1 = g_data.obs(:, indx);
    end
else
    g_res1 = [nan; nan];
end
% function s = vectorize(s)
% ss = strcmp(s, '*');
% s(ss) = {'.*'};
% ss = strcmp(s, '^');
% s(ss) = {'.^'};
% ss = strcmp(s, '/');
% s(ss) = {'./'};
% ss = strcmp(s, '&&');
% s(ss) = {'&'};
% ss = strcmp(s, '||');
% s(ss) = {'|'};

function [g_afun1, g_afun2] = checkvecfuncs(g_afun)
global g_grind
g_afun2 = g_afun;
g_afun1 = g_afun;
if g_grind.statevars.vector
    i = length(g_afun1);
    fr = [];
    fc = [];
    fcom = [];
    while i > 0
        if g_afun1(i)==''''
            i = i - 1;
            while (i>0) && g_afun1(i)~=''''
                i = i - 1;
            end
        end
        if g_afun1(i) == ')'
            fr = i;
        end
        if g_afun1(i) == ':'
            fc = i;
        end
        if g_afun1(i) == ','
            fcom = i;
        end
        if (g_afun1(i) == '(')
            fl = i;
            i = i - 1;
            while (i>0)&&(isletter(g_afun1(i))||((g_afun1(i)>='0')&&(g_afun1(i)<='9'))||(g_afun1(i)=='_'))
                i = i - 1;
            end
            i = i + 1;
            if i < fl
                inam = i_getno(g_afun1(i:fr));
                if ~isempty(inam.no)
                    if (isempty(fc)||(fc > fr))&&~inam.ispar
                        if ~isempty(fcom)&&(fcom < fr)
                            g_afun1 = sprintf('%s(:,%d)%s',g_afun1(1:fl - 1),inam.ndx,g_afun1(fr + 1:end));
                        else
                            g_afun1 = [g_afun1(1:fl) ':,' g_afun1(fl + 1:end)];
                        end
                    end
                    g_afun2 = [g_afun2(1:fl - 1) g_afun2(fr + 1:end)];
                end
            end
        end
        i = i - 1;
    end
end

function found = isoverlap(list1, list2)
%found=isempty(intersect(list1,list2));
found = 0;
for i = 1:length(list1)
    %   for j = 1:length(list2)
    %       if strcmp(list1{i}, list2{j})
    if any(strcmp(list1{i}, list2))
        found = 1;
        return;
    end
    %   end
end

function [validated_x,errmsg]=i_is_outfun_equation(x)
%just the the default check for an equation, but if the equation contains
%short cuts (as we can use in out) they are handled
errmsg='';
if nargin==0
    validated_x='equation';
    return;
end
if iscell(x)
    for i=1:length(x)
        [~,errmsg]=i_is_outfun_equation(x{i});
        if ~isempty(errmsg)
            return
        end
    end
    validated_x=x;
    return;
end
x1=regexp(x,'<param[0-9]+>','match');
if isempty(x)||~isempty(x1)
    validated_x=x1;
    errmsg='';
else
    [validated_x,errmsg]=i_validat(x,'q',{});
end

