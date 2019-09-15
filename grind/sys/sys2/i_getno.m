function iX = i_getno(avar,maytransform)
global g_grind;
if nargin<2
    maytransform=false; %transform is relatively slow
end
iX =struct('isvar', 0, 'ispar',0,'isfun', 0,'istime',0,'isext',0,'isperm',0, 'ndx', [],'no', [], 'vecno', [],'transform',[]);
if ischar(avar)
    f =  strfind(avar, '(');
    if ~isempty(f)
        f1=strfind(avar,')');
        vecvar = [avar(1:f(1) - 1) avar(f1+1:end)];
    else
        vecvar=avar;
    end
    
    [iX.no, iX.vecno] = varno(avar,vecvar);
    iX.isvar = ~isempty(iX.no);
    if ~iX.isvar
        iX.no = i_parno(vecvar);
        if ~isempty(iX.no)
            iX.ispar = 1;
        else
            if strcmp(avar,'t')
                iX.istime=1;
                iX.no=1;
                return;
            end
            %   if ~isfield(g_grind.funcnames, 'dims')
            %      i_evalfuncs;
            %   end
            
            iX.no = i_permno(vecvar);
            iX.isperm = ~isempty(iX.no);
            if ~iX.isperm
                [iX.no, iX.vecno] = i_funno(vecvar);
                iX.isfun = ~isempty(iX.no);
                if ~iX.isfun
                    iX.no = i_externno(vecvar);
                    iX.isext = ~isempty(iX.no);
                end
                
            end
            
        end
        
    end
    
    if ~iX.isfun
        %      s = '[]';
        if ~isempty(f)
            f2 = strfind(avar, ')');
            %      f1 = strfind(avar, ',');
            %        f3 = strfind(avar, ':');
            %would we need something like this?
            %         if isempty(f3)
            %             if ~isempty(f1)
            %                f3 = [f1(1), f2(1)];
            %             else
            %                f3 = f2(1);
            %             end
            
            %          end
            
            iX.ndx=sscanf(avar(f + 1:f2 - 1),'%d,%d');
            if length(iX.ndx)==2
                if iX.isvar
                    iX.ndx=sub2ind([g_grind.statevars.dims{iX.vecno}.dim1,g_grind.statevars.dims{iX.vecno}.dim2],iX.ndx(1),iX.ndx(2));
                elseif iX.ispar || iX.isperm || iX.isext
                    p=evalin('base',avar(1:f - 1));
                    iX.ndx=sub2ind(size(p),iX.ndx(1),iX.ndx(2));
                else
                    if ~isfield(g_grind.funcnames, 'dims')
                        iX.ndx = NaN;
                    elseif ~isempty(iX.no)
                        iX.ndx = sub2ind([g_grind.funcnames.dims{iX.no}.dim1,g_grind.funcnames.dims{iX.no}.dim2],iX.ndx(1),iX.ndx(2));
                    end
                    
                end
                
            end
            
        end
        
        if iX.isvar&&~isempty(iX.ndx)&&g_grind.statevars.vector
            d=g_grind.statevars.dims{iX.vecno};
            iX.no=d.from-1+iX.ndx;
        elseif iX.isfun
            iX.no = iX.no + iX.ndx - 1;
        elseif maytransform
            %is it an transformed variable?
            p=parsed_equation(avar);
            if length(p.fs)>1  %there is an equation
                vars=symvar(p);
                if length(vars)==1 %only one variabele is allowed for transformation
                    thevar=vars{1};
                    iX1=i_getno(thevar);
                    if iX1.no>0
                        iX1.transform.fun=str2func(sprintf('@(%s)%s',thevar,avar));
                        if length(p.fs)==4&&any(strcmp(p.fs{1},{'imag','real'}))
                            if strcmp(p.fs{1},'imag')
                                iX1.transform.fun=@imag;
                                iX1.transform.invfun=@(x,prevx)complex(real(prevx),x);
                            else
                                iX1.transform.fun=@real;
                                iX1.transform.invfun=@(x,prevx)complex(x,imag(prevx));
                            end
                        elseif strcmp(p.fs{1},'mod')
                            [pars,j]=p.getfunctionpars(1);
                            if length(p.fs)==j(end)+1&&strcmp(pars{1},thevar)
                                iX1.transform.invfun=str2func(sprintf('@(%s)%s',thevar,avar));
                            end
                        elseif i_hastoolbox('symbolic')
                            %use symbolic toolbox to determine the inverse
                            %function
                            try
                            if strcmp(thevar,'x')
                                invfun=sprintf('@(x1)%s',char(solve(sym(['x1=' avar]),sym(thevar))));
                            else
                                invfun=sprintf('@(x)%s',char(solve(sym(['x=' avar]),sym(thevar))));
                            end
                            iX1.transform.invfun=str2func(invfun);
                            catch
                            end
                        else
                            %this is not always working
                            funstruc=fliplr(p.structure);
                            newstruc=funstruc;
                            for i=1:length(newstruc)
                                if strcmp(funstruc(i).args{1},thevar)
                                    newstruc(i).leftvar='ans';
                                else
                                    newstruc(i).leftvar=funstruc(i).args{1};
                                end
                                if strcmp(funstruc(i).leftvar,'ans')
                                    newstruc(i).args{1}='x';
                                else
                                    newstruc(i).args{1}=funstruc(i).leftvar;
                                end
                                switch funstruc(i).oper
                                    case '+'
                                        newstruc(i).oper='-';
                                    case '.*'
                                        newstruc(i).oper='./';
                                    case '*'
                                        newstruc(i).oper='./';
                                    case '-'
                                        newstruc(i).oper='+';
                                    case './'
                                        newstruc(i).oper='.*';
                                    case '/'
                                        newstruc(i).oper='.*';
                                    case 'exp'
                                        newstruc(i).oper='log';
                                    case 'log'
                                        newstruc(i).oper='exp';
                                    case 'ln'
                                        newstruc(i).oper='exp';
                                    case 'sin'
                                        newstruc(i).oper='asin';
                                    case 'cos'
                                        newstruc(i).oper='acos';
                                    case 'tan'
                                        newstruc(i).oper='atan';
                                    otherwise
                                        return;
                                end
                            end
                            invfun=char(parsed_equation(newstruc));
                            iX1.transform.invfun=str2func(sprintf('@(x)%s',invfun));
                        end
                    end
                    
                    aval=rand(1);
                    try
                        if nargin(iX1.transform.invfun)==2
                           tval=iX1.transform.invfun(iX1.transform.fun(aval),aval);
                        else
                           tval=iX1.transform.invfun(iX1.transform.fun(aval));
                        end
                        if abs(tval-aval)<1e-10
                            iX=iX1;
                        end
                    catch %#ok<CTCH>
                    end
                end
            end
        end
        
    end
    
end

%i_varno returns the number of a statevar and [] if no statevar
function [varno, vecno] = varno(~,vecvar)
global g_grind;
varno=[];
vecno=[];
if g_grind.statevars.vector
    vecno = find(strcmp(vecvar, g_grind.statevars.vectnames));
    if ~isempty(vecno)
        varno =g_grind.statevars.dims{vecno}.from;
    end
    
else
    varno=find(strcmp(vecvar, g_grind.statevars.names));
end

function permno = i_permno(apar)
global g_grind;
permno=[];
if isfield(g_grind, 'permanent')
    for k = 1:size(g_grind.permanent, 2)
        if strcmp(apar, char(g_grind.permanent{k}.name))
            permno = k;
            return;
        end
    end
end

function extno = i_externno(apar)
global g_grind;
extno = [];
for k = 1:size(g_grind.externvars, 2)
    if strcmp(apar, char(g_grind.externvars{k}.name))
        extno = k;
        return;
    end
end

function parno = i_parno(apar)
global g_grind;
parno=find(strcmp(apar,g_grind.pars));
%
% for k = 1:length(g_grind.pars)
%    if strcmp(apar, char(g_grind.pars{k}))
%       parno = k;
%       return;
%    end
% end
function [no, funno] = i_funno(apar)
global g_grind;
funno = [];
no = [];
f =  strfind(apar, '(');
f2 = strfind(apar, ')');
if (length(f)==1) && (length(f2)==1) && (f2==length(apar))
    apar = apar(1:f(1) - 1);
elseif ~isempty(f)
    return;
end

funno=find(strcmp(apar,g_grind.funcnames.names));
if ~isempty(funno)
    if ~isfield(g_grind.funcnames, 'dims')
        no = NaN;
    else
        no = g_grind.funcnames.dims{funno}.from;
    end
end


% for k = 1:size(g_grind.funcnames.names, 2)
%    if strcmp(apar, char(g_grind.funcnames.names{k}))
%       funno = k;
%       if ~isfield(g_grind.funcnames, 'dims')
%          no = NaN;
%       else
%          no = g_grind.funcnames.dims{k}.from;
%       end

%       return;
%    end
%end
