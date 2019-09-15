function h = i_getodehandle(type, parslist, opt)
% 0=normal = normal run (incl backw) i_ru
% 1=singlestep = singlestep  +1
% -1=singleback = singlestep   -1
% 2=potential = potential
% 3=makeimplicit = make implicit
% 5=findeq = findeq with fminsearch
% 6=findeq_fsolve = findeq with fsolve
% 7=makemap = maps of makemap
% 8=timesens = test symbolic timesens one parameter needed
% 9=jac = Jacobian (t, x)
% matcont = matcont ode file
% 10=coco = coco ode file {nonautonomic not yet supported)
% 11=coco_jac = Jacobian (x, p)(coco) [parslist must be cell]
% 12=coco_jacp = Jacobianp (x, p) (coco) [parslist must be cell]
% 13=hess = Hessian (t, x) not vector ;pars=vector of parslist
% 14=coco_lyap-First lyapunov in COCO
%15 quasipot - quasipotential function
%16 delays 
global g_grind;
if nargin == 1
    parslist = '';
end

if isempty(g_grind.pars)
    allpars = '';
else
    allpars = sprintf(',%s',g_grind.pars{:});
end

h = [];
[parslist,allpars]=checkpars(parslist,allpars);
switch type
    case {0,'normal'} %normal i_ru
        %Note that i_differ handles the iters
        if ~g_grind.solver.backwards
            if g_grind.solver.isimplicit
                h=evalin('base',sprintf('@(t,g_X1,g_X3)%s(t,g_X1,g_X3%s)',g_grind.odefile,allpars));
            elseif g_grind.solver.haslags
                h=evalin('base',sprintf('@(t,g_X1,g_lags%s)%s(t,g_X1,g_lags%s)',parslist,g_grind.odefile,allpars));
            else
                h=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1%s)',g_grind.odefile,allpars));
            end

        else
            if g_grind.solver.isimplicit
                h=evalin('base',sprintf('@(t,g_X1,g_X3)-%s(t,g_X1,-g_X3%s)',g_grind.odefile,allpars));
            elseif g_grind.solver.isdiffer
                h=i_getodehandle('singleback');
            elseif g_grind.solver.haslags
                h=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1,[]%s)',g_grind.odefile,allpars));
            else
                h=evalin('base',sprintf('@(t,g_X1)-%s(t,g_X1%s)',g_grind.odefile,allpars));
            end

        end

        if ~isempty(g_grind.solver.opt.Jacobian)
            enterjac('-update');
        end

    case {1,'singlestep'} %singlestep  + 1
        if g_grind.solver.isimplicit
            evalin('base',sprintf('g_grind.solver.h=@(t,g_X1,g_X3%s)%s(t,g_X1,g_X3%s);',parslist,g_grind.odefile,allpars));
            h = @implicitstep;
        elseif g_grind.solver.haslags
            h=evalin('base',sprintf('@(t,g_X1%s)%s(t,g_X1,[]%s)',parslist,g_grind.odefile,allpars));
        elseif g_grind.solver.iters > 1
            g_hsolver=evalin('base',sprintf('@(t,g_X1%s)%s(t,g_X1%s);',parslist,g_grind.odefile,allpars));  %#ok<NASGU>
            h = eval(sprintf('@(t,g_X1%s)iterstep(t,g_X1,g_hsolver,g_grind.solver.iters%s)',parslist,parslist));
        else
            h=evalin('base',sprintf('@(t,g_X1%s)%s(t,g_X1%s)',parslist,g_grind.odefile,allpars));
        end

    case   {-1,'singleback'} %singlestep   - 1
        if g_grind.solver.isimplicit
            h=evalin('base',sprintf('@(t,g_X1,g_X3%s)%s(t,g_X1,-g_X3%s)',g_grind.odefile,allpars));
        elseif g_grind.solver.isdiffer
            TOL=1E-10;
            MAXITS=100;
            odef=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1%s)',g_grind.odefile,allpars));
            Jfun= g_grind.solver.opt.Jacobian;
            if g_grind.solver.iscomplex
                 error('grind:backw:complex','Backwards not supported for complex difference equations');
           else
                h = @(t,y)newton4differ(t,y,odef,Jfun,TOL,MAXITS,g_grind.solver.iters);
            end
        else
            h=evalin('base',sprintf('@(t,g_X1)-%s(t,g_X1%s)',g_grind.odefile,allpars));
        end

    case {2,'potential'} %potential
        h=evalin('base',sprintf('@(g_X1,t%s)-%s(1 + 1E-5,g_X1%s)',parslist,g_grind.odefile,allpars));
    case {3,'makeimplicit'} %make implicit
        g_odefun=i_getodehandle(1,parslist);  %#ok<NASGU>
        h=eval(sprintf('@(t,g_X1,g_X3%s)makeimplicit(t,g_X1,g_X3,g_odefun%s)',parslist,parslist));
        % h=evalin('base',sprintf('@(t,g_X1,g_X3%s)g_X3-%s(t,g_X1%s)',parslist,g_grind.odefile,allpars));
    case {5,'findeq'} %findeq with fminsearch
        if ~isfield(g_grind,'findeq')||isempty(g_grind.findeq.mask)
            if g_grind.solver.isimplicit
                evalin('base',sprintf('g_grind.solver.h=@(t,g_X1,g_X3)%s(t,g_X1,g_X3%s);',g_grind.odefile,allpars));
                h = @(x)norm(implicitstep(1, x));
            elseif g_grind.solver.haslags
                h=evalin('base',sprintf('@(g_X1)norm(%s(t,g_X1,[]%s))',g_grind.odefile,allpars));
            elseif g_grind.solver.isdiffer
                h1=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1%s);',g_grind.odefile,allpars));
                h = @(x)norm(iterstep(1, x, h1, g_grind.solver.iters)-x);
            else
                if isfield(g_grind,'findeq')&&~isempty(g_grind.findeq.userfun)
                   h=evalin('base',sprintf('@(g_X1)norm([%s(1,g_X1%s);g_grind.findeq.userfun(0,g_X1)])',g_grind.odefile,allpars));
                else
                   h=evalin('base',sprintf('@(g_X1)norm(%s(1,g_X1%s))',g_grind.odefile,allpars));
                end
            end

        else
            g_grind.solver.h = i_getodehandle(1, '');
            h = @(x)norm(findeqmask(x));
        end
    case 'findeq_complex' %findeq with fminsearch and complex values
        if ~isfield(g_grind,'findeq')||isempty(g_grind.findeq.mask)
            if g_grind.solver.isdiffer
                h1=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1%s);',g_grind.odefile,allpars));
                h = @(x)norm(iterstep(1, complex(x(1:end/2,:),x(end/2+1:end)), h1, g_grind.solver.iters)-complex(x(1:end/2,:),x(end/2+1:end)));
            else
                h=evalin('base',sprintf('@(g_X1)norm(%s(1,complex(g_X1(1:end/2,:),g_X1(end/2+1:end,:))%s))',g_grind.odefile,allpars));
            end

        else
            g_grind.solver.h = i_getodehandle(1, '');
            h = @(x)norm(findeqmask(complex(x(1:end/2,:),x(end/2+1:end))));
        end

    case {6,'findeq_fsolve'} %findeq with fsolve
        if ~isfield(g_grind,'findeq')||isempty(g_grind.findeq.mask)
            if g_grind.solver.isimplicit
                evalin('base',sprintf('g_grind.solver.h=@(t,g_X1,g_X3)%s(t,g_X1,g_X3%s);',g_grind.odefile,allpars));
                h = @(x)implicitstep(1, x);
            elseif g_grind.solver.haslags
                h=evalin('base',sprintf('@(g_X1)%s(t,g_X1,[]%s)',g_grind.odefile,allpars));
            elseif g_grind.solver.iters > 1
                evalin('base',sprintf('g_grind.solver.h=@(t,g_X1)%s(t,g_X1%s);',g_grind.odefile,allpars));
                h = @(x)iterstep(1, x);
            else
                h=evalin('base',sprintf('@(g_X1)%s(1,g_X1%s)',g_grind.odefile,allpars));
            end

        else
            g_grind.solver.h = i_getodehandle(1, '');
            h = @findeqmask;
        end

    case {7,'makemap'} %maps of makemap
        if g_grind.solver.map.haslags
            h=evalin('base',sprintf('@(t,g_X1,g_lags)%s(t,g_X1,g_lags%s)',g_grind.solver.map.odefile,allpars));
        else
            h=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1%s)',g_grind.solver.map.odefile,allpars));
        end

    case {8,'timesens'} %test symbolic timesens one parameter needed
        g_odehandle = i_getodehandle(1, parslist);
        g_X1 = i_initvar;
        h = @(t, apar)g_odehandle(t, g_X1, apar);
        
    case 'Jacobianp' %Jacobianp (x, p)
        if nargin<3
            opt='';
        end

        if ~isempty(parslist)
            parslist1 = regexp(parslist,'[,]','split');
            ndx=cellfun('isempty',parslist1);
            parslist1=parslist1(~ndx);
        %    parslist1(end)=[];
        else
            parslist1={};
        end

        if ~isempty(g_grind.syms.Jacobianp)&&~strcmp(opt,'numonly')
            if ~g_grind.solver.isimplicit
                [ndx1,ndx2] = ismember(g_grind.pars, parslist1);
                ndx = find(ndx1); %some tricks to preserve order of pars
                if length(ndx) ~= length(parslist1)
                    unknownpars = setdiff(parslist1, g_grind.pars);
                    error('grind:i_getodehandle','Unknown or double parameter(s): %s',sprintf('%s ',unknownpars{:}))
                end

                ndx2 = ndx2(ndx1);
                ndx = ndx(ndx2);
                if isa(g_grind.syms.Jacobianp,'function_handle')
                s=func2str(g_grind.syms.Jacobianp);
                if strncmp(s,'@(',2)
                    h=evalin('base',s);
                else
                    s= sprintf('@(t,g_X1,%s)%s(1,g_X1%s)',parslist,s,allpars);
%                     for i = 1:length(parslist1)
%                         s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',parslist{i}),sprintf('g_par(%d,:)', i));
%                     end
                    h=evalin('base',s);
                    if numel(ndx)<numel(g_grind.pars)
                        %first calculate the whole jac then selecting the
                        %ndxed part
                        h=@(g_X1,varargin)vindex(h(g_X1,varargin{:}),ndx);
                    end
                end
                return
                end
                jac = g_grind.syms.Jacobianp(:, ndx); %select parameters
                s = sprintf('%s,', jac{:});
                if strcmp(g_grind.solver.opt.Vectorized, 'on')
                    if isempty(s)
                       s = sprintf('@(t,g_X1%s)i_reshapevcat([%d,%d,max(size(g_X1,2),size(t,2))]);',parslist,g_grind.statevars.dim,numel(parslist1));
                    else
                       s = sprintf('@(t,g_X1%s)i_reshapevcat([%d,%d,max(size(g_X1,2),size(t,2))],%s);',parslist,g_grind.statevars.dim,numel(parslist1),s(1:end - 1));
                    end
%                     for i = 1:length(parslist)
%                         s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',parslist{i}),sprintf('g_par(%d,:)', i));
%                     end
                    s = vectorize(s);
                else
                    s = sprintf('@(t,g_X1%s)reshape([%s],[%d,%d]);',parslist,s,g_grind.statevars.dim,length(parslist1));
%                     for i = 1:length(parslist)
%                         s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',parslist{i}),sprintf('g_par(%d)', i));
%                     end

                end
                h = evalin('base', s);
              %  g_odefun = i_getodehandle('coco', parslist);
              %  h=@(x,p)checknan(x,p,handle,@(x,p)coco_num_DFDP(g_odefun, x, p));
            end

        end    
        
        
    case {9,'Jacobian'} %Jacobian (t, x)
        if nargin<3
            opt='';
        end
        if ~g_grind.solver.isimplicit&&strcmp(solver('name'),'ode15i')
            h=[];
            return;
        end
        if ~isempty(parslist)
            parslist1 = regexp(parslist,'[,]','split');
            ndx=cellfun('isempty',parslist1);
            parslist1=parslist1(~ndx);
        %    parslist1(end)=[];
        else
            parslist1={};
        end

        if ~isempty(g_grind.syms.Jacobian)&&~strcmp(opt,'numonly')
            if isa(g_grind.syms.Jacobian,'function_handle')
                s=func2str(g_grind.syms.Jacobian);
                if strncmp(s,'@(',2)
                    h=evalin('base',s);
                else
                    h=evalin('base', sprintf('@(t,g_X1%s)%s(t,g_X1%s)',parslist,s,allpars));
                end
            else
                if ~g_grind.solver.isimplicit
                    s = sprintf('%s,', g_grind.syms.Jacobian{:});
                    if strcmp(g_grind.solver.opt.Vectorized, 'on')
                        if g_grind.solver.backwards
                            s = sprintf('@(t,g_X1%s)-i_reshapevcat([%d,%d,max(size(g_X1,2),size(t,2))],%s);',parslist,g_grind.statevars.dim,g_grind.statevars.dim,s(1:end - 1));
                        else
                            s = sprintf('@(t,g_X1%s)i_reshapevcat([%d,%d,max(size(g_X1,2),size(t,2))],%s);',parslist,g_grind.statevars.dim,g_grind.statevars.dim,s(1:end - 1));
                        end
                        s = vectorize(s);
                    else
                        if g_grind.solver.backwards
                            s= sprintf('@(t,g_X1%s)-reshape([%s],[%d,%d,max(size(t,2),size(g_X1,2))]);',parslist,s,g_grind.statevars.dim,g_grind.statevars.dim);
                        else
                            s= sprintf('@(t,g_X1%s)reshape([%s],[%d,%d,max(size(t,2),size(g_X1,2))]);',parslist,s,g_grind.statevars.dim,g_grind.statevars.dim);
                        end
                    end
                    h=evalin('base',s);
                else
                    h = @i_implicit_Jac;
                end
            end
        elseif ~g_grind.solver.isimplicit
%            h=[];
             g_odefun = i_getodehandle('singlestep', parslist);  %#ok<NASGU>
             if strcmp(g_grind.solver.opt.Vectorized, 'on')
                 h = eval(sprintf('@(t,g_X1%s)num_DFDXv(g_odefun, t, g_X1%s);',parslist,parslist));
             else
                 h = eval(sprintf('@(t,g_X1%s)num_DFDX(g_odefun, t, g_X1%s);',parslist,parslist));
             end
        else
            h=[];
        end

    case 'vectorcheck' %parameter vector with time (for i_parcheck)
        if isempty(g_grind.pars)
            parlist='';
        else
            g_par = str2cell(sprintf('g_par(%d,:)\n',1:length(g_grind.pars)));
            parlist = sprintf(',%s',g_par{:});
        end

        h=evalin('base',sprintf('@(t,g_X1,g_par)%s(0,g_X1%s)',g_grind.odefile,parlist));
    case 'matcont'
        if g_grind.statevars.vector
            newpars=pars2coco(parslist,true);
            parlist = sprintf(',%s',newpars{:});
            
            if nargin<3
                opt=g_grind.odefile;
            elseif isa(opt,'function_handle')
                opt=funct2str(opt);
            end
            h=evalin('base',sprintf('@(t,g_X1,varargin)%s(t,g_X1%s)',opt,parlist));
        else
            if isempty(parslist)
                parslist='';
            else
                parslist = sprintf(',%s',parslist{:});
            end
            h=i_getodehandle('singlestep',parslist);
        end
    case {10,'coco'} %coco parameter vector
        if isempty(parslist)
            parslist = g_grind.pars;
        end

        if nargin<3
            opt=g_grind.odefile;
        elseif isa(opt,'function_handle')
            opt=funct2str(opt);
        end

        newpars=pars2coco(parslist);
        parlist = sprintf(',%s',newpars{:});
        if g_grind.solver.isdiffer
            if g_grind.solver.iters==1
                h=evalin('base',sprintf('@(g_X1,g_par)%s(0,g_X1%s)-g_X1',opt,parlist));
            else
                error('not yet implemented');
            end
        elseif ~g_grind.solver.isimplicit
            h=evalin('base',sprintf('@(g_X1,g_par)%s(0,g_X1%s)',opt,parlist));
        end
     case {12,'coco_jacp'} %Jacobian (x, p)
        if isempty(parslist)
            parslist = g_grind.pars;
        end

        if nargin<3
            opt=g_grind.odefile;
        elseif isa(opt,'function_handle')
            opt=funct2str(opt);
        end

        newpars=pars2coco(parslist);
        ndx= find(~strncmp(newpars,'g_par(',6));
        for i=1:length(ndx)
            newpars{ndx(i)}=sprintf('g_parlist.%s',newpars{ndx(i)});
            %we need to have the values of the non-active parameters
        end
        parlist = sprintf(',%s',newpars{:});
        
        if ~isempty(g_grind.syms.Jacobianp)&&~strcmp(opt,'numonly')
            g_odefun=i_getodehandle('Jacobianp',sprintf(',%s',g_grind.pars{:})); %#ok<NASGU>
            g_parlist=par('-v'); %#ok<NASGU>
            g_handle=eval(sprintf('@(g_X1,g_par)g_odefun(0,g_X1%s)',parlist));
            g_odefun = i_getodehandle('coco', parslist);
            h=@(x,p)checknan(x,p,g_handle,@(x,p)coco_num_DFDP(g_odefun, x, p));
        else
             g_odefun = i_getodehandle('coco', parslist);
            if strcmp(g_grind.solver.opt.Vectorized, 'on')
               h = @(x,p)coco_num_DFDPv(g_odefun, x, p);
            else
               h = @(x,p)coco_num_DFDP(g_odefun, x, p);
            end
        end 
    case {11,'coco_jac'} %Jacobian (x, p)
        if isempty(parslist)
            parslist = g_grind.pars;
        end

        if nargin<3
            opt=g_grind.odefile;
        elseif isa(opt,'function_handle')
            opt=funct2str(opt);
        end

        newpars=pars2coco(parslist);
        ndx= find(~strncmp(newpars,'g_par(',6));
        for i=1:length(ndx)
            newpars{ndx(i)}=sprintf('g_parlist.%s',newpars{ndx(i)});
            %we need to have the values of the non-active parameters
        end
        parlist = sprintf(',%s',newpars{:});
        
        if ~isempty(g_grind.syms.Jacobian)&&~strcmp(opt,'numonly')
            g_odefun=i_getodehandle('Jacobian',sprintf(',%s',g_grind.pars{:})); %#ok<NASGU>
            g_parlist=par('-v'); %#ok<NASGU>
            h=eval(sprintf('@(g_X1,g_par)g_odefun(0,g_X1%s)',parlist));
        else
             g_odefun = i_getodehandle('coco', parslist);
            if strcmp(g_grind.solver.opt.Vectorized, 'on')
               h = @(x,p)coco_num_DFDXv(g_odefun, x, p);
            else
               h = @(x,p)coco_num_DFDX(g_odefun, x, p);
            end
        end          

    case 'coco_lyap' %coco first lyapunov coefficent for Sub/Sup Hopf
        if isempty(parslist)
            parslist = g_grind.pars;
        end

     %   activepars=sprintf(',%s',parslist{:});
        allpars=sprintf(',%s',g_grind.pars{:});
        newpars=pars2coco(parslist);
        parlist = sprintf(',%s',newpars{:});
        if ~g_grind.solver.isimplicit&&~g_grind.solver.isdiffer
            % [lyap,matcontlyap] = firstlyap(t,x,odefile,Jac,Hess,varargin)
            evalin('base',sprintf('g_odehandle=i_getodehandle(''singlestep'',''%s'');',allpars));
            evalin('base',sprintf('g_odehandlejac=i_getodehandle(''Jacobian'',''%s'');',allpars));
            evalin('base',sprintf('g_odehandlehess=i_getodehandle(''Hessian'',''%s'');',allpars));
            g_firstlyap=evalin('base',sprintf('@(g_X1,g_par)firstlyap(0,g_X1,g_odehandle,g_odehandlejac,g_odehandlehess%s)',parlist));
            evalin('base','clear var g_odehandle g_odehandlejac g_odehandlehess');
            h=@(prob,data,u)fun2cocofun(prob,data,u,g_firstlyap);
        else
            h=[];
         end

%     case {12,'coco_jacp'} %Jacobianp (x, p)
%         if nargin<3
%             opt='';
%         end
% 
%         if isempty(parslist)
%             parslist = g_grind.pars;
%         end
% 
%         if ~isempty(g_grind.syms.Jacobianp)&&~strcmp(opt,'numonly')
%             if ~g_grind.solver.isimplicit
%                 [ndx1,ndx2] = ismember(g_grind.pars, parslist);
%                 ndx = find(ndx1); %some tricks to preserve order of pars
%                 if length(ndx) ~= length(parslist)
%                     unknownpars = setdiff(parslist, g_grind.pars);
%                     error('grind:i_getodehandle','Unknown or double parameter(s): %s',sprintf('%s ',unknownpars{:}))
%                 end
% 
%                 ndx2 = ndx2(ndx1);
%                 ndx = ndx(ndx2);
%                 if isa(g_grind.syms.Jacobianp,'function_handle')
%                 s=func2str(g_grind.syms.Jacobianp);
%                 if strncmp(s,'@(',2)
%                     h=evalin('base',s);
%                 else
%                     s= sprintf('@(g_X1,g_par)%s(1,g_X1%s)',s,allpars);
%                     for i = 1:length(parslist)
%                         s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',parslist{i}),sprintf('g_par(%d,:)', i));
%                     end
%                     h=evalin('base',s);
%                     if numel(ndx)<numel(g_grind.pars)
%                         %first calculate the whole jac then selecting the
%                         %ndxed part
%                         h=@(g_X1,g_par)vindex(h(g_X1,g_par),ndx);
%                     end
%                 end
%                 return
%                 end
%                 jac = g_grind.syms.Jacobianp(:, ndx); %select parameters
%                 if strcmp(g_grind.solver.opt.Vectorized, 'on')
%                     s = sprintf('%s,', jac{:});
%                     s = sprintf('@(g_X1,g_par)i_reshapevcat([%d,%d,max(size(g_X1,2),size(g_par,2))],%s);',g_grind.statevars.dim,length(parslist),s(1:end - 1));
%                     for i = 1:length(parslist)
%                         s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',parslist{i}),sprintf('g_par(%d,:)', i));
%                     end
% 
%                     s = vectorize(s);
%                 else
%                     s = sprintf('%s;', jac{:});
%                     s = sprintf('@(g_X1,g_par)reshape([%s],[%d,%d]);',s,g_grind.statevars.dim,length(parslist));
%                     for i = 1:length(parslist)
%                         s=regexprep(s,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',parslist{i}),sprintf('g_par(%d)', i));
%                     end
% 
%                 end
% 
%                 g_handle = evalin('base', s);
%                 g_odefun = i_getodehandle('coco', parslist);
%                 h=@(x,p)checknan(x,p,g_handle,@(x,p)coco_num_DFDP(g_odefun, x, p));
%             end
% 
%         else
%             g_odefun = i_getodehandle('coco', parslist);
%             if strcmp(g_grind.solver.opt.Vectorized, 'on')
%                 h = @(x,p)coco_num_DFDPv(g_odefun, x, p);
%             else
%                 h = @(x,p)coco_num_DFDP(g_odefun, x, p);
%             end
% 
%         end
    case {'jumps'}
        h=[];
        if i_split_sde&&isfield(g_grind.solver,'djump')
            s=[];
            for j=1:size(g_grind.solver.djump.jumps,1)
                s1='';
                for k=1:size(g_grind.solver.djump.jumps,2)
                    if g_grind.statevars.vector&&g_grind.statevars.dims{j}.dim2>1
                        s1=sprintf('%sreshape(%s,[%d,1]),',s1,g_grind.solver.djump.jumps{j,k},g_grind.statevars.dims{j}.dim1*g_grind.statevars.dims{j}.dim2);
                    else
                        s1=sprintf('%s%s,',s1,g_grind.solver.djump.jumps{j,k});
                    end
                end
                if ~g_grind.statevars.vector
                    s1(end)=';';
                end
                s=[s s1];
            end
            if g_grind.statevars.vector
                s=s(1:end-1);
                s=sprintf('@(t,g_X1%s)reshape(vertcat(%s),[%d,%d]).''',parslist,s,size(g_grind.solver.djump.jumps,2),g_grind.statevars.dim);
            else
                s=sprintf('@(t,g_X1%s)[%s]',parslist,s);
            end
            h=evalin('base',s);
        end
    case  {'drift','diffusion'}
        h=[];
        if i_split_sde&&isfield(g_grind.solver,'dwiener')
            s=[];
            for j=1:size(g_grind.syms.(type),1)
                s1='';
                for k=1:size(g_grind.syms.(type),2)
                    if g_grind.statevars.vector&&g_grind.statevars.dims{j}.dim2>1
                        s1=sprintf('%sreshape(%s,[%d,1]),',s1,g_grind.syms.(type){j,k},g_grind.statevars.dims{j}.dim1*g_grind.statevars.dims{j}.dim2);
                    else
                        s1=sprintf('%s%s,',s1,g_grind.syms.(type){j,k});
                    end
                end
                if ~g_grind.statevars.vector
                    s1(end)=';';
                end
                s=[s s1];
            end
            if g_grind.statevars.vector
                s=s(1:end-1);
                s=sprintf('@(t,g_X1%s)reshape(vertcat(%s),[%d,%d]).''',parslist,s,size(g_grind.syms.(type),2),g_grind.statevars.dim);
            else
                s=sprintf('@(t,g_X1%s)[%s]',parslist,s);
            end
            h=evalin('base',s);
        end
%     case {'drift','diffusion'}
%         if i_split_sde
%             s=[];
%             for j=1:size(g_grind.syms.(type),1)
%                 s1=sprintf('%s,',g_grind.syms.(type){j,:});
%                 if ~g_grind.statevars.vector
%                    s1(end)=';';
%                 end;
%                 s=[s s1];
%             end
%             if g_grind.statevars.vector
%                s=s(1:end-1);
%                s=sprintf('@(t,g_X1%s)reshape(vertcat(%s),[%d,%d])',parslist,s,g_grind.statevars.dim,size(g_grind.syms.(type),2));
%             else
%                s=sprintf('@(t,g_X1%s)[%s]',parslist,s);
%             end
%             h=evalin('base',s);
%         end        
    case 'ExtendedSyst'
        odehandle=i_getodehandle('normal');
        jachandle=i_getodehandle('Jacobian');
        h=@(t,g_X1)runextended(t,g_X1,odehandle,jachandle,g_grind.statevars.dim);
    case {13,'Hessian'} %Hessian (t, x)
        if nargin<3
            opt='';
        end
        if ~isempty(g_grind.syms.Jacobian)&&~isempty(g_grind.syms.Hessian)&&~strcmp(opt,'numonly')
            if ~g_grind.solver.isimplicit
                if isstruct(g_grind.syms.Hessian)
                    h=getsymhandle('Hessian',opt,parslist);
                else
                    s = sprintf('%s,', g_grind.syms.Hessian{:});
                    h=evalin('base', sprintf('@(t,g_X1%s)i_reshapevcat([%d,%d,%d,max(size(t,2),size(g_X1,2))],%s);',parslist,g_grind.statevars.dim,g_grind.statevars.dim,g_grind.statevars.dim,s(1:end-1)));
                end
            else
                error('grind:hess','Hessian not implemented for implicit models');
                %  h = @i_implicit_Jac;
            end

       else
            g_odefunjac = i_getodehandle('Jacobian', parslist,opt);  %#ok<NASGU>
            if strcmp(g_grind.solver.opt.Vectorized, 'on')
                h = eval(sprintf('@(t,g_X1%s)num_DFDXDXv(g_odefunjac, t, g_X1%s);',parslist,parslist));
            else
                h = eval(sprintf('@(t,g_X1%s)num_DFDXDX(g_odefunjac, t, g_X1%s);',parslist,parslist));
            end

        end
    case {'Hessianp'} %der3 (t, x)
        if nargin<3
            opt='';
        end
        h=getsymhandle('Hessianp',opt,parslist);
    case {'der3'} %der3 (t, x)
        if nargin<3
            opt='';
        end
        h=getsymhandle('der3',opt,parslist);
   case {'der4'} %der4 (t, x)
        if nargin<3
            opt='';
        end
        h=getsymhandle('der4',opt,parslist);
    case {'der5'} %der5 (t, x)
        if nargin<3
            opt='';
        end
        h=getsymhandle('der5',opt,parslist);
    case 'polar2cart'
        odehandle=g_grind.solver.polar.odehandle;
        iradius=g_grind.solver.polar.radius;
        itheta=g_grind.solver.polar.theta;
        h=@(t,g_X,varargin)polar2cart(odehandle,iradius,itheta,t,g_X,varargin{:});
    case {15,'quasipot'} %Hessian (t, x)
        g_odefun=evalin('base',sprintf('@(t,g_X1)%s(t,g_X1(1:%d,:)%s)',g_grind.odefile,g_grind.statevars.dim,allpars)); %#ok<NASGU>
      %  g_odefun = i_getodehandle('normal', parslist);
        h = eval(sprintf('@(t,g_X1%s)quasipotfun(g_odefun, t, g_X1%s);',parslist,parslist));
    case {'delays'} %delay function for ddesd
         s=sprintf('%s;',g_grind.dde.lags{:});
         h=evalin('base',sprintf('@(t,g_X1)t-[%s]',s(1:end-1))); 
    case '?' %get help text
        h =sprintf(['0=normal = normal run (incl backw) i_ru\n'...
            '1=singlestep = singlestep  +1\n'...
            '-1=singleback = singlestep   -1\n'...
            '2=potential = potential\n'...
            '3=makeimplicit = make implicit\n'...
            '5=findeq = findeq with fminsearch\n'...
            '6=findeq_fsolve = findeq with fsolve\n'...
            '7=makemap = maps of makemap\n'...
            '8=timesens = test symbolic timesens one parameter needed\n'...
            '9=jac = Jacobian (t, x)\n'...
            'vectorcheck = function to check vectorised'...
            '10=coco = coco ode file {nonautonomic not yet supported)\n'...
            '11=coco_jac = Jacobian (x, p)(coco) [parslist must be cell]\n'...
            'coco_lyap = function for lyapunov exponent \n'...
            '12=coco_jacp = Jacobianp (x, p) (coco) [parslist must be cell]\n'...
            '13=hess = Hessian (t, x) not vector ;pars=vector of parslist\n'...
            '15=quasipot = function for quasipot\n'...
            'delays = delays for ddesd\n']);
        
%             case {0,'normal'} %normal i_ru
%     case {1,'singlestep'} %singlestep  + 1
%     case   {-1,'singleback'} %singlestep   - 1
%     case {2,'potential'} %potential
%     case {3,'makeimplicit'} %make implicit
%     case {5,'findeq'} %findeq with fminsearch
%      case {6,'findeq_fsolve'} %findeq with fsolve
%     case {7,'makemap'} %maps of makema
%     case {8,'timesens'} %test symbolic timesens one parameter needed
%     case {9,'Jacobian'} %Jacobian (t, x)
%     case 'vectorcheck' %parameter vector with time (for i_parcheck)
%     case {10,'coco'} %coco parameter vector
%     case {11, 'coco_jac'} %Jacobian (x, p)
%     case 'coco_lyap' %coco first lyapunov coefficent for Sub/Sup Hopf
%     case {12,'coco_jacp'} %Jacobianp (x, p)
%     case {13,'Hessian'} %Hessian (t, x)
%     case {15,'quasipot'} %Hessian (t, x)
%     case {16,'delays'} %delay function for ddesd
%         

        
        
    otherwise
        error('grind:odehandle','Unknown flag in i_getodehandle');
end
end

function ydot=polar2cart(odehandle,iradius,itheta,t,ycart,varargin)
%transforms a model with polar coordinates to cartesian
%is vectorizable if the odehandle is vectorizable
%(used in null for instance and is very fast)
%At least 2 dimensions are requiered and coordinates of radius and theta
%should be supplied + handle to odefile
%Using g_grind is slow but complete binding to anonymus function is more complex
%to implement
x=ycart(iradius,:);
y=ycart(itheta,:);
theta = atan2(y,x);
radius =sqrt(x.^2+y.^2);
ycart(iradius,:)=radius;
ycart(itheta,:)=theta;
ypoldot=odehandle(t,ycart,varargin{:});
ydot=ypoldot;
ydot(iradius,:)=x./radius.*ypoldot(iradius,:)-y.*ypoldot(itheta,:);
ydot(itheta,:)=y./radius.*ypoldot(iradius,:)+x.*ypoldot(itheta,:);
end

function res=runextended(t,x,odehandle,jachandle,n)
x1=x(1:n);
res=zeros(size(x));
res(1:n)=odehandle(t,x1);
jac=jachandle(t,x1);
Y=reshape(x(n+1:end),n,n);
res(n+1:end)=jac*Y;
end


function jac=checknan(x,p,handle1,numhandle)
%try to replace nan values in the analytic Jacp with numeric
jac=handle1(x,p);
if any(isnan(jac(:)))
   jacn=numhandle(x,p);
   jac(isnan(jac))=jacn(isnan(jac));
%   disp('Nans replaced')
end
end
function  h=getsymhandle(symfield,opt,parslist)
global g_grind;
if ~isempty(g_grind.syms.Jacobian)&&~isempty(g_grind.syms.(symfield))&&~strcmp(opt,'numonly')
    if ~g_grind.solver.isimplicit
        if isfield(g_grind.syms.(symfield),'unique')
           c={g_grind.syms.(symfield).unique(:).equation};
           s = sprintf('%s,', c{:});
        else
            s='[] ';
        end
        h=evalin('base', sprintf('@(t,g_X1%s)i_runsymstruc(g_grind.syms.%s,max(size(t,2),size(g_X1,2)),%s);',parslist,symfield,s(1:end-1)));
    else
        error(sprintf('grind:%s',symfield),'%s not implemented for implicit models',symfield);
        %  h = @i_implicit_Jac;
    end
    
else
    h=[];
end
end

function A=vindex(A,ndx)
A=A(:,ndx);
end

function [data, y] = fun2cocofun(~, data, u, g_funhandle)
%make function g_funhandle(x,p) suitable to use in coco
x = u(data.x_idx);
p = u(data.p_idx);
y=g_funhandle(x,p);
end

function pars=pars2coco(parslist,matcont)
%quite complex for vector models
global g_grind;
if nargin==1
    matcont=false;
end
if isempty(parslist)
    parslist = g_grind.pars;
end

pars = g_grind.pars;
hasbracks=strcontains(parslist,'('); 
if any(hasbracks)
    %check if the parslist contains all elements and removes whole
    %matrices/vectors from the list
    if length(parslist)~=length(unique(parslist))
        error('grind:pars2coco','the parameterlist contains double elements');
    end

    
    %split parameter names and arguments (x)
    parslist2=regexp(parslist,'[a-zA-Z_][a-zA-Z0-9_]*','match','once');
    args=regexp(parslist,'[\(].*','match','once');

    parslist1=unique(parslist2);
    for i=1:length(parslist1)
        siz=evalin('base',sprintf('size(%s)',parslist1{i}));
        ndx1=strcmp(parslist2,parslist1{i});
        fndx=find(ndx1);
        nelem=length(fndx);
        ind=zeros(nelem,1);
        for j=1:nelem
              arg=sscanf(args{fndx(j)},'(%d,%d)');%str2num(args{fndx(j)}); 
              if numel(arg)==1
                  ind(j)=arg;
              elseif ~isempty(arg)
                  ind(j)=sub2ind(siz,arg(1),arg(2));
              end

        end

        diffind=diff(ind);
        if nelem>1&&nelem==prod(siz)&&~isempty(diffind)&&all(diffind==1)
            fndx1=find(ndx1,1);
            parslist(fndx1)=parslist1(i);
            args{fndx1}='';
            ndx1=~ndx1;
            ndx1(fndx1)=true;
            parslist=parslist(ndx1);
            args=args(ndx1);
            parslist2=regexp(parslist,'[a-zA-Z_][a-zA-Z0-9_]*','match','once');
        end

    end
    args=regexp(parslist,'[\(].*','match','once');
%    args=strrep(strrep(args,'(','['),')',']');
    hasbracks=strcontains(parslist,'('); %   
end

if any(hasbracks)&&g_grind.statevars.vector
    numels=ones(size(parslist));
    for i=1:length(parslist)
        if ~hasbracks(i)
           numels(i)=evalin('base',sprintf('numel(%s)',parslist{i}));
        end

    end

    ks=cumsum(numels);
 %   parslist1=unique(parslist2);
    allpars=pars;
    %g_par=cell(size(parslist1));
    for i=1:length(allpars)
        siz=evalin('base',sprintf('size(%s)',allpars{i}));
        ndx1=strcmp(allpars{i},parslist2);
        ndx2=ndx1&hasbracks;
        if any(ndx2)
           fndx=find(ndx2);
           ind=zeros(size(fndx));
           for j=1:length(fndx)
               arg=sscanf(args{fndx(j)},'(%d,%d)');
               if numel(arg)==1
                   ind(j)=arg;
               else
                  ind(j)=sub2ind(siz,arg(1),arg(2));
               end

           end
           repind=sprintf('%d;',ind);
           k1=min(ks(ndx1));
           k2=max(ks(ndx1));
           if k2+1-k1~=size(ind)
               error('grind:getodehandle','Error generating function handle for coco, elements of vector should be next to eachother');
           end
           if matcont
                pars{i}=sprintf('rep_elem(%s,{varargin{%d:%d}],[%s])',allpars{i},k1,k2,repind(1:end-1));
           else
               pars{i}=sprintf('rep_elem(%s,g_par(%d:%d),[%s])',allpars{i},k1,k2,repind(1:end-1));
           end
        elseif ~any(ndx1)
            pars{i}=allpars{i};
        else
            fndx=find(ndx1,1);
            k2=ks(fndx);
            if fndx==1
                k1=1;
            else
                k1=ks(fndx-1)+1;
            end
            if prod(siz)>1
                if matcont
                    pars{i}=sprintf('reshape([varargin{%d:%d}],[%d,%d])',k1,k2,siz(1),siz(2));
                else
                    pars{i}=sprintf('reshape(g_par(%d:%d),[%d,%d])',k1,k2,siz(1),siz(2));
                end
            else
                if matcont
                   pars{i}=sprintf('varargin{%d}',k1);
                else
                   pars{i}=sprintf('g_par(%d)',k1);
                end
            end
        end
    end
else
    %No brackets, simpler
    [ndx1, ndx2] = ismember(pars, parslist);
    ndx = find(ndx1); %some tricks to preserve order of pars
    if length(ndx) ~= length(parslist)
        unknownpars = setdiff(parslist, g_grind.pars);
        error('grind:i_getodehandle','Unknown or double parameter(s): %s',sprintf('%s ',unknownpars{:}))
    end

    ndx2 = ndx2(ndx1);
    ndx = ndx(ndx2);
    if ~g_grind.statevars.vector
        if matcont
            g_par = str2cell(sprintf('[varargin{%d,:}]\n',1:length(parslist)));
        else
           g_par = str2cell(sprintf('g_par(%d,:)\n',1:length(parslist)));
        end
        pars(ndx) = g_par(:);
    else
        g_par=cell(size(parslist));
        k=1;
        for i=1:length(parslist)
            siz=evalin('base',sprintf('size(%s)',parslist{i}));
            k2=k+prod(siz);
            if prod(siz)>1
                if matcont
                   g_par{i}=sprintf('reshape([varargin{%d:%d}],[%d,%d])',k,k2-1,siz(1),siz(2));
                else
                   g_par{i}=sprintf('reshape(g_par(%d:%d),[%d,%d])',k,k2-1,siz(1),siz(2));
                end
            else
                if matcont
                    g_par{i}=sprintf('varargin{%d}',k);
                else
                    g_par{i}=sprintf('g_par(%d)',k);
                end
            end
            k=k2;
        end
        pars(ndx) = g_par(:);
    end
end
end
function [parslist,allpars]=checkpars(parslist,allpars)
%replace A(2,3),A(4,4) in the parslist with A_2_3_ ; A_4_4_
%replace A in allpars with rep_elem(A,[A_2_3_;A_4_4_],[ndx1ndx2])
if ischar(parslist)&& any(strcontains(parslist,'('))
    vecpars= regexp(parslist,'[a-zA-Z_][a-zA-Z0-9_]*[\(][0-9]*([,][0-9])*[\)]','match');
    parslist=regexprep(vecpars,'[\(,\)]','_');
    vecpars1=regexp(vecpars,'[a-zA-Z_][a-zA-Z0-9_]*','match','once');
    args=regexp(vecpars,'[\(].*','match','once');
    %args=strrep(strrep(args,'(','['),')',']');
    allpars1=regexp(allpars,',','split');
    allpars1(1)=[];
    allpars2=allpars1;
    for i=1:length(allpars1)
        ndx=find(strcmp(vecpars1,allpars1{i}));
        if ~isempty(ndx)
           reppars=sprintf('%s;',parslist{ndx});
           siz=evalin('base',sprintf('size(%s);',allpars1{i}));
           ind=zeros(size(ndx));
           for j=1:length(ndx)
               if sum(siz>1)==1
                  ind(j)=sscanf(args{ndx(j)},'(%d)');
               %   ind(j)=sub2ind(siz,arg(1),1);
               else
                  arg=sscanf(args{ndx(j)},'(%d,%d)');
                 % arg=str2num(args{ndx(j)}); 
                  ind(j)=sub2ind(siz,arg(1),arg(2));
               end
           end
           repind=sprintf('%d;',ind);
           allpars1{i}=sprintf('rep_elem(%s,[%s],[%s])',allpars1{i},reppars(1:end-1),repind(1:end-1));
        else
           allpars1{i}=allpars2{i};
        end
    end
    allpars=sprintf(',%s',allpars1{:});
    parslist=sprintf(',%s',parslist{:});
end
end
function g_X2=quasipotfun(g_odefun, t, g_X1, varargin)  %#ok<DEFNU>
%g_h=g_odefun(t,g_X1(1:end-1,:),varargin{:});
g_h=g_odefun(t,g_X1,varargin{:});
g_X2=[g_h;-sum(g_h.^2,1)];
end
function y=newton4differ(t,y,odef,Jfun,TOL,MAXITS,iters)
% [Y,isConverged]=newton4euler(f,x,ytranspose,Y,h)
% special function to evaluate Newton's method for back_euler
yold=y;
for n=1:MAXITS
  ynew=y;
  for i=1:iters
      ynew=feval( odef, t, ynew);
  end

  fValue = ynew-yold;
  if isempty(Jfun)||iters>1
     fPartial = i_calcjac(1, iters, y);%-eye(length(Y));
  else
     fPartial = Jfun(t,y);%-eye(length(Y));
  end

  increment=fPartial\fValue;
  y = y - increment;
  if norm(increment,inf) < TOL*norm(y,inf)
    return
  end
end
warning('grind:backdiffer','Backwards difference equation: Newton method failed to converge');
end

function yy = implicitstep(t, N0, varargin)
global g_grind;
if nargin > 2
    h = @(t, g_X1,g_X2)g_grind.solver.h(t, g_X1,g_X2, varargin{:});
else
    h = g_grind.solver.h;
end

try
    [~, yy] = decic(h, t, N0, ones(size(N0)), rand(size(N0)), zeros(size(N0)));
catch
    yy = zeros(size(N0)) + NaN;
end
end
function N1 = iterstep(t, N0, handl, iters, varargin)
N1=N0;
for i = 1:iters
    N1 = handl(t, N1,varargin{:});
end
end
%N1=N1-N0;
function Nres = findeqmask(x)
global g_grind;
N0 = g_grind.findeq.n0;
N0(g_grind.findeq.mask(:)) = x(g_grind.findeq.mask(:));
Nres = feval(g_grind.solver.h, 1,N0);
Nres = Nres(g_grind.findeq.mask(:));
% function parslist=makecocoparlist(parslist,pars)
% g_par=str2cell(sprintf('g_par(%d)\n',1:length(parslist)));
% ndx=ismember(pars,parslist);
% pars(ndx)=g_par(:);
% parslist=sprintf(',%s',pars{:});
%
end
function g_X2=makeimplicit(t,g_X1,g_X3,g_odefun) %#ok<DEFNU>
%some problems if vectorized
g_X2=g_odefun(t,g_X1);
size3=size(g_X3,2);
size2=size(g_X2,2);
if size3>size2
    g_X2=g_X3-repmat(g_X2,1,1+size(g_X3,2)-size(g_X2,2));
elseif size2>size3
    g_X2=repmat(g_X3,1,1+size(g_X2,2)-size(g_X3,2))-g_X2;
else
    g_X2=g_X3-g_X2;
end
end
function J = num_DFDX(F,t, x,varargin)  %#ok<DEFNU>
%Numerical Jacobian non-vectorized
[m, n] = size(x);
fr = F(t,x(:, 1),varargin{:});
l  = size(fr, 1);
J  = zeros(l, m, n);
for j = 1:n
    x0 = x(:, j);
    h  = 1.0e-8 * ( 1.0 + abs(x0) );
    hi = 0.5 ./ h;
    for i = 1:m
        xx       = x0;
        xx(i)    = x0(i) + h(i);
        fr       = F(t,xx,varargin{:});
        xx(i)    = x0(i) - h(i);
        fl       = F(t,xx,varargin{:});
        J(:, i, j) = hi(i) * (fr - fl);
    end
end
end

function J = num_DFDXv(F,t, x,varargin)  %#ok<DEFNU>
%Numerical Jacobian vectorized
x=x(:,:);
%p=p(:,:);
[m, n] = size(x);
idx   = repmat(1:n, [m 1]);
x0    = x(:, idx);
%p0    = p(:, idx);
idx = repmat(1:m, [1 n]);
idx = sub2ind([m m * n], idx, 1:m * n);
h = 1.0e-8 * (1.0 + abs(x0(idx)));
x = x0;
x(idx) = x0(idx) + h;
fr     = F(t,x,varargin{:});
x(idx) = x0(idx) - h;
fl     = F(t,x,varargin{:});
l  = size(fr, 1);
hi = repmat(0.5 ./ h, [l 1]);
J  = reshape(hi .* (fr - fl), [l m n]);
end
function J = num_DFDXDXv(F,t, x,varargin)  %#ok<DEFNU>
%numerical Hessian vectorized
x=x(:,:);
%p=p(:,:);
[m, n] = size(x);
idx   = repmat(1:n, [m 1]);
x0    = x(:, idx);
%p0    = p(:, idx);
idx = repmat(1:m, [1 n]);
idx = sub2ind([m m * n], idx, 1:m * n);
h = 1.0e-4 * (1.0 + abs(x0(idx)));
x = x0;
x(idx) = x0(idx) + h;
fr     = F(t,x,varargin{:});
x(idx) = x0(idx) - h;
fl     = F(t,x,varargin{:});
l  = size(fr, 1);
hi = repmat(0.5 ./ h, [l*l, 1]);
hi= reshape(hi,[l l m*n]);
J  = reshape(hi .* (fr - fl), [l l m n]);
end
function [J,hs] = num_DFDXDX(F,t, x,varargin)  %#ok<DEFNU>
%numerical Hessian non-vectorized
[m, n] = size(x);
fr = F(t,x(:, 1),varargin{:});
l  = size(fr, 1);
J  = zeros(l,l, m, n);
hs=J;
for j = 1:n
    x0 = x(:, j);
    h  = 1.0e-4 * ( 1.0 + abs(x0) );
    hi = 0.5 ./ h;
    for i = 1:m
        xx       = x0;
        xx(i)    = x0(i) + h(i);
        fr       = F(t,xx,varargin{:});
        xx(i)    = x0(i) - h(i);
        fl       = F(t,xx,varargin{:});
        J(:,:, i, j) = hi(i) * (fr - fl);
        hs(:,:,i,j)=hi(i);
    end
end
end
