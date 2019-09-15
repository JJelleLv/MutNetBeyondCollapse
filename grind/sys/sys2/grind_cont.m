classdef grind_cont < handle
    %abstract class
    properties
        settings = struct('grind',[],'derived',[]);
        curves = struct('freepars',[],'ctype',{},'frompoint','','data',[],'color',[],'propndx',[],'results',struct('startndx',1,'s',[],'stabil',[]));
        points = struct('id', {},'p0',[],'x0',[],'data',[],'propndx',[]);
    end
    methods (Abstract)
        args=set(obj,varargin);
        [args,found]=get(obj,varargin);
        %cont(obj,from,ctype,varargin);
        %    add_points(obj,newpoints,p0);
    end
    methods (Abstract,Access = protected)
        [point,pntndx]=create_point(obj,spoint,x0,p0)
        init_engine(obj,close);
        expand_curve(obj,frompoint,ctype);
        run_point(obj,frompoint,ctype);
    end
    
    methods (Static, Access = public)
        %% load a previously saved sessions
        function obj=load_session(filename)
            global g_grind;
            if nargin<1||isempty(filename)
                [filename, pathname] = uigetfile({'*.mat';'*.*'}, 'Load');
                if ~isempty(pathname)
                    cd(pathname)
                end
            end
            [~,nam,ext]=fileparts(g_grind.inifile);
            curr_inifile=[nam ext];
            s=load(filename,'cont_session');
            if ~isfield(s,'cont_session')
                error('grind:conteq','File does not contain a conteq session')
            end
            [~,nam,ext]=fileparts(s.cont_session.settings.derived.inifile);
            new_inifile=[nam ext];
            if ~strcmp(curr_inifile,new_inifile)
                btn=questdlg(sprintf('This file uses another inifile (%s), Ok to open this inifile?',new_inifile),'Load conteq session','Yes','No','Yes');
                if strcmp(btn,'Yes')
                    use(s.cont_session.settings.derived.inifile);
                    obj=grind_cont.load_session(filename);
                    return;
                else
                    error('grind:conteq','File incompatible with current session: other ini-file (%s)',...
                        s.cont_session.settings.derived.inifile);
                end
            end
            switch s.cont_session.settings.derived.engine
                case 'matcont'
                    obj=grind_matcont;
                case 'coco'
                    obj=grind_coco;
                otherwise
                    error('grind:conteq','Unknown continuation engine')
            end
            obj.settings=s.cont_session.settings;
            obj.curves=s.cont_session.curves;
            obj.points=s.cont_session.points;
        end
    end
    methods (Access = public)
        
        function obj = grind_cont
            global g_grind;
            %default settings
            obj.settings.derived.inifile=fullfile(cd,g_grind.inifile);
            obj.settings.derived.opened=false;
            obj.settings.derived.ndim=g_grind.statevars.dim;
            obj.settings.derived.frompoint='';
            obj.settings.derived.ctype='';
            obj.settings.derived.ischanged=false;
            obj.settings.derived.statevars=transpose(i_statevars_names);
            obj.settings.grind.symbolic=~g_grind.statevars.vector;
            obj.set('-defaults');
        end
        
        %% open_matcont add search path
        function saveas(obj,filename)
            
            oldpath=cd;
            if nargin<2||isempty(filename)
                [filename, pathname] = uiputfile({'*.mat';'*.m';'*.*'}, 'Save session as');
                if ~ischar(filename)
                    return;
                end
                if ~isempty(pathname)
                    cd(pathname)
                end
            end
            cd(oldpath)
            if strcontains(filename,'.mat')
                cont_session.settings=obj.settings;
                cont_session.curves=obj.curves;
                cont_session.points=obj.points; 
                save(filename,'cont_session');
                if ~strcontains(filename,'.')
                    fprintf('%s session saved to "%s.mat"\n',obj.settings.derived.engine,filename);
                else
                    fprintf('%s session saved to "%s"\n',obj.settings.derived.engine,filename);
                end
            elseif strcontains(filename,'.m')
                fid=fopen(filename,'w');
                fprintf(fid,'%%This script is created by GRIND\n%%Date: %s\n\n%%open the model\nuse(''%s'');\n\n',datestr(now()),obj.settings.derived.inifile);
                fprintf(fid,'%%create object "g_cont" that uses the %s engine\nglobal g_cont;\ng_cont=%s;\n\n%%first add all initial points\n',obj.settings.derived.engine,class(obj));
                ndx=find(obj.getndx('points','label','P')|obj.getndx('points','label','EP'));
                if ~isempty(ndx)
                    g_set{1}=obj.points(ndx(1)).p0;
                    parnr=ones(size(ndx));
                    for i=2:length(ndx)
                        j=1;
                        while j<=length(g_set)
                            if struccmp(obj.points(ndx(i)).p0,g_set{j})
                                parnr(i)=j;
                                break;
                            end
                            j=j+1;
                        end
                        if j>length(g_set)
                            parnr(i)=j;
                            g_set{j}=obj.points(ndx(i)).p0;
                        end
                    end
                end
                for j=1:length(g_set)
                    fprintf(fid,'g_set{%d}=%s;\n',j,mat2str(g_set{j},10));
                end
                for i=1:length(ndx)
                    pnt=obj.points(ndx(i));
                    fprintf(fid,'g_cont.add_points(''id'',''%s'',''msg'',''%s'',''p0'',g_set{%d},...\n  ''x0'',%s);\n',pnt.id,pnt.data.msg,parnr(i),mat2str(pnt.x0,10));
                end
                fprintf(fid,'\n');
                ndx=find(obj.getndx('points','label','P'));
                fprintf(fid,'g_cont.select_point(''%s'',true); %%select the last parameter values in GRIND\n\n',obj.points(ndx(end)).id);
                for i=1:length(obj.curves)
                    if i==1
                        sett=obj.curves(i).data.settings;
                        opts=fieldnames(obj.curves(i).data.settings);
                        changed=true;
                    else
                        [changed,opts]=struccmp(sett,obj.curves(i).data.settings);
                        changed=~changed;
                        sett=obj.curves(i).data.settings;
                        opts=opts';
                    end
                    if changed
                        vals=opts;
                        for j=1:length(opts)
                            if ~isfield(sett,opts{j})
                                vals{j}=obj.get('-properties',opts{j}).default;
                            else
                                vals{j}=sett.(opts{j});
                            end
                            vals{j}=any2str(vals{j});
                        end
                        setts=[opts';vals'];
                        setts=strtrim(sprintf('''%s'',%s,...\n    ',setts{:}));
                        fprintf(fid,'%%change the %s/GRIND settings\ng_cont.set(%s);\n\n',obj.settings.derived.engine,setts(1:end-4));
                    end
                    %id=regexp(obj.curves(i).frompoint,'[A-Za-z+]*','match','once');
                    ndx=obj.getndx('points','id',obj.curves(i).frompoint);
                    frompnt=obj.points(find(ndx,1));
                    fprintf(fid,'%%Continue %s to a %s curve\n',obj.pointprops(frompnt.propndx).descr,obj.curveprops(obj.curves(i).propndx).descr);
                    fprintf(fid,'g_cont.cont(''%s'',''%s'');\n\n', obj.curves(i).frompoint,obj.curves(i).ctype);
                    if isfield(obj.curves(i).data,'expanded')
                        for k=1:obj.curves(i).data.expanded
                            fprintf(fid,'%%Expand the previous curve\ng_cont.cont(''%s'',''%s'',''expand'');\n\n', obj.curves(i).frompoint,obj.curves(i).ctype);
                        end
                    end
                end
                fprintf(fid,'g_cont.show;\n\ng_cont.plot;\n\ndisp(''You can use <a href="matlab:conteq">conteq</a> to continue working with this session'')\n\n');
                fclose(fid);
            end
        end
        
        %% open continuation engine add search path
        function open(obj)
            init_engine(obj);
            obj.settings.derived.opened=true;
        end
        
        %% close_matcont remove search path
        function close(obj)
            init_engine(obj,true)
            obj.settings.derived.opened=false;
        end
        
        
        %%
        function cont(obj,from,ctype,varargin)
            %use this preferrably varargin are additional arguments for
            %initialization as property value pairs.
            %If information is missing, the user is prompted
            %varargin1 = 'expand' for expanding the (existing) curve
            if nargin==1&&~isempty(obj.settings.derived.frompoint)
                from=obj.settings.derived.frompoint;
                if ~isempty(obj.settings.derived.ctype)
                    ctype=obj.settings.derived.ctype;
                end
            end
            
            if iscell(from)
                if nargin<3
                    ctype=[];
                end
                for i=1:length(from)
                    cont(obj,from{i},ctype,varargin{:});
                end
                return;
            end
            obj.open;
            if isempty(obj.settings.derived.freepars)
                error('grind:cont:nopars','No free parameters selected');
            end
            frompoint=obj.points(getndx(obj,'points','id',from));
            if isempty(frompoint)
                error('grind:matcont','Special point "%s" not found',from);
            end
            if (nargin<3||isempty(ctype))&&isempty(obj.settings.derived.ctype)
                %the first is the default curve type
                ctypes=obj.pointprops(frompoint.propndx).ctypes;
                ctype=ctypes{1};
            end
            i_waitbar(0, 2, sprintf('running %s',obj.settings.derived.engine),sprintf('Continuing "%s"->"%s" \n(cannot show progress)',frompoint.id,ctype))
            obj.settings.derived.frompoint=frompoint.id;
            obj.settings.derived.ctype=ctype;
            doexpand= ~isempty(varargin)&& strcmp(varargin{1},'expand');
            if length(varargin)>1&&doexpand
                obj.set(varargin{2:end});
            elseif ~isempty(varargin)&&~doexpand
                obj.set(varargin{:});
            end
            if doexpand
                expand_curve(obj,frompoint,ctype);
            else
                run_point(obj,frompoint,ctype)
            end
            i_waitbar([]);
            obj.settings.derived.ischanged=true;
            obj.close;
        end
        %%
        function changed=update_points(obj)
            changed=false;
            ndxp=find(obj.getndx('points','ptype','P'),1,'last');
            point=obj.points(ndxp);
            if isempty(point)
                x0_changed=true;
                p0_changed=true;
                id='P1';
            else
                id=point.id;
                p0_changed=sum((point.p0-par('-vector')).^2,1);
                x0_changed=sum((point.x0-i_initvar).^2,1);
            end
            hasEP=any(obj.getndx('points','id','EP1'));
            if p0_changed||x0_changed
                changed=true;
                if isempty(ndxp)
                    N0=i_initvar;
                    s2 = sprintf('%g ', round(N0 * 10000) / 10000);
                    obj.add_points(struct('id',id,'label','P','msg',sprintf('initial point: %s',s2),'data',[],'N0',N0));
                    ndxp=obj.getndx('points','id',id);
                end
                if ~obj.points_in_use(id)
                    obj.points(ndxp).x0=i_initvar;
                    obj.points(ndxp).data.msg=sprintf('initial point: %s',sprintf('%g ', obj.points(ndxp).x0));
                    obj.points(ndxp).p0=par('-vector');
                else
                    N0=i_initvar;
                    s2 = sprintf('%g ', round(N0 * 10000) / 10000);
                    obj.add_points(struct('id',id,'label','P','msg',sprintf('initial point: %s',s2),'data',[],'N0',N0));
                end
            end
            if p0_changed&&hasEP
                changed=true;
                epndx=obj.getndx('points','ptype','EP');
                inuse=obj.points_in_use({obj.points.id});
                obj.points=obj.points(epndx&inuse|~epndx);
                %     EPs={obj.points(obj.points_in_use(id)).id};
                obj.add_points(findeqs('maxtime',5));
            end
            if p0_changed
                epndx=obj.getndx('points','ptype','EP')|obj.getndx('points','ptype','P');
                p0_values=horzcat(obj.points(epndx).p0);
                changedpars=sum(p0_values-repmat(p0_values(:,1),1,size(p0_values,2)),2)>0;
                for i=find(epndx)
                    s=sprintf('%s=%%g ',obj.settings.derived.allpars{changedpars});
                    s1=sprintf(s,obj.points(i).p0(changedpars));
                    if any(changedpars)
                        s1=sprintf(': %s,',strtrim(s1));
                    else
                        s1=': ';
                    end
                    obj.points(i).data.msg= regexprep(obj.points(i).data.msg,'(:[^,]*)(=?,)|:', s1);
                end
            end
        end
        
        
        %% Show the points and the curves
        function [points,curves,par_sets]=show(obj,dims,filter)
            function [str]=mymat2str(aval)
                if numel(aval)<10
                    str= ['[' strtrim(sprintf('%.4g ',aval)) ']'];
                else
                    str= ['[' strtrim(sprintf('%.4g ',aval(1:9))) sprintf(' ... %.4g]',aval(end))];
                end
            end
            if nargin<3
                filter=[];
            end
            if nargin<2
                dims=[];
            end
            %get the name of the calling object
            varname=inputname(1); %this only works if show is called from the workspace
            vars=evalin('base','whos');
            if ~any(strcmp({vars.name},varname))
                f=find(strcmp({vars.class},class(obj)));
                if ~isempty(f)
                    if length(f)==2
                        f1=f;
                        f=f(1);
                        for k=1:length(f1)
                            if strcmp(vars(f1(k)).name,'g_cont')
                                %if there are more instances of the class
                                %give priority to g_cont as name
                                f=f1(k);
                                break;
                            end
                        end
                    end
                    varname=vars(f).name;
                else
                    varname=[];
                end
            end
            
            if ~isempty(dims)
                fpoint=find(obj.getndx('points','npars',dims));
            else
                fpoint=find(true(size(obj.points)));
            end
            if ~isempty(filter)
                names=fields(filter);
               [~,filterpars]=ismember(names,obj.settings.derived.allpars);
               ndx=true(size(fpoint));
               ppars=[obj.points(fpoint).p0];
               for i=1:length(filterpars)
                    ndx=ndx & ppars(filterpars(i),:)==filter.(names{i});
               end
               fpoint=fpoint(ndx);
            end
                
            res1=cell(length(fpoint),1);
            
            curveid=zeros(size(obj.points));
            for i=1:length(curveid)
                f=find(obj.getndx('curves','id',obj.points(i).id),1);
                if ~isempty(f)
                    curveid(i)=f;
                end
            end
            %get changed parameters others than freepars
            ppars=[obj.points(:).p0];
            isfreepar=false(size(ppars));
            for i=1:length(curveid)
                if curveid(i)~=0
                    freepars=obj.curves(curveid(i)).freepars;
                    isfreepar(freepars,i)=true;
                end
            end
            changedpars=[];
            for i=1:size(ppars,1)
                upar=unique(ppars(i,~isfreepar(i,:)));
                if length(upar)>1
                    changedpars(end+1)=i; %#ok<AGROW>
                end
            end
            if nargout>2
                %calculate the unique sets of parameters used
                if isempty(changedpars)
                    par_sets={};
                else
                    ppars1=ppars;
                    ppars1(isfreepar)=Inf;
                    ppars1=unique(transpose(ppars1(changedpars,:)),'rows');
                    ppars1=ppars1(~any(isinf(ppars1),2),:);
                    parnames=obj.settings.derived.allpars(changedpars);
                    par_sets=struct(parnames{1},num2cell(ppars1(:,1)));
                    for i=2:length(parnames)
                        for j=length(par_sets):-1:1
                            par_sets(j).(parnames{i})=ppars1(j,i);
                        end
                    end
                end
            end
            for i=1:length(fpoint)
                if ~isempty(varname)&&(nargout==0)
                    %the name of the calling object is needed here to
                    %let the user select a special point in a hyperlink
                    id=sprintf('<a href="matlab:%s.select_point(''%s'')">%s</a>',varname,obj.points(fpoint(i)).id,obj.points(fpoint(i)).id);
                else
                    %if there is not varname we cannot make a link
                    id=obj.points(fpoint(i)).id;
                end
                
                %                 if any(strcmp(obj.pointprops(obj.points(fpoint(i)).propndx(1)).ptype,{'P','EP'}))
                %                     res1{i}=sprintf('%s - %s',id,obj.points(fpoint(i)).data.msg);
                %                 else
                vals=obj.points(fpoint(i)).x0;
                if curveid(fpoint(i))==0
                    freepars=changedpars;
                else
                    freepars=union(changedpars,obj.curves(curveid(fpoint(i))).freepars);
                end
                if isempty(freepars)
                    res1{i}=sprintf('%s - %s, %s',id,obj.points(fpoint(i)).data.msg,mymat2str(vals));
                elseif length(freepars)==1
                    parval=obj.points(fpoint(i)).p0(freepars);
                    res1{i}=sprintf('%s - %s, %s=%.4g, %s',id,obj.points(fpoint(i)).data.msg,obj.settings.derived.allpars{freepars(1)},parval,mymat2str(vals));
                else
                    parvals=sprintf('%.4g,',obj.points(fpoint(i)).p0(freepars));
                    pars=sprintf('%s,',obj.settings.derived.allpars{freepars});
                    res1{i}=sprintf('%s - %s, [%s]=[%s], %s',id,obj.points(fpoint(i)).data.msg,pars(1:end-1),...
                        parvals(1:end-1),mymat2str(vals));
                end
                %                end
            end
            if nargout==0
                disp('Points:')
                fprintf('%s\n',res1{:})
            else
                points=res1;
            end
            if ~isempty(dims)
                fcurve=find(obj.getndx('curves','npars',dims));
            else
                fcurve=find(true(size(obj.curves)));
            end
            res1=cell(length(fcurve),1);
            for i=1:length(fcurve)
                curve=obj.curves(fcurve(i));
                pars=outfun(obj,'curvendx',fcurve(i),'fun',obj.settings.derived.allpars(curve.freepars));
                if ~isempty(varname)&&(nargout==0)
                    %the name of the calling object is needed here to
                    %let the user select a special point in a hyperlink
                    No=sprintf('<a href="matlab:i_curve_cont_dlg(%s,%d)">Curve %d</a>',varname,fcurve(i),fcurve(i));
                else
                    %if there is not varname we cannot make a link
                    No=sprintf('Curve %d',fcurve(i));
                end
                
                if length(curve.freepars)==1
                    res1{i}=sprintf('%s - %s->%s [%dx%d] %s=[%g to %g]',No,curve.frompoint,curve.ctype,size(pars,1),obj.settings.derived.ndim,obj.settings.derived.allpars{curve.freepars},min(pars),max(pars));
                else
                    res1{i}=sprintf('%s - %s->%s [%dx%d] %s=[%g to %g], %s=[%g to %g]',No,curve.frompoint,curve.ctype,size(pars,1),...
                        obj.settings.derived.ndim,obj.settings.derived.allpars{curve.freepars(1)},min(pars(:,1)),max(pars(:,1)),obj.settings.derived.allpars{curve.freepars(2)},min(pars(:,2)),max(pars(:,2)));
                end
            end
            if nargout==0
                if ~isempty(res1)
                    disp(' ')
                    disp('Curves:')
                    fprintf('%s\n',res1{:})
                end
            else
                curves=res1;
            end
        end
        
        %% clear all curves and points or settings
        function clear(obj,target)
            %obj.settings = struct('matcont',[],'grind',[]);
            if nargin<2
                target='';
            end
            if isempty(target)||strcmpi(target,'settings')
                obj.set('-defaults');
            end
            if isempty(target)||strcmpi(target,'curves')
                if isfield(obj.settings.derived,'frompoint')
                    obj.settings.derived.frompoint='';
                    obj.settings.derived.ctype='';
                end
                obj.curves = struct('freepars',[],'ctype',{},'frompoint','','data',[],'color',[],'propndx',[],'results',struct('startndx',1,'s',[],'stabil',[]));
                obj.points = struct('id',{},'p0',[],'x0',[],'data',[],'propndx',[]);
            end
        end
        %% flexible function to get data out of the session
        function [out,xtra]=outfun(obj,varargin)
            args=i_parseargs('fun(s#c),freepars(i>0#c),id(s#c),filter(l),curvendx(i#l),pointndx(i#l)','fun,id,freepars','',varargin);
            if ~isfield(args,'freepars')
                args.freepars=obj.settings.derived.freepars;
            end
            if iscell(args.freepars)
                if length(args.freepars)==2
                    freepars(2) =find(strcmp(args.freepars{2},obj.settings.derived.allpars));
                end
                freepars(1) =find(strcmp(args.freepars{1},obj.settings.derived.allpars));
                args.freepars=freepars;
            end
            if ~isfield(args,'fun')||isempty(args.fun)
                %select default functions
                if numel(args.freepars)==1
                    args.fun=[obj.settings.derived.allpars(args.freepars) i_statevars_names(1:2)];
                else
                    args.fun=obj.settings.derived.allpars(args.freepars);
                end
            end
            out=[];
            xtra=[];
            crvs=[];
            pnts=[];
            if isfield(args,'id')&&~isempty(args.id)
                ndx=obj.getndx('points','id',args.id);
                if any(ndx)
                    pnts=obj.points(ndx);
                else
                    ndx=obj.getndx('curves','ctype',args.id);
                    crvs=obj.curves(ndx);
                end
            else
                if isfield(args,'curvendx')
                    crvs=obj.curves(args.curvendx);
                    if ~isempty(crvs)
                        args.freepars=crvs(1).freepars;
                    end
                end
                if isfield(args,'pointndx')
                    pnts=obj.points(args.pointndx);
                end
            end
            %the id is assumed to be from a curve (only tested for EP)
            s=[];
            if ~isempty(crvs)
                s=[];
                xtra=[];
                for i=1:length(crvs)
                    %all parameters inculding the free pars
                    [s1,xtra1]=obj.curveprops(crvs(i).propndx).resfun(obj,crvs(i));
                    s=catstruc(s,s1,'pars');
                    xtra=catstruc(xtra,xtra1);
                end
            elseif ~isempty(pnts)
                %the id is from one or more points;
                xtra.ids={pnts.id};
                x0s=cell(size(pnts));
                for i=1:length(pnts)
                    x0s{i}=pnts(i).x0.';
                end
                %                 numelx0=cellfun(@(x)numel(x),x0s);
                %                 if all(numelx0==numelx0(1))
                %                    siz2=size(pnts(1).x0,2);
                %                    parvals=permute(reshape(repmat([pnts.p0],[siz2,1]),[numel(pnts(1).p0),numel(x0s),siz2]),[2,1,3]);
                %                    parvals=reshape(permute(parvals,[1 3 2]),[numel(x0s)*siz2,numel(pnts(1).p0)]);
                %                    Y=transpose([pnts(:).x0]);
                %                    s=struct('parvalues',parvals,'pars',{obj.settings.derived.allpars},'Y',permute(Y,[3 2 1]),'t',zeros(size(Y,1),1),'perm',[]);
                %                 else
                Y=transpose(x0s);
                s=struct('parvalues',transpose([pnts.p0]),'pars',{obj.settings.derived.allpars},'Y',{Y},'t',zeros(size(Y,1),1),'perm',[]);
                %              end
            end
            if isfield(args,'filter')&&args.filter
                s=obj.filter_par_var(s,args.freepars);
            end
            if ~isempty(s)
                if iscell(s.Y)
                    snew=s;
                    out=cell(size(s.Y));
                    for i=1:length(s.Y)
                        snew.Y=s.Y{i};
                        snew.t=repmat(s.t(i),[size(snew.Y,1),1]);
                        snew.parvalues=s.parvalues(i,:);
                        out{i}=outfun(args.fun,'datastruc',snew);
                    end
                else
                    out=outfun(args.fun,'datastruc',s);
                    xtra.iscurve=false;
                    if isfield(xtra,'stabil')&&~isempty(xtra.stabil)&&size(xtra.stabil,3)>1
                        xtra.iscurve=true;
                        xtra.stabil=transpose(reshape(permute(xtra.stabil, [2, 1, 3]),...
                            size(xtra.stabil, 2), size(xtra.stabil, 1) * size(xtra.stabil, 3)));
                    end
                end
            end
        end
        %% select data of a point
        function select_point(obj,id,silent)
            if nargin<3
                silent=false;
            end
            ndx=getndx(obj,'points','id',id);
            if any(ndx)
                point=obj.points(ndx);
                i_keep(point.x0);
                ndx=1;
                p0=point.p0;
                allpars=obj.allpars_vector;
                for i=1:length(allpars)
                    siz=evalin('base',sprintf('size(%s);',allpars{i}));
                    assignin('base',allpars{i},reshape(p0(ndx:ndx+prod(siz)-1),siz(1),siz(2)));
                    ndx=ndx+prod(siz);
                end
                descr='';
                if isfield(point.data,'msg')
                    descr=sprintf(' (%s)',point.data.msg);
                end
                if ~silent
                    fprintf('Selected "%s"%s\n',id,descr);
                end
            else
                error('grind:cont:notfoud','Point "%s" not found', id);
            end
            ndx=obj.getndx('pointprops','ptype',regexp(id,'([A-Za-z+0-9]*_)|([A-Za-z+]*)','match','once'));
            if ~silent
                fprintf('%s = %s (codim %d)\n',obj.pointprops(ndx).ptype,obj.pointprops(ndx).descr,obj.pointprops(ndx).codim);
            end
        end
        
        %% flexible function to select  points or curves
        function [ndx,differ]=getndx(obj,point_or_curve,fieldname,value)
            function x0=get_x0(x0)
                if size(x0,2)==1
                    return;
                else
                    x0=max(x0(:,1:end-1),[],2);
                end
            end
            differ=[];
            strct=obj.(point_or_curve);
            if isempty(strct)
                ndx=[];
                return;
            end
            
            switch fieldname
                case {'id','ctype','frompoint'}
                    %you can use a simple wildcard that replaces the number
                    %H* finds H1 H2 .. H10 etc. (but not HH1)
                    if strcmp(point_or_curve,'curves')&&strcmp(fieldname,'id')
                        %find the curve(s) where point id is created
                        ndx=false(size(obj.curves));
                        for i=1:length(obj.curves)
                            for j=1:length(obj.curves(i).results.s)
                                if strcmp(obj.curves(i).results.s(j).data.id,value)
                                    ndx(i)=true;
                                    break;
                                end
                            end
                        end
                    else
                        vals={strct(:).(fieldname)};
                        if isempty(value)
                            ndx=false(size(vals));
                        elseif ischar(value)
                            ndx=strcmp_w(vals,value);
                        else
                            ndx=strcmp_w(vals,value{1});
                            for i=2:length(value)
                                ndx=ndx|strcmp_w(vals,value{i});
                            end
                        end
                    end
                case {'label','ptype','clabel'}
                    if strcmp(point_or_curve,'curves')
                        vals= {obj.curveprops([obj.curves.propndx]).(fieldname)};
                    elseif strcmp(point_or_curve,'points')
                        vals= {obj.pointprops([obj.points.propndx]).(fieldname)};
                    else %if searching in pointprops or curveprops
                        vals={strct(:).(fieldname)};
                    end
                    if isempty(value)
                        ndx=false(size(vals));
                    elseif ischar(value)
                        ndx=strcmp_w(vals,value);
                    else
                        ndx=strcmp_w(vals,value{1});
                        for i=2:length(value)
                            ndx=ndx|strcmp_w(vals,value{i});
                        end
                    end
                case 'dim'
                    frepars={strct.freepars};
                    codims=cellfun('length',frepars);
                    ndx=codims==value;
                case 'freepars'
                    freepars= {strct(:).freepars};
                    ndx=cellfun(@numel,freepars);
                    ndx=ndx==numel(value);
                    freep= vertcat(freepars{ndx});
                    ndx1=all(freep==repmat(value,size(freep,1),1),2).';
                    ndx(ndx)=ndx1;
                    %                     for i=f
                    %                         ndx(i)=all(freepars{i}==value);
                    %                     end
                    %                     ndx=false(size(strct));
                    %                     for i=1:length(strct)
                    %                         if numel(strct(i).freepars)==numel(value)
                    %                             ndx(i)=all(strct(i).freepars==value);
                    %                         end
                    %                     end
                case 'ctypes'
                    ndx=false(size(strct));
                    for i=1:length(strct)
                        ndx(i)=any(strcmp(strct(i).ctypes,value));
                    end
                case 'samecurve'
                    %works only for curves
                    %same frompoint,same ctype, same freepars same p0
                    %(except freepars)
                    ndx=obj.getsamecurve(value);
                case 'x0+p0'
                    siz=size(strct(1).x0,1)+length(strct(1).p0);
                    vals=zeros(siz,length(strct));
                    for i=1:length(strct)
                        vals(:,i)=[get_x0(strct(i).x0);strct(i).p0];
                    end
                    if iscell(value)
                        value=[get_x0(value{1});value{2}];
                    end
                    differ=sum((vals-repmat(value,1,size(vals,2))).^2,1);
                    if size(strct(1).x0,2)==1
                        ndx=differ<obj.settings.grind.mindist;
                    else
                        %allow more distance for curves
                        ndx=differ<obj.settings.grind.mindist*10;
                    end
                case 'npars'
                    %works only for points
                    if strcmp(point_or_curve,'curves')
                        vals= [obj.curveprops([obj.curves.propndx]).(fieldname)];
                        ndx=vals==value;
                    else
                        ndx=false(size(strct));
                        % labs=regexp({strct.id},'[A-Z]*','match','once');
                        for i=1:length(strct)
                            ctypes=obj.pointprops(strct(i).propndx).ctypes;
                            if ~isempty(ctypes)
                                ndx(i)=any([obj.curveprops(ismember({obj.curveprops.ctype},ctypes)).npars]==value);
                            end
                        end
                    end
                otherwise
                    ndx=[];
                    %error('grind:matcont','Searching fieldname "%s" unknown or not supported', fieldname)
            end
            % ndx=find(ndx);
        end
        
        %% add point to the list of points (g_cont.add_points(findeqs) works to add the equilibria
        function ids=add_points(obj,varargin)
            if ~isempty(varargin)&&length(varargin)<3&&isstruct(varargin{1})
                args.newpoints=varargin{1};
                if length(varargin)>1
                    args.p0=varargin{2};
                end
            else
                args=i_parseargs('newpoints(s#c),id(s#c),x0(n),p0(n),label(s#c),ptype(s#c),data(s),msg(s#c)','if nargs==1,deffields=''x0'';else,deffields=''label,x0'';end;','',varargin);
            end
            if ~isfield(args,'newpoints')
                if ~isfield(args,'x0')
                    args.x0=i_initvar;
                end
                if ~isfield(args,'p0')
                    args.p0=par('-vector');
                end
                if ~isfield(args,'id')&&~isfield(args,'label')&&~isfield(args,'ptype')
                    [~,iseq]=findeq('display','off','maxiter',1,'N0',args.x0,'p0',args.p0);
                    if iseq
                        args.label='EP';
                        args.msg='User defined EP';
                    else
                        args.label='P';
                        args.msg='Initial value';
                    end
                    %                     s2 = sprintf('%g ', round(args.x0 * 10000) / 10000);
                    %                     if ~isfield(args,'msg')
                    %                         args.msg=sprintf('%s: %s',msg,s2);
                    %                     end
                end
                if isfield(args,'ptype')
                    ndx=obj.getndx('pointprops','ptype',args.ptype);
                    args.label=obj.pointprops(ndx).label;
                end
                if ~isfield(args,'msg')
                    args.msg='user defined';
                else
                    args.msg=regexp(args.msg,'[^:]*','match','once');
                end
                if isfield(args,'id')
                    %note that if the id is given no checks are done and it
                    %is assumed that the type is correct
                    ndxp=find(obj.getndx('points','id',args.id),1,'last');
                    obj.points(ndxp)=[];
                    [point,ndx]=create_point(obj,struct('id',args.id,'data',struct(),'msg',args.msg),args.x0,args.p0);
                else
                    if strcmp(args.label,'P')
                        ndxp=find(obj.getndx('points','ptype','P'),1,'last');
                        if ~isempty(ndxp)&&~obj.points_in_use({obj.points(ndxp).id})
                            obj.points(ndxp)=[];
                        end
                    end
                    [point,ndx]=create_point(obj,struct('label',args.label,'data',struct(),'msg',args.msg),args.x0,args.p0);
                end
                if isempty(ndx)
                    obj.points(end+1)=point;
                else
                    point.id=obj.points(ndx).id;
                    obj.points(ndx)=point;
                end
                ids{1}=point.id;
            else
                if isstruct(args.newpoints)&&isfield(args.newpoints,'results')
                    %newpoints is a curve
                    curve=args.newpoints;
                    s=obj.curveprops(curve.propndx).resfun(obj,curve);
                    ids=cell(size(curve.results.s));
                    [ids{:}]=deal('');
                    for i=1:numel(curve.results.s)
                        if ~any(strcmp(curve.results.s(i).label,{'00','99'}))&&...
                                ~(strncmpi(curve.results.s(i).msg,'Neutral saddle',14)||...
                                strcmp(curve.results.s(i).msg,'Zero-Hopf point: neutral saddle'))
                            ndx=curve.results.s(i).index;
                            x0=transpose(s.Y(:,:,ndx));
                            p0=transpose(s.parvalues(ndx,:));
                            [point,ndx1]=create_point(obj,curve.results.s(i),x0,p0);
                            if ~isempty(point)
                                if isempty(ndx1)
                                    obj.points(end+1)=point;
                                else
                                    point.id=obj.points(ndx1).id;
                                    obj.points(ndx1)=point;
                                end
                                ids{i}=point.id;
                            end
                        end
                    end
                else
                    %newpoints are equilibria/initial conditions
                    if ~isfield(args,'p0')
                        args.p0=par('-vector');
                    end
                    %                 if any(strcmp({newpoints.label},'EP'))
                    %                     ndx=obj.getndx('points','ptype','EP');
                    %                     if any(ndx)
                    %                         obj.points=obj.points(~ndx);
                    %                     end
                    %                 end
                    if isfield(args.newpoints,'N0')
                        %if the points are from findeqs, the id field
                        %should be removed as the id needs to be assigned
                        %by the obj
                        args.newpoints=rmfield(args.newpoints, 'id');
                    end
                    if isfield(args.newpoints,'msg')
                        m=regexp({args.newpoints.msg},'[^:]*','match','once');
                        for i=1:length(args.newpoints)
                            args.newpoints(i).msg=m{i};
                        end
                    end
                    
                    ids=cell(size(args.newpoints));
                    for i=1:numel(args.newpoints)
                        [point,ndx]=create_point(obj,args.newpoints(i),args.newpoints(i).N0,args.p0);
                        if isempty(ndx)
                            obj.points(end+1)=point;
                        else
                            point.id=obj.points(ndx).id;
                            obj.points(ndx)=point;
                        end
                        ids{i}=point.id;
                    end
                end
            end
            ids1={obj.points.id};
            [~,ndx]=sort(ids1);
            obj.points=obj.points(ndx);
            %            drawnow;
        end
        %% open the dialog
        function gui(obj,npars)
            obj.settings.derived.ischanged=false;
            if nargin==1
                i_conteq_cont_dlg(obj);
            else
                i_conteq_cont_dlg(obj,npars);
            end
        end
        %% plot the results
        function plot(obj,varargin)
            function freepars=getfreepars(fndx)
                freepars1= {obj.curves(fndx).freepars};
                numpars=cellfun(@numel,freepars1);
                freepars={};
                for i3=1:3
                    freep= vertcat(freepars1{numpars==i3});
                    freep=unique(freep,'rows');
                    for i2=1:size(freep,1)
                        freepars=[freepars,freep(i2,:)];
                    end
                end
            end
            if nargin==1
                freepars=getfreepars(1:length(obj.curves));
                for i1=1:length(freepars)
                    obj.plot('freepars',freepars{i1});
                end
                return;
            end
            args=i_parseargs('fun(s#c),hfig,hax,freepars(i>0#c),curvendx(i),npars(i>0)','npars','',varargin);
            if ~isfield(args,'freepars')
                args.freepars=obj.settings.derived.freepars;
            end
            if isfield(args,'hfig')
                %                 if ishandle(args.hfig)
                %                     delete(args.hfig);
                %                 end
                figure(args.hfig);
            elseif isfield(args,'hax')
                args.hfig=get(args.hax,'Parent');
            else
                args.hfig=figure;
            end
            if ~isfield(args,'hax')
                args.hax=gca;
            end
            if iscell(args.freepars)
                if length(args.freepars)==2
                    freepars(2) =find(strcmp(args.freepars{2},obj.settings.derived.allpars));
                end
                freepars(1) =find(strcmp(args.freepars{1},obj.settings.derived.allpars));
                args.freepars=freepars;
            end
            
            if isfield(args,'curvendx')
                ndx=args.curvendx;
                args.freepars=obj.curves(ndx(1)).freepars;
            elseif isfield(args,'npars')
                ndx=find(getndx(obj,'curves','npars',args.npars));
                if ~isempty(ndx)
                    args.freepars=obj.curves(ndx(1)).freepars(1:args.npars);
                end
            else
                ndx=find(getndx(obj,'curves','freepars',args.freepars));
            end
            crvs=obj.curves(ndx);
            if ~isfield(args,'fun')||isempty(args.fun)
                %select default functions
                if numel(args.freepars)==1
                    v2=i_statevars_names(2);
                    if isempty(v2)
                        args.fun=[obj.settings.derived.allpars(args.freepars) i_statevars_names(1)];
                    else
                        args.fun=[obj.settings.derived.allpars(args.freepars) i_statevars_names(1:2)];
                    end
                elseif ~isempty(args.freepars)
                    args.fun=obj.settings.derived.allpars(args.freepars);
                else
                    args.fun={};
                end
            else
                for j1=1:length(args.fun)
                    ndx2=strcmp('<param1>',args.fun{j1});
                    if any(ndx2)
                        args.fun{j1}(ndx2)=obj.settings.derived.allpars(args.freepars(1));
                    end
                    ndx2=strcmp('<param2>',args.fun{j1});
                    if any(ndx2)
                        args.fun{j1}(ndx2)=obj.settings.derived.allpars(args.freepars(2));
                    end
                end
                if iscell(args.fun{1})
                    ls=cellfun('length',args.fun);
                    if all(ls==1)
                        args.fun=[args.fun{:}];
                    else
                        %not implemented
                        xfun=args.fun{1};
                        yfun=args.fun{2};
                        if numel(args.fun)==3
                            zfun=args.fun{3};
                        else
                            zfun=[];
                        end
                        for i=1:max(ls)
                            if i<=length(xfun)
                                args.fun{1}=xfun{i};
                            else
                                args.fun{1}=xfun{1};
                            end
                            if i<=length(yfun)
                                args.fun{2}=yfun{i};
                            else
                                args.fun{2}=yfun{1};
                            end
                            if ~isempty(zfun)
                                if i<=length(zfun)
                                    args.fun{3}=zfun{i};
                                else
                                    args.fun{3}=zfun{1};
                                end
                            end
                            obj.plot(args);
                            hold on;
                        end
                        return;
                    end
                end
            end
            pnts={};
            for i=1:length(crvs)
                pnts1=cell(size(crvs(i).results.s));
                for j=1:length(crvs(i).results.s)
                    pnts1{j}=crvs(i).results.s(j).data.id;
                end
                pnts=[pnts; pnts1];
            end
            pnts=unique(pnts);
            hax=args.hax;
            hold(hax,'on');
            if numel(args.fun)==3&&~isempty(args.fun{3})
                set(hax, 'View', [8.5 10]);
            else
                set(hax, 'View', [0 90]);
            end
            if ~isoctave&&verLessThan('matlab','8.4.0')
                set(hax, 'drawmode','fast');
            else
                set(hax,'SortMethod','depth');
            end
            i_plotdefaults(args.hfig);
            s=struct('fun',{args.fun},'freepars',args.freepars,'id','','filter',true);
            if isfield(args,'curvendx')
                s.curvendx=args.curvendx;
                [out,xtra]=obj.outfun(s);
                plotstabline(args.hax,out,xtra.stabil,obj.curves(s.curvendx(1)),args.fun);
            else
                for i1=1:length(crvs)
                    s.curvendx=ndx(i1);
                    [out,xtra]=obj.outfun(s);
                    plotstabline(args.hax,out,xtra.stabil,crvs(i1),args.fun);
                end
            end
            % pnts=obj.points();
            s=struct('fun',{args.fun},'freepars',args.freepars,'pointndx',getndx(obj,'points','id',pnts),'filter',true);
            [Y,xtra]=obj.outfun(s);
            if isempty(xtra)
                xtra.ids={};
            end
            plotlabels(args.hfig,Y,xtra.ids)
            hax=args.hax;
            ranpars={obj.settings.derived.parranges.par};
            %set the ranges of the axes if they are parameters with ranges
            if ~isempty(args.fun)
                f=find(strcmp(args.fun{1},ranpars));
            else
                f=[];
            end
            if ~isempty(f)
                ran=obj.settings.derived.parranges(f).range;
                xlims=get(hax,'xlim');
                ran(isnan(ran))=xlims(isnan(ran));
                xlim(args.hax,ran);
            end
            if numel(args.fun)>=2
                f=find(strcmp(args.fun{2},ranpars));
                if ~isempty(f)
                    ran=obj.settings.derived.parranges(f).range;
                    ylims=get(hax,'ylim');
                    ran(isnan(ran))=ylims(isnan(ran));
                    ylim(args.hax,ran);
                end
            end
            if numel(args.fun)>=3
                f=find(strcmp(args.fun{3},ranpars));
                if ~isempty(f)
                    ran=obj.settings.derived.parranges(f).range;
                    zlims=get(hax,'zlim');
                    ran(isnan(ran))=zlims(isnan(ran));
                    zlim(args.hax,ran);
                end
            end
        end
        
        
        
    end
    methods (Access = protected)  %for testing public must become private
        function ndx=getsamecurve(obj,curve)
            %same frompoint,same ctype, same freepars same p0
            %(except freepars)
            %can be overwritten to check other essential parameters
            if isfield(curve,'ctype')
                ndx=obj.getndx('curves','frompoint',curve.frompoint)&obj.getndx('curves','ctype',curve.ctype);
            else
                ndx=obj.getndx('curves','frompoint',curve.frompoint);
            end
            ndx1=obj.getndx('points','id',curve.frompoint);
            if any(ndx1)
                fndx=find(ndx);
                f_p0=obj.points(ndx1).p0;
                mask=true(size(f_p0));
                mask(curve.freepars)=false;
                for i=fndx
                    %              if ndx(i)
                    if length(curve.freepars)==length(obj.curves(i).freepars)&&all(curve.freepars==obj.curves(i).freepars)
                        t_p0=obj.points(obj.getndx('points','id',obj.curves(i).frompoint)).p0;
                        ndx(i)=all(t_p0(mask)==f_p0(mask));
                    else
                        ndx(i)=false;
                    end
                    %            end
                end
            end
        end
        
        function ndx1=all2active(obj,ndx)
            ndx1=ndx;
            for i=1:length(ndx)
                ndx1(i)=sum(obj.settings.grind.activepars(1:ndx(i)));
            end
        end
        function res=translate_option(obj,theoption,toclass)
            %translate a name of a matcont option to coco and vice versa
            opts=obj.get;
            if any(strcmp(theoption,opts(:,1)))
                %same name must be same option
                res=theoption;
                return;
            else
                res='';
            end
            matcont2coco={'MaxStepsize','h_max';...
                'MinStepsize','h_min';...
                'InitStepsize','h0';...
                'MaxCorrIters','ItMX';...
                'MaxNumPoints','PtMX';...
                'FunTolerance','TOL'};
            if isa(obj,'grind_coco')
                if nargin<3
                    toclass='grind_matcont';
                end
                if strcmp(toclass,'grind_matcont')
                    f=strcmp(theoption,matcont2coco(:,1));
                    if any(f)
                        res=matcont2coco{f,2};
                    end
                end
            elseif isa(obj,'grind_matcont')
                if nargin<3
                    toclass='grind_coco';
                end
                if strcmp(toclass,'grind_coco')
                    f=strcmp(theoption,matcont2coco(:,2));
                    if any(f)
                        res=matcont2coco{f,1};
                    end
                end
            end
        end
        %% create vector with all parametere
        function inuse=points_in_use(obj,point_ids) %point_ids should be a cell with ids
            if isempty(obj.curves)
                inuse=false(size(point_ids));
                return;
            end
            usedpnts={obj.curves.frompoint}';
            if isfield(obj.curves(1).results,'s')
                for i=1:length(obj.curves)
                    ids=cell(size(obj.curves(i).results.s));
                    for j=1:length(obj.curves(i).results.s)
                        ids{j}=obj.curves(i).results.s(j).data.id;
                    end
                    usedpnts=[usedpnts; ids];
                end
            end
            % usedpnts=unique(usedpnts);
            inuse=ismember(point_ids,usedpnts);
        end
        function s=allpars_vector(obj)
            pp=regexp(obj.settings.derived.allpars,'[^\(]*','match','once');
            ndx1=false(size(obj.settings.derived.allpars));
            [~,ndx]=unique(regexp(obj.settings.derived.allpars,'[^\(]*','match','once'));
            ndx1(ndx)=true;
            s=transpose(pp(ndx1));
        end
        
        %% filter the parameters and variables
        function s=filter_par_var(obj,s,freepars)
            %apply stateranges
            if ~isempty(s)
                if ~isempty(obj.settings.grind.stateranges)&&any(~isnan(obj.settings.grind.stateranges(:)))
                    ran=obj.settings.grind.stateranges;
                    if size(ran,1)==1
                        %one row is for all state variables
                        if ~isnan(ran(1))
                            if iscell(s.Y)
                                for k=1:numel(s.Y)
                                    s.Y{k}(s.Y{k}<ran(1))=NaN;
                                end
                            else
                                s.Y(s.Y<ran(1))=NaN;
                            end
                        end
                        if  ~isnan(ran(2))
                            if iscell(s.Y)
                                for k=1:numel(s.Y)
                                    s.Y{k}(s.Y{k}>ran(2))=NaN;
                                end
                            else
                                s.Y(s.Y>ran(2))=NaN;
                            end
                        end
                    else
                        %more rows
                        if iscell(s.Y)
                            for k=1:numel(s.Y)
                                yy=s.Y{k};
                                for i=1:size(ran,1)
                                    aa=yy(:,i,:);
                                    aa(~isnan(ran(i,1))&aa<ran(i,1))=NaN;
                                    aa(~isnan(ran(i,2))&aa>ran(i,2))=NaN;
                                    yy(:,i,:)=aa;
                                end
                                s.Y{k}=yy;
                            end
                        else
                            for i=1:size(ran,1)
                                aa=s.Y(:,i,:);
                                aa(~isnan(ran(i,1))&aa<ran(i,1))=NaN;
                                aa(~isnan(ran(i,2))&aa>ran(i,2))=NaN;
                                s.Y(:,i,:)=aa;
                            end
                        end
                    end
                end
                if numel(s.pars)==numel(obj.settings.derived.allpars)
                    fpars=obj.settings.derived.allpars(freepars);
                    ranpars={obj.settings.derived.parranges.par};
                    for i=1:length(freepars)
                        f=find(strcmp(fpars{i},ranpars));
                        if ~isempty(f)
                            ran=obj.settings.derived.parranges(f).range;
                            if ~isnan(ran(1))
                                %makenan_one adds one point at the border
                                %of the parameter range to fill the whole
                                %parameter range
                                s.parvalues(:,freepars(i))=makenan_one(s.parvalues(:,freepars(i)),s.parvalues(:,freepars(i))<ran(1));
                                %  s.parvalues(s.parvalues(:,freepars(i))<ran(1),freepars(i))= NaN;
                            end
                            if ~isnan(ran(2))
                                %  s.parvalues(s.parvalues(:,freepars(i))>ran(2),freepars(i))= NaN;
                                s.parvalues(:,freepars(i))=makenan_one(s.parvalues(:,freepars(i)),s.parvalues(:,freepars(i))>ran(2));
                            end
                        end
                    end
                    
                end
            end
            %          s=struct('parvalues',transpose([pnts.p0]),'pars',{obj.settings.derived.allpars},'Y',permute(Y,[3 2 1]),'t',zeros(size(Y,1),1),'perm',[]);
            %          end
        end
        
        
    end
    
    
end

%% plot line and indicate stability
function plotstabline(hax,Y,stabil,curve,funs)
    function doplot(hax,Y,color,linestyle,displayname)
        if isempty(color)
            color='b';
        end
        if size(Y,2)==1
            plot(hax,Y,'color',color,'linestyle',linestyle,'LineWidth', 1,'tag','conteq','DisplayName',displayname);
        elseif size(Y,2)==2
            plot(hax,Y(:,1),Y(:,2),'color',color,'linestyle',linestyle,'LineWidth', 1,'tag','conteq','DisplayName',displayname);
        else
            h=plot3(hax,Y(:,1),Y(:,2),Y(:,3));
            set(h,'color',color,'linestyle',linestyle,'LineWidth', 1,'tag','conteq','DisplayName',displayname);
        end
    end
%stabil codes:
%1 stable point
%2 stable spiral
%3 saddle point
%4 unstable node
%5 unstable spiral
if ~isempty(funs)
    xlabel(hax,funs{1});
end
if length(funs)>1
    ylabel(hax,funs{2});
end
if length(funs)>2
    zlabel(hax,funs{3});
end
if numel(stabil)~=size(Y,1)
    doplot(hax,Y,curve.color,'-',curve.ctype);
    return;
end
ndx1=stabil==1|stabil==2|stabil==7;
Yres=Y;
Yres(~ndx1,:)=NaN;
if strcmp(curve.ctype,'H')
    doplot(hax,Yres,curve.color,'-',sprintf('supercritical %s',curve.ctype));
else
    doplot(hax,Yres,curve.color,'-',sprintf('stable %s',curve.ctype));
end
ndx1=stabil==3|stabil==6;
Yres=Y;
Yres(~ndx1,:)=NaN;
if strcmp(curve.ctype,'H')
    doplot(hax,Yres,curve.color,':',sprintf('subcritical %s+',curve.ctype));
else
    doplot(hax,Yres,curve.color,':',sprintf('saddle of %s',curve.ctype));
end
ndx1=stabil==4|stabil==5;
Yres=Y;
Yres(~ndx1,:)=NaN;
doplot(hax,Yres,curve.color,'-.',sprintf('unstable %s',curve.ctype));

end
%% plot adjustable lables
function plotlabels(hfig,Y,ids)
%plot labels and adjust to make them visible, optimizing positions after rotating
if iscell(Y)
    sizes=cellfun(@(X)size(X,1),Y);
    ndx=sizes==1;
    Y1=vertcat(Y{ndx});
    ids1=ids(ndx);
    locplotlabels(hfig,Y1,ids1);
    if any(~ndx)
        for i=1:length(Y)
            if ~ndx(i)
                locplotlabels(hfig,Y{i},ids(i));
            end
        end
    end
    i_adjustlabels(hfig,[],'label-text');
    return
else
    locplotlabels(hfig,Y,ids)
    i_adjustlabels(hfig,[],'label-text');
end

    function locplotlabels(hfig,Y,ids)
        hrot=rotate3d(hfig);
        hax=get(hfig,'CurrentAxes');
        ud=get(hfig,'userdata');
        ud.allowedoverlap=0.15; %allowed overlap between the labels, can be made user adaptable
        set(hfig,'userdata',ud);
        hrot.ActionPostCallback = @(hfig,event_obj)i_adjustlabels(hfig,event_obj,'label-text');
        hrot.ActionPreCallback = @(hfig,event_obj)i_resetlabels(hfig,event_obj,'label-text');
        hz= zoom(hfig);
        hz.ActionPostCallback=@(hfig,event_obj)i_adjustlabels(hfig,event_obj,'label-text');
        set(hfig,'ResizeFcn',@(hfig,event_obj)i_adjustlabels(hfig,event_obj,'label-text'));
        
        %offset = (max(Y(:, 1)) - min(Y(:, 1))) * 0.01;
        %labs=regexp(ids,'[A-Z]*','match','once');
        % hopf=strcmp(labs,'H');
        % fold=strcmp(labs,'F');
        if length(ids)==1
            ids1=ids{1};
        else
            ids1=ids;
        end
        if size(Y,2)==2
            %     plot(hax,Y(hopf,1), Y(hopf,2),'r.');
            %     plot(hax,Y(fold,1), Y(fold,2),'b.');
            %     plot(hax,Y(~fold&~hopf,1), Y(~fold&~hopf,2),'k.');
            if length(ids)==size(Y,1)
                h=plot(hax,Y(:,1), Y(:,2),'k.');
            else
                %bifurcation of cycles
                h=plot(hax,Y(:,1), Y(:,2),'r.');
            end
            
            set(h,'tag','conteq', 'DisplayName','special point','MarkerSize',16);
            
            h=text(Y(1:length(ids),1), Y(1:length(ids),2), ids1,'parent',hax);
            set(h,'tag','label-text');
            set(h, 'fontsize', 10);
        elseif size(Y,2)==3
            %     plot3(Y(hopf,1), Y(hopf,2), Y(hopf,3),'r.','parent',hax)
            %     plot3(Y(fold,1), Y(fold,2), Y(fold,3),'b.','parent',hax);
            %     plot3(Y(~fold&~hopf,1), Y(~fold&~hopf,2),Y(~fold&~hopf,3),'k.','parent',hax);
            if length(ids)==size(Y,1)
                plot3(Y(:,1), Y(:,2),Y(:,3),'k.','parent',hax,'tag','conteq', 'DisplayName','special point','MarkerSize',16);
            else
                %bifurcation of cycles
                plot3(Y(:,1), Y(:,2),Y(:,3),'r-','parent',hax,'tag','conteq', 'DisplayName','special point');
            end
            h=text(Y(1:length(ids),1), Y(1:length(ids),2), Y(1:length(ids),3), ids1,'parent',hax);
            set(h,'tag','label-text');
            set(h, 'fontsize', 10);
        end
        %set(hfig,'CurrentAxes',findobj(hfig,'type','axes'));
    end
end
%% replace nan with one?
function res=makenan_one(x,ndx)
res=x;
res(ndx)=NaN;
a1=[0; diff(ndx)]==1;
a2=[diff(ndx); 0]==-1;
res(a1)=x(a1);
res(a2)=x(a2);
end
%%
function res=strcmp_w(vals,value)
%with wildcard * means all
f= strfind(value,'*');
if isempty(f)
    res=strcmp(vals,value);
else
    res=strncmp(regexp(vals,'([A-Za-z+0-9]*_)|([A-Za-z+]*)','match','once'),value,f(1)-1);
end

end
%%
function s=catstruc(s,s1,varargin)
if isempty(s)
    s=s1;
elseif ~isempty(s1)
    f=fieldnames(s1);
    ndx=~ismember(f,varargin);
    f=f(ndx);
    for i=1:length(f)
        x=s1.(f{i});
        if length(size(x))>2
            s.(f{i})=cat(3,s.(f{i}),x);
        else
            s.(f{i})=cat(1,s.(f{i}),x);
        end
    end
end
end
