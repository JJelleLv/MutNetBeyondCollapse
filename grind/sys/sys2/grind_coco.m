classdef grind_coco < grind_cont
    
    
    methods (Access = public)
        function obj = grind_coco
            global g_grind
            %find path for matcont
            if g_grind.solver.haslags
                error('grind:coco:lags','COCO cannot handle delay differential equations')
            end
            if g_grind.solver.isimplicit
                error('grind:coco:dae','COCO cannot handle implicit differential equations (DAE)')
            end
            obj = obj@grind_cont;
            obj.settings.derived.engine_path='';
            obj.settings.derived.engine='coco';
            init_engine(obj);
            init_engine(obj,true);
            verfile = fullfile(obj.settings.derived.engine_path, 'cocoversion.txt');
            if exist(verfile, 'file')
                fid = fopen(verfile, 'r');
                obj.settings.derived.version = fgetl(fid);
                obj.settings.derived.version = obj.settings.derived.version(6:end-4);
                fclose(fid);
            else
                obj.settings.derived.version = 'unknown';
            end
 
            
            crveprops={'EP',1,'EP - Equilibrium points',false,'b';...
                'H',2,'H - Hopf bifurcation',false,'r';...
                'LC',1,'LC - Limit cycle',true,'b';...
                'T',2,'T - Transcritical bifurcation',false,'b';...
                'F',2,'F - Fold bifurcation',false,'g'};
            
            obj.curveprops=struct('ctype',crveprops(:,1),'clabel',crveprops(:,1),'npars',crveprops(:,2),'iscycle',crveprops(:,4),...
                'descr',crveprops(:,3),'color',crveprops(:,5),'resfun',@extract_coco_curve_eq);
            ndx=obj.getndx('curveprops','ctype','F');
            obj.curveprops(ndx).clabel='SN';
            ndx=obj.getndx('curveprops','ctype','T');
            obj.curveprops(ndx).clabel='SN';
            pntprops={'BT','Bogdanov-Takens',2,'';...
                'CP','Cusp',2,'';...
                'EP','Equilibrium',0,'EP';...
                'F','Fold',1,'F';...
                'GH','Generalized Hopf',2,'';...
                'H','Hopf',1,'H,LC';...
                'HH','Double Hopf',2,'';...
                'LC','Limit cycle',0,'';...
                'P','Any initial conditions',0,'EP,LC';...
                'T','Transcritical bif.',1,'EP,T';...
                'ZH','Zero-Hopf',2,''};
            obj.pointprops=struct('ptype',pntprops(:,1),'label',pntprops(:,1),'descr',pntprops(:,2),...
                'codim',pntprops(:,3),'ctypes',pntprops(:,4));
            for i=1:length(obj.pointprops)
                obj.pointprops(i).ctypes=regexp(obj.pointprops(i).ctypes,'[\,]','split');
            end
            ndx=obj.getndx('pointprops','ptype','F');
            obj.pointprops(ndx).label='SN';
            ndx=obj.getndx('pointprops','ptype','T');
            obj.pointprops(ndx).label='SN';
            ndx=obj.getndx('pointprops','ptype','H');
            obj.pointprops(ndx).label='HB';
            N0=i_initvar;
            s2 = sprintf('%g ', round(N0 * 10000) / 10000);
            obj.add_points(struct('id','P1','label','P','msg',sprintf('initial point: %s',s2),'data',[],'N0',N0));
            
        end
        
        function help(obj)
            open(obj);
            open('COCO_ShortRef.pdf')
            close(obj);
        end
        
        
        function [args,found]=get(obj,varargin)
            %get without arguments gives full cell matrix with option names
            %and description and values
            %
            %get('option') gives the result of the current option
            %get('keyopts',pointid,curveid) - gives key options for this
            %combination
            found=true;
            if nargin==1||(nargin==2&&strcmp('-nondefault',varargin{1}))||(nargin==3&&strcmp('-properties',varargin{1}))
                descr=[obj.cont_opt;obj.corr_opt;obj.grind_opt];
                allopts=[[fields(obj.settings.cont);fields(obj.settings.corr);fields(obj.settings.grind)] [struct2cell(obj.settings.cont);struct2cell(obj.settings.corr);struct2cell(obj.settings.grind)]];
                %%% NOTE Here I assume that the options structure and grind_opt have the same order!!
                if nargin==2&&strcmp('-nondefault',varargin{1})
                    args=[descr(:,1),allopts(:,end)];                 
                    for i=size(args,1):-1:1
                        if struccmp(descr{i,4},args{i,end})
                            args(i,:)=[];
                        elseif strcmp(args{i,1},'activepars')&&(islogical(args{i,end})&&all(args{i,end}))
                            args(i,:)=[];                   
                        end
                    end
                    args=args';
                    args=struct(args{:});
                 elseif nargin==3&&strcmp('-properties',varargin{1})
                    ndx=strcmp(descr(:,1),varargin{2});
                    if ~any(ndx)
                        args=[];
                        found=false;
                        return;
                    else
                        args = struct('name',varargin{2},'type',descr{ndx,2},'descr',descr{ndx,3},'default',descr{ndx,4},'value',obj.get(varargin{2}));
                        return;
                    end
                else
                    args=[descr(:,[1 3]),allopts(:,end)];
                    
                    for i=1:size(args,1)
                        if (isnumeric(args{i,end})||islogical(args{i,end}))&& (numel(args{i,end})>1)
                            args{i,end}=mat2str(double(args{i,end}));
                        end
                    end
                end
                return;
            end
            if nargin==2
                opt=varargin{1};
                if isfield(obj.settings.corr,opt)
                    args=obj.settings.corr.(opt);
                elseif isfield(obj.settings.cont,opt)
                    args=obj.settings.cont.(opt);
                elseif isfield(obj.settings.grind,opt)
                    args=obj.settings.grind.(opt);
                elseif isfield(obj.settings.derived,opt)
                    args=obj.settings.derived.(opt);
                else
                    args=[];
                    found=false;
                    if nargout<2
                        fprintf('Unknown option "%s"\n',opt);
                    end
                end
                return;
            end
            if nargin==4 &&any(strcmp(varargin{1},{'alloptions','keyoptions'}))
                if strcmp(varargin{1},'alloptions')
                    args=obj.get;
                    return;
                end
                funargs=obj.get_func_args(varargin{2},varargin{3});
                if isempty(funargs)
                    found=false;
                end
                f=[fieldnames(obj.settings.cont);fieldnames(obj.settings.corr);fieldnames(obj.settings.grind)];
                ndx=ismember(funargs,f);
                keyopts=[funargs(ndx) {'relh'}];
                args=obj.get;
                ndx=ismember(args(:,1),keyopts);
                args=args(ndx,:);
            end
        end
        function args=set(obj,varargin)
            %% set options
            if nargin==1
                %without arguments a list of current settings is given
                disp('Current options:');
                f1={'grind','corr','cont'};
                for j=1:length(f1)
                    nsiz=18;
                    f=fieldnames(obj.settings.(f1{j}));
                    for i=1:length(f)
                        if isempty(obj.settings.(f1{j}).(f{i}))
                            fprintf('%s%s []\n',f{i},repmat(' ',1,nsiz-length(f{i})));
                        elseif ischar(obj.settings.(f1{j}).(f{i}))
                            fprintf('%s%s ''%s''\n',f{i},repmat(' ',1,nsiz-length(f{i})),obj.settings.(f1{j}).(f{i}));
                        elseif iscell(obj.settings.(f1{j}).(f{i}))
                            fprintf('%s%s %dx%d cell array\n',f{i},repmat(' ',1,nsiz-length(f{i})),size(obj.settings.(f1{j}).(f{i})));
                        else
                            fprintf('%s%s %s\n',f{i},repmat(' ',1,nsiz-length(f{i})),mat2str(obj.settings.(f1{j}).(f{i})));
                        end
                    end
                end
                return;
            end
%             pairs1=obj.corr_opt(:,1);
%             for i=1:length(pairs1)
%                 if ~isempty(obj.corr_opt(i,2))
%                     pairs1{i}=sprintf('%s(%s)',pairs1{i},obj.corr_opt{i,2});
%                 end
%             end
%             pairs2=obj.cont_opt(:,1);
%             for i=1:length(pairs2)
%                 if ~isempty(obj.cont_opt(i,2))
%                     pairs2{i}=sprintf('%s(%s)',pairs2{i},obj.cont_opt{i,2});
%                 end
%             end
%             pairs3=obj.grind_opt(:,1);
%             for i=1:length(pairs3)
%                 if ~isempty(obj.grind_opt(i,2))
%                     pairs3{i}=sprintf('%s(%s)',pairs3{i},obj.grind_opt{i,2});
%                 end
%             end
%             pairs4={'frompoint(c#s)';'ctype(c#s)';'file(s)'};
             ext_opts={'frompoint','c#s','code of the point from which the continuation starts (for instance ''EP1'') (default: empty)','';...
                        'ctype','c#s','type of curve for continuations (for instance EP)','';...
                         'file','s','open a previously saved conteq session',''};
            args= i_parseargs([obj.corr_opt;obj.cont_opt;obj.grind_opt;ext_opts].','if nargs==1,deffields=''file'';else,deffields=''par1,frompoint,ctype'';end;','-defaults,-coco,-matcont,-c,-out,-l',varargin);
            if any(strcmp(args.opts,'-defaults'))
                obj.setdefaults;
            end
            if isfield(args,'frompoint')  %if we only define frompoint ctype should be emptied
                obj.settings.derived.ctype='';
            end
            if isfield(args,'activepars')
                if length(args.activepars)>length(obj.settings.derived.allpars)
                    args.activepars=args.activepars(1:length(obj.settings.derived.allpars));
                elseif length(args.activepars)<length(obj.settings.derived.allpars)
                    actpar=args.activepars;
                    args.activepars=obj.settings.grind.activepars;
                    args.activepars(1:length(actpar))=actpar(:);
                end
            end
            
            args=rmfield(args,'opts');
            f=fieldnames(args);
            grindfields=obj.grind_opt(:,1);
            corrfields=obj.corr_opt(:,1);
            derfields={'frompoint','ctype'};
            for i=1:length(f)
                if any(strcmp(f{i},grindfields))
                    obj.settings.grind.(f{i})=args.(f{i});
                elseif any(strcmp(f{i},corrfields))
                    obj.settings.corr.(f{i})=args.(f{i});
                elseif any(strcmp(f{i},derfields))
                    obj.settings.derived.(f{i})=args.(f{i});
                else
                    obj.settings.cont.(f{i})=args.(f{i});
                end
            end
            if ~isempty(obj.settings.grind.par2)&&strcmp(obj.settings.grind.par1,obj.settings.grind.par2)
                obj.settings.grind.par2='';
                error('grind:coco','The second free parameter must be different from the first')
            end
            siz=size(obj.settings.grind.stateranges);
            if siz(1)==2&&siz(2)~=2
                obj.settings.grind.stateranges=tranpose(obj.settings.grind.stateranges);
                siz=size(obj.settings.grind.stateranges);
            end
            if siz(1)~=1&&siz(1)~=obj.settings.derived.ndim||siz(2)~=2
                error('grind:matcont:stateranges','Size of "stateranges" is not correct (either one row or a row for each state)');
            end
            %update derived options
            %obj.settings.derived.allpars; (allways all parameters)
            %obj.settings.derived.freepars update with par1 and par2
            %update with parranges
            %obj.settings.derived.parranges=struct('par',{},'range',[NaN NaN]);
            %not yet implemented
            %obj.settings.derived.varranges=struct('varno',{},'range',[NaN NaN]);
            ipar1=find(strcmp(obj.settings.grind.par1, obj.settings.derived.allpars),1);
            if ~obj.settings.grind.activepars(ipar1)
                warning('grindcoco:notactive','parameter %s is not active',obj.settings.grind.par1)
                ipar1=find(obj.settings.grind.activepars,1);
                args.par1=obj.settings.derived.allpars{ipar1};
                obj.settings.grind.par1=args.par1;
            end
            ipar2=find(strcmp(obj.settings.grind.par2, obj.settings.derived.allpars),1);
            if ~isempty(ipar2)&&~obj.settings.grind.activepars(ipar2)
                warning('grindcoco:notactive','parameter %s is not active',obj.settings.grind.par2)
                pars=obj.settings.grind.activepars;
                pars(ipar1)=false;
                ipar2=find(pars,1);
                args.par2=obj.settings.derived.allpars{ipar2};
                obj.settings.grind.par2=args.par2;
            elseif isempty(ipar2)&&~isempty(obj.settings.grind.par2)
                error('grind:coco:nopar','%s is not a parameter',args.par2)
            end
            obj.settings.derived.freepars=[ipar1,ipar2];
            pars=obj.settings.derived.allpars(obj.settings.derived.freepars);
            existpars={obj.settings.derived.parranges.par};
            for i=1:length(pars)
                f= find(strcmp(pars{i},existpars));
                if isempty(f)
                    f=length(obj.settings.derived.parranges)+1;
                end
                if i==1
                    obj.settings.derived.parranges(f)=struct('par',pars{i},'range',obj.settings.grind.parranges1);
                elseif i==2
                    obj.settings.derived.parranges(f)=struct('par',pars{i},'range',obj.settings.grind.parranges2);
                end
                
            end
            
        end
        
        
        
        function saveas(obj,filename)
            if nargin<2
                filename='';
            end
            [~,~,ext]=fileparts(filename);
            if ~strcmp(ext,'.m')
                prob=obj.settings.derived.prob;
                obj.settings.derived.prob=[];
                path=cd;
                obj.open;
                warning off MATLAB:dispatcher:UnresolvedFunctionHandle
                for i=1:length(obj.curves)
                    for k=1:length(obj.curves(i).data.runids)
                        cd(obj.curves(i).data.runids{k});
                        dd=dir('*.mat');
                        for j=1:length(dd)
                            [~,fname]=fileparts(dd(j).name);
                            obj.curves(i).data.files{k}.(fname)=load(dd(j).name);
                        end
                        dd=dir('*.txt');
                        for j=1:length(dd)
                            [~,fname]=fileparts(dd(j).name);
                            fid=fopen(dd(j).name,'r');
                            obj.curves(i).data.files{k}.(fname)=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                            fclose(fid);
                        end
                    end
                end
                cd(path);
            end
            saveas@grind_cont(obj,filename);
            if isempty(obj.settings.derived.prob)
                for i=1:length(obj.curves)
                    if isfield(obj.curves(i).data,'files')
                        obj.curves(i).data=rmfield(obj.curves(i).data,'files');
                    end
                end  
                obj.settings.derived.prob=prob;
                warning on MATLAB:dispatcher:UnresolvedFunctionHandle
            end
        end
        
        
        function out=handles(obj)
            global  g_grind;
            %  Matcont style:
            %             out{2} = ode_handle;
            %             out{3} = jac_handle;
            %             out{4} = jacp_handle;
            %             out{4} = hessian;
            out{9}=[];
            if g_grind.statevars.vector
                if all(obj.settings.grind.activepars)
                    out{2} = i_getodehandle('coco', obj.settings.derived.allpars);
                    out{3} = [];%i_getodehandle('coco_jac', obj.settings.derived.allpars, 'numonly');
                    out{4} = [];%i_getodehandle('coco_jacp', obj.settings.derived.allpars, 'numonly');
                else
                    out{2} = i_getodehandle('coco', obj.settings.derived.allpars(obj.settings.grind.activepars));
                    out{3} = [];%i_getodehandle('coco_jac', obj.settings.derived.allpars(obj.settings.grind.activepars), 'numonly');
                    out{4} = [];%i_getodehandle('coco_jacp', obj.settings.derived.allpars(obj.settings.grind.activepars), 'numonly');
                end
                
                return;
            end
            out{2} = i_getodehandle('coco', obj.settings.derived.allpars(obj.settings.grind.activepars));
            if ~obj.settings.grind.symbolic
                out{3} = i_getodehandle('coco_jac', obj.settings.derived.allpars(obj.settings.grind.activepars), 'numonly');
                out{4} = i_getodehandle('coco_jacp', obj.settings.derived.allpars(obj.settings.grind.activepars), 'numonly');
            else
                out{3} = i_getodehandle('coco_jac', obj.settings.derived.allpars(obj.settings.grind.activepars));
                out{4} = i_getodehandle('coco_jacp', obj.settings.derived.allpars(obj.settings.grind.activepars));
            end
            allpars=sprintf(',%s',obj.settings.derived.allpars{obj.settings.grind.activepars});
            if ~isempty(g_grind.syms.Hessian)
                out{5} = i_getodehandle('Hessian',allpars);
            end
        end
        
    end
    methods (Static, Access = public)
        
        function newpars = translatepars(pars, allpars)
            if nargin == 1
                allpars = {pars};
            end
            
            newpars = pars;
            if ischar(pars)
                cocoidentifiers={'PT','StepSize','TIME','SLAB','LAB','TYPE', 'x'};
                if any(strcmp(pars, cocoidentifiers))
                    i = 1;
                    newpars = sprintf('%s%d', pars, i);
                    while any(strcmp(newpars, allpars))
                        i = i + 1;
                        newpars = sprintf('%s%d', pars, i);
                    end
                    
                end
                
            else
                for i = 1:length(pars)
                    newpars{i} = grind_coco.translatepars(pars{i}, allpars);
                end
                
            end
            
        end
        
        function  argnames=get_func_args(label,clabel)
            
            init_functs={'BP_BP',{'odefile','x','p','ap','bp'};...
                'BPC_BPC',{'odefile','x','s','ap','ntst','NCOL','bp'};...
                'BP_EP',{'odefile','x','p','s','amp'};...
                'EP_EP',{'odefile','x','p','ap'};...
                'H_EP',{'odefile','x','p','ap'};...
                'LP_EP',{'odefile','x','p','ap','varargin'};...
                'HTHet_Het',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'Het_Het',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'BT_Hom',{'odefile','x','s','p','ap','NTST','NCOL','TTolerance ','amplitude','extravec'};...
                'HSN_Hom',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'HTHom_Hom',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'Hom_Hom',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'LC_Hom',{'odefile','x','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'NCH_Hom',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'HSN_HSN',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'HTHSN_HSN',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'Hom_HSN',{'odefile','x','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'LC_HSN',{'odefile','x','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'NCH_HSN',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                'HTHet_HTHet',{'odefile','x','v','s','p','ap','up','aup','sp','asp','NTST','NCOL','T','eps1','eps1tol'};...
                'HTHom_HTHom',{'odefile','x','v','s','p','ap','up','aup','sp','asp','NTST','NCOL','T','eps1','eps1tol'};...
                'HTHSN_HTHSN',{'odefile','x','v','s','p','up','aup','sp','asp','NTST','NCOL','T','eps1','eps1tol'};...
                'BT_H',{'odefile','x','p','ap'};...
                'GH_H',{'odefile','x','p','ap'};...
                'HH_H',{'odefile','x','p','ap'};...
                'H_H',{'odefile','x','p','ap'};...
                'ZH_H',{'odefile','x','p','ap'};...
                'BPC_LC',{'odefile','x','v','s','NTST','NCOL','h'};...
                'H_LC',{'odefile','x','p','ap','amp','NTST','NCOL'};...
                'LC_LC',{'odefile','x','v','s','par','ap','NTST','NCOL'};...
                'LPC_LC',{'odefile','x','v','s','par','ap','NTST','NCOL'};...
                'NS_LC',{'odefile','x','v','s','par','ap','NTST','NCOL'};...
                'P_EP',{'odefile','x','p','ap','ndays'};...
                'P_LC',{'odefile','x','p','ap','NTST','NCOL','ndays','cycletol'};...
                'PD_LC',{'odefile','x','s','NTST','NCOL','h'};...
                'PD_LC2',{'odefile','x','s','ap','NTST','NCOL','h'};...
                'BP_LP',{'odefile','x','p','ap'};...
                'BT_LP',{'odefile','x','p','ap'};...
                'CP_LP',{'odefile','x','p','ap'};...
                'LP_LP',{'odefile','x','p','ap','varargin'};...
                'ZH_LP',{'odefile','x','p','ap'};...
                'BPC_LPC',{'odefile','x','s','ap','NTST','NCOL'};...
                'CPC_LPC',{'odefile','x','s','ap','NTST','NCOL'};...
                'GH_LPC',{'odefile','x','p','s','ap','NTST','NCOL','eps','varargin'};...
                'LPC_LPC',{'odefile','x','s','ap','NTST','NCOL','varargin'};...
                'LPNS_LPC',{'odefile','x','s','ap','NTST','NCOL'};...
                'LPPD_LPC',{'odefile','x','s','ap','NTST','NCOL'};...
                'R1_LPC',{'odefile','x','s','ap','NTST','NCOL'};...
                'CHNS_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'HH_NS1',{'odefile','x','p','s','ap','NTST','NCOL','eps'};...
                'HH_NS2',{'odefile','x','p','s','ap','NTST','NCOL','eps'};...
                'LPNS_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'NS_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'PDNS_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'R1_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'R2_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'R3_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'R4_NS',{'odefile','x','s','ap','NTST','NCOL'};...
                'ZH_NS',{'odefile','x','p','s','ap','NTST','NCOL','eps'};...
                'GPD_PD',{'odefile','x','s','ap','NTST','NCOL'};...
                'LPPD_PD',{'odefile','x','s','ap','NTST','NCOL'};...
                'PDNS_PD',{'odefile','x','s','ap','NTST','NCOL'};...
                'PD_PD',{'odefile','x','s','ap','NTST','NCOL'};...
                'R2_PD',{'odefile','x','s','ap','NTST','NCOL'}};
            ndx=strcmp([label '_' clabel],init_functs(:,1));
            argnames=init_functs(ndx,2);
            if ~isempty(argnames)
                argnames=argnames{1};
            end
        end
        function engine_path=findupdate(checkupdates, url)
            function writecfg(engine_path)
                if exist('grind.cfg','file')
                    fid=fopen('grind.cfg','r');
                    lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                    lines=lines{1};
                    fclose(fid);
                    ndx=~cellfun('isempty',(strfind(lines,'cocodir=')));
                    if ~any(ndx)
                        ndx=length(lines)+1;
                    end
                    lines{ndx}=sprintf('cocodir=%s',engine_path);
                else
                    lines={sprintf('cocodir=%s',engine_path)};
                end
                fid=fopen(fullfile(grindpath,'grind.cfg'),'w');
                fprintf(fid,'%s\n',lines{:});
                fclose(fid);
            end
            %find COCO on the computer
            %if not available install it from the url
            %add coco toolbox directories to the search path
            %
            %if doupdate is true it checks whether there is a newer version available
            if nargin == 0
                checkupdates = [];
            end
            engine_path='';
            if nargin < 2
                url = 'http://content.alterra.wur.nl/webdocs/internet/aew/downloads/';
            end
            if exist('grind.cfg','file')
                fid=fopen('grind.cfg','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                ndx=~cellfun('isempty',(strfind(lines,'cocodir=')));
                if any(ndx)
                    line=lines{ndx};
                    engine_path=regexp(line,'(?<=[=]).*','match', 'once');
                    if ~exist(engine_path,'dir')||~exist(fullfile(engine_path,'core','toolbox','coco.m'),'file')
                        engine_path='';
                    end
                end
            else
                cocotoolboxdir = which('coco.m');
                if ~isempty(cocotoolboxdir)
                    f=strfind(cocotoolboxdir,fullfile('core','toolbox','coco.m'));
                    if ~isempty(f)
                        engine_path = cocotoolboxdir(1:f(end) - 2);
                    end
                end
            end
            if isempty(engine_path)
                cocotoolboxdir = findgrindfile('coco.m');
                if ~isempty(cocotoolboxdir)
                    f=strfind(cocotoolboxdir,fullfile('core','toolbox'));
                    if ~isempty(f)
                        engine_path = cocotoolboxdir(1:f(end) - 2);
                        writecfg(engine_path)
                    end
                end
            end
            
            if ~isempty(engine_path)
                cocofound = true;
            else
                %if ~checkupdates
                %    msgbox('COCO (Continuation Core and Toolboxes) not available, installing last version','GRIND depends on COCO');
                %end
                if isempty(checkupdates)
                    checkupdates = true;
                end
                cocofound = false;
            end
            
            if ~isempty(checkupdates)&& checkupdates
                [newversion, status] = urlread([url 'cocoversion.txt']);
                if newversion(end) + 0 == 10
                    newversion = newversion(1:end - 1);
                end
                
                if status == 1
                    if ~isempty(engine_path)
                        verfile = fullfile(engine_path, 'cocoversion.txt');
                        if exist(verfile, 'file')
                            fid = fopen(verfile, 'r');
                            oldversion = fgetl(fid);
                            fclose(fid);
                        else
                            oldversion = '';
                        end
                        
                        if strcmp(oldversion, newversion)
                            disp('No newer version of COCO (Continuation Core and Toolboxes) available');
                            return;
                        else
                            if ~isempty(oldversion)|| strcmp( questdlg('OK to remove old version of COCO? (recommended)', ...
                                    'Remove coco', 'Yes', 'No', 'Yes'),'Yes')
                                disp('Removing old version of COCO');
                                warning('off','MATLAB:RMDIR:RemovedFromPath');
                                rmdir(engine_path, 's');
                                warning('on','MATLAB:RMDIR:RemovedFromPath');
                                oldversion = 'removed';
                            end
                            
                        end
                        
                    else
                        oldversion = '';
                    end
                    
                    if isempty(engine_path)
                        engine_path = grindpath;
                        f = strfind(engine_path, 'grind');
                        if ~isempty(f)
                            engine_path = engine_path(1:f(1) - 2);
                        end
                        
                    else
                        f = strfind(engine_path, 'coco');
                        if ~isempty(f)
                            engine_path = engine_path(1:f(1) - 2);
                        end
                        
                    end
                    
                    fprintf('Downloading %s\n', newversion)
                    unzip([url newversion], engine_path);
                    engine_path = fullfile(engine_path, 'coco');
                    if ~isempty(oldversion)||~cocofound
                        verfile = fullfile(engine_path, 'cocoversion.txt');
                        fid = fopen(verfile, 'w');
                        fprintf(fid, '%s\n', newversion);
                        fclose(fid);
                    end
                    writecfg(engine_path)
                    disp('Updated Continuation Core and Toolboxes (COCO)');
                end
                
            end
            
        end  %method
    end
    methods (Access = protected)  %for testing public must become private
        
        
        function parrange = makerange(obj,parrange,freepar)
            if nargin<3
                freepar=obj.settings.derived.freepars(1);
            end
            parrange = parrange(~isnan(parrange));
            currval = evalin('base', obj.settings.derived.allpars{freepar});
            if isempty(parrange)
                parrange = [-1E100 1E100];
            elseif length(parrange) == 1
                if parrange < currval
                    parrange = [parrange 1E100];
                else
                    parrange = [-1E100, parrange];
                end
            elseif length(parrange) == 2
                if parrange(1)>parrange(2)
                    parrange=flipud(parrange);
                end
                if parrange(1)>currval
                    parrange(1)=currval;
                elseif parrange(2)<currval
                    parrange(2)=currval;
                end
            end
        end
        function extract_curve_files(obj)
            global g_grind;
            No=0;
            warning off MATLAB:dispatcher:UnresolvedFunctionHandle
            olddir=cd;
            for i=1:length(obj.curves)
                if isfield(obj.curves(i).data,'files')
                    for k=1:length(obj.curves(i).data.runids)
                        No=No+1;
                        therunid=fullfile(grindpath, 'tmp',g_grind.odefile,int2str(No));
                        obj.curves(i).data.runids{k}=therunid;
                        dir_success = mkdir(therunid);
                        if dir_success
                            cd(therunid);
                            fnames=fieldnames(obj.curves(i).data.files{k});
                            for j=1:length(fnames)
                                if isstruct(obj.curves(i).data.files{k}.(fnames{j}))
                                    g_s=obj.curves(i).data.files{k}.(fnames{j}); %#ok<NASGU>
                                    save(fnames{j},'-struct','g_s');
                                else
                                    fid=fopen([fnames{j} '.txt'],'w');
                                    fprintf(fid,'%s\n',obj.curves(i).data.files{k}.(fnames{j}){1}{:});
                                    fclose(fid);
                                end
                            end
                        end
                    end
                    obj.curves(i).data=rmfield(obj.curves(i).data,'files');
                end
            end
            cd(olddir);
           % warning on MATLAB:dispatcher:UnresolvedFunctionHandle
        end
        function setdefaults(obj)
            global g_grind;
            sett=obj.corr_opt(:,[1 4]).';
            obj.settings.corr=struct(sett{:});
            sett=obj.cont_opt(:,[1 4]).';
            obj.settings.cont=struct(sett{:});
            %             obj.settings.matcont.InitStepsize='0.01*relh*range';
            %             obj.settings.matcont.MinStepsize='0.001*relh*range';
            %             obj.settings.matcont.MaxStepsize='relh*range';
            sett=obj.grind_opt(:,[1 4]).';
            obj.settings.grind=struct(sett{:});
            if g_grind.statevars.vector
                pars = {};
                for i = 1:length(g_grind.pars)
                    siz= evalin('base',sprintf('size(%s);',g_grind.pars{i}));
                    elems = allelems(g_grind.pars{i}, siz);
                    pars = [pars ,transpose(elems(:))]; %#ok<AGROW>
                end
                
            else
                pars = g_grind.pars;
            end
            obj.settings.grind.symbolic=~g_grind.statevars.vector;
            obj.settings.grind.activepars=true(size(pars));
            obj.settings.derived.allpars=pars;
            obj.settings.derived.direction=[true true];
            obj.settings.derived.freepars=[];
            obj.settings.derived.prob=[];
            obj.settings.derived.parranges=struct('par',{},'range',[NaN NaN]);
        end
        
        function run_point(obj,frompoint,ctype,varargin)
            global g_grind;
            %              ndx=obj.getndx('curveprops','ctype',ctype);
            %             if any(strcmp(ctype,{'F','T'}))
            %                 ctype='SN';
            %             end
            fprintf(['\nRunning <a href="https://sourceforge.net/projects/cocotools/">COCO[%s] (Continuation Core and Toolboxes)</a> ' ...
                '\n(Dankowicz, H. & Schilder, F. 2013 Recipes for Continuation. \nSIAM, Philadelphia. 584 pp.)\n'],obj.settings.derived.version);
            extract_curve_files(obj)
            if obj.settings.grind.symbolic&&(isempty(g_grind.syms.Jacobian)||isempty(g_grind.syms.Jacobianp))
                try
                    if i_hastoolbox('symbolic')&&~g_grind.statevars.vector
                        enterjac('-s');
                    end
                catch %#ok<CTCH>
                    obj.settings.grind.symbolic = false;
                    warning('conteq:symfail', 'Cannot determine symbolic differentials, using numeric approximations instead');
                    g_grind.syms.Jacobian = {};
                    g_grind.syms.Jacobianp = {};
                end
            end
            npars= obj.curveprops(obj.getndx('curveprops','ctype',ctype)).npars;
            if numel(obj.settings.derived.freepars)<npars
                error('grind:conteq:freepars','Not enough parameters selected: curve "%s" needs %d free parameters',ctype,npars);
            end
            obj.updateprob(true);
            No = length(obj.curves) + 1;
            therunid=fullfile(grindpath, 'tmp',g_grind.odefile,int2str(No));
            
            contpar1 = obj.translatepars(obj.settings.derived.allpars(obj.settings.derived.freepars(1:npars)), obj.settings.derived.allpars);
            cocopars = obj.translatepars(obj.settings.derived.allpars(obj.settings.grind.activepars), obj.settings.derived.allpars);
            %[ode_handle, jac_handle, jacp_handle]
            hands= obj.handles;
            numjac= isempty(hands{3})||strcontains(func2str(hands{3}),')coco_num_');
            numjacp= isempty(hands{4})||~strcontains(func2str(hands{4}),')checknan(');
            if numjac&&numjacp
                disp('All Jacobians evaluated numerically');
            elseif numjac
                disp('Jacobian evaluated numerically');
            elseif numjacp
                disp('Jacobianp evaluated numerically');
            end
            
            parrange = obj.makerange(obj.settings.grind.parranges1);
            findingGH=false;
            % label=obj.pointprops(frompoint.propndx).label;
            id=obj.pointprops(frompoint.propndx).ptype;
            action=sprintf('%s_%s',id,ctype);
            switch action
                case {'EP_EP','P_EP'}
                    if isempty(frompoint.x0)
                        frompoint.x0 = i_initvar;
                    end
                    frompoint.x0(frompoint.x0 == 0)=1E-10;
                    if ~isempty(hands{4})
                        try
                            Jacp = hands{4}(frompoint.x0, frompoint.p0(obj.settings.grind.activepars));
                        catch %#ok<CTCH>
                            hands{4} = i_getodehandle(12, obj.settings.derived.allpars, 'numonly');
                            Jacp = hands{4}(frompoint.x0, frompoint.p0(obj.settings.grind.activepars));
                        end
                        
                        if ~any(isnan(Jacp(:)))
                            x0 = zeros(size(frompoint.x0));
                            if any(isnan((hands{2}(x0, frompoint.p0(obj.settings.grind.activepars)))))
                                warning('grind:conteq','The model is not defined if state variables are zero')
                                x0 = x0 + 1E-10;
                            end
                            Jacp = hands{4}(x0, frompoint.p0(obj.settings.grind.activepars));
                        end
                        
                        if any(isnan(Jacp(:)))
                            warning('grind:conteq','Jacp has NaN values, removed some active parameters, alternatively you can set "symbolic" false')
                            activepars = obj.settings.grind.activepars & ~any(isnan(Jacp), 1);
                            if ~activepars(obj.settings.derived.freepars(1))
                                obj.settings.grind.symbolic = false;
                            else
                                obj.settings.grind.activepars = activepars;
                            end
                            hands = obj.handles;
                            contpar1 = obj.translatepars(obj.settings.derived.allpars(obj.settings.derived.freepars(1)), obj.settings.derived.allpars);
                            cocopars = obj.translatepars(obj.settings.derived.allpars(obj.settings.grind.activepars), obj.settings.derived.allpars);
                        end
                    end
                    obj.settings.derived.prob = ode_isol2ep(obj.settings.derived.prob, '', hands{2:4}, frompoint.x0, cocopars, frompoint.p0(obj.settings.grind.activepars));
                case 'T_EP'
                    obj.settings.derived.prob = ode_BP2ep(obj.settings.derived.prob,'',frompoint.data.runid, frompoint.data.labnrBP);
                case 'H_LC'
                    obj.settings.derived.prob = ode_HB2po(obj.settings.derived.prob,'',frompoint.data.runid, frompoint.data.labnr);
                case 'P_LC'  %needs to be adapted
                    ispo = true;
                    obj.settings.derived.prob = ode_isol2po(obj.settings.derived.prob, '', hands{2:4}, obj.settings.IO.t, obj.settings.IO.Y,cocopars, obj.settings.p0(obj.settings.grind.activepars));
                case 'H_H'
                    obj.settings.derived.prob = ode_HB2HB(obj.settings.derived.prob,'',frompoint.data.runid, frompoint.data.labnr);
                    findingGH = obj.addfuncsHB;
                case 'F_F'
                    obj.settings.derived.prob = ode_SN2SN(obj.settings.derived.prob,'',frompoint.data.runid, frompoint.data.labnr);
                case 'T_T'
                    obj.settings.derived.prob = ode_SN2SN(obj.settings.derived.prob,'',frompoint.data.runid, frompoint.data.labnr);
                otherwise
                    error('grind:conteq:unknownlab','Unknown label %s',action);
            end
            if ~all(isnan(obj.settings.grind.stateranges))
                if ~isvalidrange(frompoint.x0, obj.settings.grind.stateranges)
                    error('grind:conteq','Cannot start continuation outside the defined range for state variables');
                end
                
                [data, uidx] = coco_get_func_data(obj.settings.derived.prob, 'ep', 'data', 'uidx');
                %                    minstates=obj.settings.N0ranges(:,1);%this can become adjustable mask with minimal values (0) or NAN for state variables
                %     maxstates=obj.settings.N0ranges(:,2);
                obj.settings.derived.prob = coco_add_func(obj.settings.derived.prob, 'coco_state_boundary', @(prob, data, u)coco_state_boundary(prob,data,u,obj.settings.grind.stateranges), data.ep_eqn, ...
                    'regular', 'x.min', 'uidx', uidx);
                obj.settings.derived.prob = coco_add_event(obj.settings.derived.prob, 'EPS', 'boundary','x.min', 0);
            end
            
            %run continuation problem
            %             bd3 = coco(obj.prob,runid_bif, [], contpars, makerange(obj.settings.parranges{1})); % cont toolbox arguments
            
            bd1 = coco(obj.settings.derived.prob,therunid, [],1, contpar1, parrange); % cont toolbox arguments
            i_waitbar(1);
            %save the results in the same way as in MATCONT
            curve=struct('freepars',obj.settings.derived.freepars(1:npars),'ctype',ctype,'frompoint',frompoint.id,...
                'data',[],'color',[],'propndx',find(obj.getndx('curveprops','ctype',ctype)),'results',struct('stabil',[]));
            curve.color=obj.curveprops(curve.propndx).color;
            curve.data.haslyapunov=findingGH;
            curve.data.runids={therunid};
            curve.data.settings=obj.get('-nondefault');
            curve=obj.extract_curve_data(bd1,curve,cocopars);
            
            add_curve(obj,curve);
            
        end
        
        function curve=extract_curve_data(obj,bd1,curve,cocopars)
            %replace x as it may lead to confusion;
            bd1{1,strcmp(bd1(1,:),'x')}='results.x';
            therunid=curve.data.runids{end};
            fcol=transpose(coco_bd_col(bd1,cocopars(obj.all2active(curve.freepars)),true));
            fcol=num2cell(fcol,2);
            fcol=[{'results.freepars'};fcol];
            if size(fcol,1)<size(bd1,1)
                fcol(size(fcol,1)+1:size(bd1,1),1)={[]};
            end
            bd1=[bd1,fcol];
            parndx=ismember( bd1(1,:),cocopars);
            curve.results.bd=bd1(:,~parndx);
            frompoint=obj.points(obj.getndx('points','id',curve.frompoint));
            p0=frompoint.p0;
            labnr=coco_bd_col(bd1,'LAB',false);
            types=coco_bd_col(bd1,'TYPE',false);
            x=coco_bd_col(bd1,'results.x',true);
            fpars=coco_bd_col(bd1,'results.freepars',true);
            types{1}='00';
            types{end}='99';
            ndx=~cellfun('isempty',types);
            indexs=double(ndx);
            indexs(ndx)=find(ndx);
            ndx=~(cellfun('isempty',types)|strcmp(types,'RO')|strcmp(types,'EP'));
            curve.results.s=struct('index',num2cell(indexs(ndx)),'label',types(ndx),'msg','','data',[]);
            labnr=labnr(ndx);
            for i=1:length(curve.results.s)
                switch curve.results.s(i).label
                    case 'HB'
                        parslist = sprintf(',%s',obj.settings.derived.allpars{:});
                        if obj.settings.grind.symbolic
                            odefun = i_getodehandle('singlestep', parslist);
                            jac = i_getodehandle('Jacobian', parslist);
                            hess = i_getodehandle('Hessian', parslist);
                        else
                            odefun = i_getodehandle('singlestep', parslist,'numonly');
                            jac = i_getodehandle('Jacobian', parslist,'numonly');
                            hess = i_getodehandle('Hessian', parslist,'numonly');
                        end
                        x0 = x(1:obj.settings.derived.ndim,curve.results.s(i).index);
                        p0(curve.freepars) = fpars(:,curve.results.s(i).index);
                        p0C=num2cell(p0);
                        try
                            [curve.results.s(i).data.lyapunov, curve.results.s(i).data.matcontlyap] = firstlyap(0, x0, odefun, jac, hess, p0C{:});
                        catch %#ok<CTCH>
                            warning('grind:cont:nolyap','Cannot calculate Lyapunov coefficient')
                            curve.results.s(i).data.lyapunov = NaN;
                            curve.results.s(i).data.matcontlyap = NaN;
                        end
                        %  curve.results.s(i).msg='Hopf';
                    case 'BTP'
                         curve.results.s(i).label='BT';
                    case 'FP'
                        if strcmp(curve.ctype,'F')
                            curve.results.s(i).label='CP';
                        end
                    case 'SN'
                        ind=curve.results.s(i).index;
                        curve.results.s(i).data.id='F';
                        if i>1&&curve.results.s(i-1).index==ind-1
                            if strcmp(curve.results.s(i-1).label,'BP')
                                curve.results.s(i).data.id='T';
                                curve.results.s(i).data.labnrBP=labnr{i-1};
                                %   curve.results.s(i).msg='Transcritical';
                            elseif strcmp(curve.results.s(i-1).label,'FP')
                                curve.results.s(i).data.id='F';
                                %     curve.results.s(i).msg='Saddle-Node';
                            end
                        end
                        if i<length(curve.results.s)&&curve.results.s(i+1).index==ind+1
                            if strcmp(curve.results.s(i+1).label,'BP')
                                curve.results.s(i).data.id='T';
                                curve.results.s(i).data.labnrBP=labnr{i+1};
                                %   curve.results.s(i).msg='Transcritical';
                            elseif strcmp(curve.results.s(i+1).label,'FP')
                                curve.results.s(i).data.id='F';
                                %     curve.results.s(i).msg='Saddle-Node';
                            end
                        end
                    case 'BP'
                        %   curve.results.s(i).msg='Branchpoint';
                end
                curve.results.s(i).data.labnr=labnr{i};
                curve.results.s(i).data.runid=therunid;
            end
        end
        function obj = updateprob(obj, refresh)
            global g_grind;
            
            if isempty(obj.settings.derived.prob)||refresh
                obj.settings.derived.prob = coco_prob();
                obj.settings.derived.prob = coco_set(obj.settings.derived.prob, 'ode', 'vectorized', strcmp(g_grind.solver.opt.Vectorized,'on'));
                obj.settings.derived.prob = coco_set(obj.settings.derived.prob, 'all', 'data_dir',''); %just use full file names as id
                obj.settings.derived.prob = coco_set(obj.settings.derived.prob,'cont','LogLevel',1); %show results in command window
                obj.settings.derived.prob = coco_set(obj.settings.derived.prob,'corr','LogLevel',1);%show results in command window
                obj.settings.derived.prob = coco_set(obj.settings.derived.prob,'ep', 'NSA',true,'BTP',true);%detect neutral saddle points and Bogdanov-Takens points
            end
            
            %ueer settings:
            settings1 = eval_settings(obj);
            tbs = fieldnames(settings1);
            for i = 1:length(tbs)
                tb = tbs{i};
                sett = transpose([fieldnames(settings1.(tb)), struct2cell(settings1.(tb))]);
                obj.settings.derived.prob = coco_set(obj.settings.derived.prob, tb, sett{:});
            end
            
            
            % s = sprintf('%s(:);', obj.settings.grind.activepars{:});
            % obj.settings.p0 = full(par('-vector'));
        end
        function  findingGH = addfuncsHB(obj)
            if obj.settings.grind.HH && (obj.settings.derived.ndim>=4)
                [data, uidx] = coco_get_func_data(obj.settings.derived.prob, 'ep', 'data', 'uidx');
                obj.settings.derived.prob = coco_add_func(obj.settings.derived.prob, 'coco_test_HH', @coco_test_HH, data, ...
                    'regular', 'test.HH', 'uidx', uidx);
                obj.settings.derived.prob = coco_add_event(obj.settings.derived.prob, 'HH', 'test.HH', 0);
            end
            
            if obj.settings.grind.FH&& (obj.settings.derived.ndim>=3)
                [data, uidx] = coco_get_func_data(obj.settings.derived.prob, 'ep', 'data', 'uidx');
                obj.settings.derived.prob = coco_add_func(obj.settings.derived.prob, 'coco_test_FH', @coco_test_FH, data, ...
                    'regular', 'test.FH', 'uidx', uidx);
                obj.settings.derived.prob = coco_add_event(obj.settings.derived.prob, 'FH', 'test.FH', 0);
            end
            
            if obj.settings.grind.GH
                findingGH=true;
                [data, uidx] = coco_get_func_data(obj.settings.derived.prob, 'ep', 'data', 'uidx');
                obj.settings.derived.prob = coco_add_func(obj.settings.derived.prob, 'lyap', i_getodehandle('coco_lyap',obj.settings.derived.allpars(obj.settings.grind.activepars)), data.ep_eqn, ...
                    'regular', 'test.L1', 'uidx', uidx);
                obj.settings.derived.prob = coco_add_event(obj.settings.derived.prob, 'GH', 'test.L1', 0);
            end
            
        end
        
        function expand_curve(obj,frompoint,ctype,backward)
            function fullcurve=append_curve(oldcurve,addedcurve)
                fullcurve=oldcurve;
                n=size(fullcurve.results.bd,1);
                for i=1:length(addedcurve.results.s)
                    addedcurve.results.s(i).index=addedcurve.results.s(i).index+n+1;
                end
                if length(addedcurve.results.s)>1
                    fullcurve.results.bd=[fullcurve.results.bd;addedcurve.results.bd(2:end,:)];
                    fullcurve.results.s=[fullcurve.results.s(1:end-1);addedcurve.results.s(2:end)];
                end
                fullcurve.data.runids=[fullcurve.data.runids;addedcurve.data.runids];
            end
            global g_grind;
            %              ndx=obj.getndx('curveprops','ctype',ctype);
            %             if any(strcmp(ctype,{'F','T'}))
            %                 ctype='SN';
            %             end
            if nargin<4
                expand_curve(obj,frompoint,ctype,true);
                expand_curve(obj,frompoint,ctype,false);
                return;
            end
            curvendx=obj.getndx('curves','frompoint',frompoint.id)&obj.getndx('curves','ctype',ctype);
            curve=obj.curves(curvendx);
            if ~length(curve)==1
                error('grind:coco:nosuchcurve','Curve (from %s to %s) does not exist',frompoint.id,ctype);
            end
            par1=obj.outfun('fun',obj.settings.derived.allpars{curve.freepars(1)},'curvendx',curvendx);
            parrange = obj.makerange(obj.settings.grind.parranges1,curve.freepars(1));
            if (min(par1)<=parrange(1)&&backward)||(max(par1)>=parrange(2)&&~backward)
                %do not expand if the curve is already filling the parrange
                return;
            end
            
            obj.updateprob(true);
            No = length(obj.curves) + 1;
            therunid=fullfile(grindpath, 'tmp',g_grind.odefile,int2str(No));
            
            contpar1 = obj.translatepars(obj.settings.derived.allpars(curve.freepars), obj.settings.derived.allpars);
            cocopars = obj.translatepars(obj.settings.derived.allpars(obj.settings.grind.activepars), obj.settings.derived.allpars);
            
            findingGH=false;
            if backward
                obj.settings.derived.prob  = coco_set(obj.settings.derived.prob, 'cont','PtMX',[0, obj.settings.cont.PtMX]);
                point0=curve.results.s(1);
            else
                obj.settings.derived.prob  = coco_set(obj.settings.derived.prob, 'cont','PtMX',[obj.settings.cont.PtMX, 0]);
                point0=curve.results.s(end);
                if isfield(curve.data,'expanded')
                   curve.data.expanded=curve.data.expanded+1;
                else
                   curve.data.expanded=1;
                end
            end
            % label=obj.pointprops(frompoint.propndx).label;
            switch ctype
                case 'EP'
                    %remove EP
                    obj.settings.derived.prob = ode_ep2ep(obj.settings.derived.prob,'', point0.data.runid, point0.data.labnr);
                    
                case 'H'
                    obj.settings.derived.prob = ode_HB2HB(obj.settings.derived.prob,'', point0.data.runid, point0.data.labnr);
                    findingGH = obj.addfuncsHB;
                case {'F','T'}
                    obj.settings.derived.prob = ode_SN2SN(obj.settings.derived.prob,'', point0.data.runid, point0.data.labnr);
                otherwise
                    return;
            end
            bd1 = coco(obj.settings.derived.prob,therunid, [], 1,contpar1, parrange); % cont toolbox arguments
            newcurve=struct('freepars',obj.settings.derived.freepars,'ctype',ctype,'frompoint',frompoint.id,...
                'data',[],'color',curve.color,'propndx',find(obj.getndx('curveprops','ctype',ctype)),'results',struct('stabil',[]));
            i_waitbar([]);
            newcurve.data=curve.data;
            newcurve.data.haslyapunov=findingGH;
            newcurve.data.runids={therunid};
            newcurve=obj.extract_curve_data(bd1,newcurve,cocopars);
            if backward
                obj.add_curve(append_curve(newcurve,curve));
            else
                obj.add_curve(append_curve(curve,newcurve));
            end
        end
        
        
        
        function [tspan,y0,options] = my_init(obj)
            global g_grind;
            han = obj.handles;
            y0=zeros(1,obj.settings.derived.ndim);
            options = odeset('Jacobian',han(3),'JacobianP',han(4),'Hessians',han(5),'HessiansP',han(6),'Vectorized',g_grind.solver.opt.Vectorized);
            tspan = [0 10];
        end
        
        function settings1 = eval_settings(obj)
            relh = obj.settings.grind.relh; %#ok<NASGU>
            if isempty(obj.settings.grind.parranges1)
                range = NaN;
            else
                range = obj.settings.grind.parranges1;
                range = abs(range(2) - range(1));
            end
            
            if isnan(range)&&~isempty(obj.settings.derived.freepars)
                range = abs(evalin('base', obj.settings.derived.allpars{obj.settings.derived.freepars(1)}));
            end
            
            if isnan(range)||range==0
                range = 1; %#ok<NASGU>
            end
            dim=numel(obj.settings.derived.freepars);
            defaultNAdapt = 0; %has to be adapted for PO
            %    if ndim==1&&strncmp(obj.settings.id0,'IO',2)||strncmp(obj.settings.id0,'H',1)
            %       defaultNAdapt = 1;
            %   end
            NCOL = obj.settings.corr.NCOL;
            if ischar(NCOL)
                NCOL = eval(NCOL);
            end
            NTST = obj.settings.corr.NTST;
            %          sett=obj.matcont_opt(:,[1 4]).';
            %           defsettings=struct(sett{:});
            %             if ~isempty(obj.settings.cont.Backward)
            %                 if obj.settings.matcont.Backward
            %                     obj.settings.derived.direction=obj.settings.derived.direction & [true false];
            %                 else
            %                     obj.settings.derived.direction=obj.settings.derived.direction & [false, true];
            %                 end
            %             end
            settings1=struct('cont',obj.settings.cont,'corr',obj.settings.corr);
            f=fieldnames(settings1.cont);
            for i=1:length(f)
                if ischar(settings1.cont.(f{i}))
                    if isempty(settings1.cont.(f{i}))
                        v=NaN;
                    else
                        v=eval(settings1.cont.(f{i}));
                    end
                    if ~isnan(v)
                        settings1.cont.(f{i})=v;
                    end
                end
            end
            f=fieldnames(settings1.corr);
            for i=1:length(f)
                if ischar(settings1.corr.(f{i}))
                    v=eval(settings1.corr.(f{i}));
                    if ~isnan(v)
                        settings1.corr.(f{i})=v;
                    end
                end
            end
            
        end
        
        
        
        function add_curve(obj,curve)
            if strcmp(curve.ctype,'EP')
                p0=obj.points(obj.getndx('points','id',curve.frompoint)).p0;
                curve=add_stability(obj,curve);
            end
            for i=1:length(curve.results.s)
                curve.results.s(i).label=strtrim(curve.results.s(i).label);
            end
            ids=obj.add_points(curve);
            for i=1:length(ids)
                curve.results.s(i).data.id=ids{i};
            end
            ndx=find(obj.getndx('curves','frompoint',curve.frompoint)&obj.getndx('curves','ctype',curve.ctype),1);
            if ~isempty(ndx)
                obj.curves(ndx)=curve;
            else
                obj.curves(end+1)=curve;
            end
        end
        
        
        function [point,pntndx]=create_point(obj,spoint,x0,p0)
            point = struct('id', '','p0',p0,'x0',x0,'data',spoint.data,'propndx',[]);
            if isfield(spoint,'id')
                point.id=spoint.id;
                point.propndx=find(obj.getndx('pointprops','label',regexp(spoint.id,'(R[1-4])|([A-Za-z+]*)','match','once')));
                lab='';
                pntndx=[];
            else
                ndx=find(obj.getndx('points','x0+p0',[x0;p0]));
                pntndx=[];
                if ~isempty(ndx)
                    labs=obj.pointprops([obj.points(ndx).propndx]).label;
                    pntndx=find(strcmp(labs,spoint.label),1);
                    if ~isempty(pntndx)
                        pntndx=ndx(pntndx);
                        point=obj.points(pntndx);
                        return;
                    end
                end
                point = struct('id', '','p0',p0,'x0',x0,'data',spoint.data,'propndx',[]);
                oldndx=find(obj.getndx('points','label',spoint.label));
                point.propndx=find(obj.getndx('pointprops','label',spoint.label));
                if length(point.propndx)==2&&isfield(spoint.data,'id')
                    point.propndx=find(obj.getndx('pointprops','ptype',spoint.data.id));
                end
                %find translation for label For instance LP->F
                if isempty(point.propndx)
                    point=[];
                    pntndx=-1;
                    return;
                end
                lab=obj.pointprops(point.propndx).ptype;
                
                point.id=[lab int2str(length(oldndx)+1)];
            end
            if isfield(point,'propndx')
                point.data.msg=obj.pointprops(point.propndx).descr;
            end
            switch lab %spoint.label
                 case {'EP','P'}
                    point.data.msg=spoint.msg;
                 case {'H','HB'}
                    if spoint.data.lyapunov>0
                        point.id=sprintf('H+%d',length(oldndx)+1);
                    else
                        point.id=sprintf('H%d',length(oldndx)+1);
                    end
            end
            
        end
        
        
        
        function init_engine(obj,remove)
            if nargin==1
                remove=false;
            end
            if xor(obj.settings.derived.opened,remove)
                return;
            end
            coco_in_path=exist('coco.m','file');
            if ~coco_in_path&&~isempty(obj.settings.derived.engine_path)||~exist(obj.settings.derived.engine_path,'dir')
                obj.settings.derived.engine_path=obj.findupdate;
            end
            cocodir=obj.settings.derived.engine_path;
            fullpath={fullfile(cocodir, 'core', 'toolbox'),fullfile(cocodir, 'covering', 'toolbox'),...
                fullfile(cocodir, 'ep', 'toolbox'),fullfile(cocodir, 'coll', 'toolbox'),...
                fullfile(cocodir, 'po', 'toolbox'),fullfile(cocodir, 'recipes'),...
                fullfile(cocodir, 'continex', 'toolbox')};
            
            if ~remove
                for i=1:length(fullpath)
                    addpath(fullpath{i});
                end
            elseif remove
                for i=1:length(fullpath)
                    rmpath(fullpath{i});
                end
            end
        end
    end
    
    properties (Hidden= true,Access=public)
        %Options:
%         %name,validation,description,matcont default
%         pairs= {'h0','h_min','h_max','h_fac_min(n>0)','h_fac_max(n>0)','PtMX(i>0)',...
%             'NPR(n)','NSV(n)','NAdapt','RMMX(n)','ItMX(i>0)','SubItMX(i>0)',...
%             'TOL(n>0)','ResTOL(n>0)','NTST(n)','NCOL(n>0)','NTSTMN(n)','NTSTMX(n)',...
%             'engine(e[coco|cocogrind|grind])','par(s)','par1(p)','par2(p)','parranges(n)','parranges1(n)','parranges2(n)','id0(s)','p0(n)','N0(n)',...
%             'relh(n>0)','allpars(c)','activepars(l)','symbolic(l)','silent(l)',...
%             'posonly(l)','N0ranges(n)','GH(l)','HH(l)','FH(l)','ngrid(i>0)'};
        
        cont_opt={'h_min','','Minimum step size','0.001*relh*range';...
            'h0','','Initial step size','0.01*relh*range';...
            'h_max','','Maximum step size','relh*range';...
            'h_fac_min','n>0','Minimum step size adaptation factor',0.5;...
            'h_fac_max','n>0','Maximum step size adaptation factor',2;...
            'PtMX','i>0','Maximum number of continuation steps',300;...
            'NPR','n','diagnostic output every NPR steps (default 10).',10;...
            'NSV','n','save solution every NSV steps, empty is the same as NPR (default []).','';...
            'NAdapt','','remesh frequency, (default 0(ep) or 1(po)).','defaultNAdapt';...
            'RMMX','n','maximum number of remesh loops (default 10)).',10};
        
        corr_opt={'ItMX','i>0','Maximum number of iterations before step is reduced',10;...
            'SubItMX','i>0','number of damping steps (default 7)',7;...
            'TOL','n>0','Tolerance of Newton correction',1e-006;...
            'ResTOL','n>0','converge criterion of norm of the residium	(default 1.00E-06).',1e-006;...
            'NTST','n','number of mesh intervals',50;...
            'NCOL','n>0','number of collocation nodes',5;...
            'NTSTMN','n','minimum number of mesh intervals','min(NTST,5)';...
            'NTSTMX','n','maximum number of mesh intervals','max(NTST,NCOL*dim)'};

        grind_opt={'relh','n>0','Step size relative to the range of the active parameter (default 0.1).',1;...
            'par1','p','first free parameter (default: '''')','';...
            'parranges1','n&length(n)==2','range for the first parameter (default: [NaN NaN])',[NaN NaN];...
            'par2','p#E','second free parameter (default: '''')','';...
            'parranges2','n&length(n)==2','range for the second parameter (default: [NaN NaN])',[NaN NaN];...
            'par3','p#E','third free parameter (default: '''')','';...
            'mindist','n>0','minimum distance between different special points (default: 1E-5)',1E-5;...
            'stateranges','n','min/max value for each of the state variables, or one row if all are the same (default: [Nan NaN])',[NaN NaN];...
            'symbolic','l','if available use the symbolic toolbox for Jacobians',true;...
            'activepars','l#c','indices of the active parameters (default: all 1)',true;...
            'GH','l','detect generalized Hopf or Bautin bifurcations (default true).',true;...
            'HH','l','detect Hopf-Hopf bifurcations (default true)',true;...
            'FH','l','detect fold-Hopf bifurcations (default true).',true};
        
        %ctype npars args handle descr
        curveprops=struct('ctype',{},'npars',[],'color',[],'handle',[],'descr',[])
        % id descr codim ctypes
        pointprops=struct('ptype',{},'label',[],'descr',[],'codim',[],'ctypes',[]);
        
        
    end
end

function curve=add_stability(obj,curve)
%recalculate eigenvalues to have all information
global g_grind;
s=obj.curveprops(curve.propndx).resfun(obj,curve);
x=permute(s.Y,[3,2,1]);
jac_handle=i_getodehandle('coco_jac');
eigenvals = zeros(size(x))+NaN;
if strcmp(g_grind.solver.opt.Vectorized, 'on')
    Js = jac_handle(transpose(x),transpose(s.parvalues));
    for i = 1:size(x, 1)-2
        eigenvals(i, :) = eig(Js(:, :, i));
    end
else
    for i = 1:size(x, 2)
        J = jac_handle(x(i, :).',s.parvalues(i, :).');
        eigenvals(i, :) = eig(J);
    end
end
stabil=zeros(size(x,1),1,'int16')-1;
iscomplex=any(imag(eigenvals)>0,2);
isstable=max(real(eigenvals),[],2)<=0;
issaddle=~isstable&min(real(eigenvals),[],2)<0;
stabil(isstable&~iscomplex)=1; %stable point
stabil(isstable&iscomplex)=2; %stable spiral
stabil(issaddle)=3;      %saddle point
stabil(~isstable&~issaddle&~iscomplex)=4; %unstable node
stabil(~isstable&~issaddle&iscomplex)=5; %unstable spiral
curve.results.stabil=stabil;
end
% --------------------------------------------------------------------------


function [s,xtra]=extract_coco_curve_eq(obj,curve)
%all parameters inculding the free pars
p0new=obj.points(obj.getndx('points','id',curve.frompoint)).p0.';
p0=p0new(ones(1,size(curve.results.bd,1)),:);
ndx=strcmp(curve.results.bd(1,:),'results.x');
Y=transpose(cat(2,curve.results.bd{2:end,ndx}));
Y(end+1,:)=NaN; %we end with nan to keep different curves separate
xtra.stabil=curve.results.stabil;
%xtra.stabil(end+1,1)=-1;
ndx=strcmp(curve.results.bd(1,:),'results.freepars');
fp=cat(1,curve.results.bd{2:end,ndx});
fp(end+1,:)=nan;
p0(:,curve.freepars)=fp;
%xtra.stabil=xtra.stabil;
s=struct('parvalues',p0,'pars',{obj.settings.derived.allpars},'Y',permute(Y, [3 2 1]),'t',zeros(1,1,size(Y,1)),'perm',[]);
end

function ok = isvalidrange(N0, N0ranges)
if isempty(N0ranges)
    ok = true;
else
    ndx = ~isnan(N0ranges(:, 1));
    if ~isempty(ndx)&&any(ndx)
        xlow =  min(N0(ndx) - N0ranges(ndx, 1));
    else
        xlow = inf;
    end
    
    ndx = ~isnan(N0ranges(:, 2));
    if ~isempty(ndx)&&any(ndx)
        xhigh =  min(N0ranges(ndx, 2) - N0(ndx));
    else
        xhigh = inf;
    end
    
    ok=min(xlow, xhigh) > -1e-5;
end

end

