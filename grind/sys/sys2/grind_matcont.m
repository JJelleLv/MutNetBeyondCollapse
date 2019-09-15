classdef grind_matcont < grind_cont
    
    
    methods (Access = public)
        
        function obj = grind_matcont
            global g_grind;
            %find path for matcont
            
            obj = obj@grind_cont;
            
            obj.settings.derived.ismap=g_grind.solver.isdiffer;
            
            obj.settings.derived.engine_path='';
            
            if obj.settings.derived.ismap
                obj.settings.derived.engine='matcontm';
            else
                obj.settings.derived.engine='matcont';
            end
            obj.open;
            obj.close;
            is=strfind(lower(obj.settings.derived.engine_path),obj.settings.derived.engine);
            if ~isempty(is)
                obj.settings.derived.version=obj.settings.derived.engine_path(is(end)+length(obj.settings.derived.engine):end);
            else
                obj.settings.derived.version='unknown';
            end
            if ~obj.settings.derived.ismap
                %label, npars, handle to curve,descr,iscycle,color,singularities
                %           '
                %
                
                crveprops={'1DMan',2,[],'1DMan - 1D manifold',false,'','b','';...
                    'BP',3,@branchpoint,'BP - Branchpoint of fold bifurcation',false,'bpds','b','';...
                    'BPC',3,@branchpointcycle,'BPC - Branchpoint of cycles',true,'lds','b','';...
                    'EP',1,@equilibrium,'EP - Equilibrium points',false,'eds','b','T,H,F';...
                    'H',2,@hopf,'H - Hopf bifurcation',false,'hds','r','BT,ZH,HH,GH';...
                    'HSN',2,@homoclinicsaddlenode,'HSN - Homoclinic saddle node',true,'','b','NCH';...
                    'HTHSN',2,@homotopysaddlenode,'HTHSN - Homotopy saddle node',true,'HTHSNds','b','';...
                    'HTHet',2,@homotopyHet,'HTHet - Homotopy heteroclinic',true,'HTHetds','b','';...
                    'HTHom',2,@homotopysaddle,'HTHom - Homotopy homoclinic',true,'HTHomds','b','';...
                    'Het',-1,@heteroclinic,'Het - Heteroclinic bifurcation',true,'hetds','b','';...
                    'Hom',2,@homoclinic,'Hom - Homoclinic bifurcation',true,'homds','b','NS,DRS,DRU,NDS,NDU,TLS,TLU,SH,NCH,BT,OFS,OFU,IFS,IFU';...
                    'LC',1,@limitcycle,'LC - Limit cycle',true,'lds',[0.4 0.4 0.4],'BPC,PD,LPC,NS';...
                    'LC2',2,@limitcycle,'LC2 - Limit cycle (codim 2)',true,'lds','b','';...
                    'F',2,@limitpoint,'F - Fold bifurcation',false,'lpds','g','BT,ZH,CP,BP';...
                    'T',2,@limitpoint,'T - Transcritical bifurcation',false,'lpds','b','BT,ZH,CP,BP';...
                    'LPC',2,@limitpointcycle,'LPC - Limit point cycle',true,'lds','b','R1_,CPC,LPNS,LPPD';...
                    'NS',2,@neimarksacker,'NS - Neimark sacker curve',true,'lds','b','R1_,R2_,R3_,R4_,LPNS,CH,PDNS,NSNS';...
                    'NS1',2,@neimarksacker,'NS1 - Neimark sacker curve',true,'lds','b','R1_,R2_,R3_,R4_,LPNS,CH,PDNS,NSNS';...
                    'NS2',2,@neimarksacker,'NS2 - Neimark sacker curve',true,'lds','b','R1_,R2_,R3_,R4_,LPNS,CH,PDNS,NSNS';...
                    'PD',2,@perioddoubling,'PD - Period doubling curve',true,'lds','b','R2_,LPPD,GPD,PDNS';...
                    'PD2',2,@perioddoubling,'PD2 - Period doubling x2 curve',true,'lds','b','R2_,LPPD,GPD,PDNS'};
            else
                %
                %
                crveprops={'FP2',1,@fixedpointmap,'FP2 - Fixed point 2 (maps)',false,'fpmds','b','NS,PD,F,T';...
                    'FP',1,@fixedpointmap,'FP - Fixed point (maps)',false,'fpmds','b','NS,PD,F,T';...
                    'Het',-1,@heteroclinic,'Het - Heteroclinic bifurcation',true,'hetds','b','F,T';...
                    'Hom',2,@homoclinic,'Hom - Homoclinic bifurcation',true,'homds','b','F,T';...
                    'T',2,@limitpointmap,'T - Transcritical bifurcation',false,'lpmds','b','R1_,LPPD,LPNS,CP';...
                    'F',2,@limitpointmap,'F - Fold bifurcation',false,'lpmds','g','R1_,LPPD,LPNS,CP';...
                    'LP2',2,@limitpointmap,'LP2 - Limit point map',false,'','b','R1_,LPPD,LPNS,CP';...
                    'LP',2,@limitpointmap,'LP Limit point map',false,'lpmds','b','R1_,LPPD,LPNS,CP';...
                    'LP4m1',2,@limitpointmap,'LP4m1 - Limit point map',false,'lpmds','b','R1_,LPPD,LPNS,CP';...
                    'LP4m2',2,@limitpointmap,'LP4m2 - Limit point map',false,'lpmds','b','R1_,LPPD,LPNS,CP';...
                    'NS2',2,@neimarksackermap,'NS2 - Neimark sacker curve',false,'nsmds','b','CH,PDNS,LPNS,R1_,NSNS,R2_,R3_,R4_';...
                    'NS3',2,@neimarksackermap,'NS3 - Neimark sacker curve',false,'nsmds','b','CH,PDNS,LPNS,R1_,NSNS,R2_,R3_,R4_';...
                    'NS4',2,@neimarksackermap,'NS4 - Neimark sacker curve',false,'nsmds','b','CH,PDNS,LPNS,R1_,NSNS,R2_,R3_,R4_';...
                    'NS',2,@neimarksackermap,'NS - Neimark sacker curve',false,'nsmds','b','CH,PDNS,LPNS,R1_,NSNS,R2_,R3_,R4_';...
                    'NS_Other',2,@neimarksackermap,'NS_Other',false,'nsmds','b','CH,PDNS,LPNS,R1_,NSNS,R2_,R3_,R4_';...
                    'NS_Same',2,@neimarksackermap,'NS_Same',false,'nsmds','b','CH,PDNS,LPNS,R1_,NSNS,R2_,R3_,R4_';...
                    'PD',2,@perioddoublingmap,'PD - Period doubling curve',false,'pdmds','b','R2_,LPPD,PDNS,GPD'};
            end
            obj.curveprops=struct('ctype',crveprops(:,1),'clabel',crveprops(:,1),'npars',crveprops(:,2),'iscycle',crveprops(:,5),...
                'handle',crveprops(:,3),'descr',crveprops(:,4),'curveds',crveprops(:,6),'color',crveprops(:,7),'resfun',@extract_curve_eq,'extrapars',[]);
            singul=regexp(crveprops(:,8),'[,]','split')';
            [obj.curveprops.singularity]=deal(singul{:});
            ndx=obj.getndx('curveprops','ctype','LC2');
            if any(ndx)
                obj.curveprops(ndx).extrapars={'[period]'};
                obj.curveprops(ndx).clabel='LC';
            end
            ndx=obj.getndx('curveprops','ctype','F');
            obj.curveprops(ndx).clabel='LP';
            ndx=obj.getndx('curveprops','ctype','T');
            obj.curveprops(ndx).clabel='LP';
            for i=1:length(obj.curveprops)
                if obj.curveprops(i).iscycle
                    obj.curveprops(i).resfun=@extract_curve_cycle;
                end
            end
            if ~obj.settings.derived.ismap
                %label, descr,   codims,   mindim, possible branches
                %(mindim between brackets = not yet checked
                pntprops={'BPC','Branch point of cycles',2,2,'BPC,LC,LPC';...
                    'BP','Branch point of fold curves',2,1,'BP';...
                    'BT','Bogdanov-Takens',2,2,'H,Hom,F';...
                    'CH','Chenciner (generalized NS)',2,3,'NS';...
                    'CP','Cusp',2,1,'F,EP';...
                    'CPC','Cusp of cycles',2,2,'LPC';...
                    'DRS','Double real stable leading eigenvalue',2,2,'';...
                    'DRU','Double real unstable leading eigenvalue',2,2,'';...
                    'EP','Equilibrium',0,1,'EP';...
                    'F','Fold',1,1,'EP,F';...
                    'GH','Generalized Hopf',2,2,'H,LPC';...
                    'GPD','Generalized Period Doubling',2,2,'PD';...
                    'H','Hopf',1,2,'H,EP,LC';...
                    'HH','Double Hopf',2,4,'H,NS1,NS2';...
                    'HHS','Homoclinic to Hyperbolic Saddle',1,(1),'';...
                    'HSN','Homoclinic to Saddle-Node',1,(1),'HSN,Hom';...
                    'IFS','Inclination-Flip with respect to the Stable manifold',2,(1),'';...
                    'IFU','Inclination-Flip with respect to the Unstable manifold',2,(1),'';...
                    'LC','Limit cycle',0,2,'HSN,Hom,LC2';...
                    'LPC','Limit Point of cycles',1,2,'LC,LPC';...
                    'LPNS','Fold-Neimark-Sacker',2,(1),'LPC,NS';...
                    'LPPD','Fold-period doubling',2,2,'PD';...
                    'NCH','Non-Central Homoclinic to saddle-node',2,(1),'HSN,Hom';...
                    'NDS','Neutrally-Divergent saddle-focus (Stable)',2,3,'';...
                    'NDU','Neutrally-Divergent saddle-focus (Unstable)',2,3,'';...
                    'NFF','Neutral Bi-Focus',2,(1),'';...
                    'NS','Neimark-Sacker (torus)',1,3,'LC,NS';...
                    'NSF','Neutral saddle-focus',2,(1),'';...
                    'NSNS','Double Neimark-Sacker',2,(1),'';...
                    'NSS','Neutral saddle',2,(1),'';...
                    'OFS','Orbit-Flip with respect to the Stable manifold',2,(1),'';...
                    'OFU','Orbit-Flip with respect to the Unstable manifold',2,(1),'';...
                    'P','Any initial conditions',0,1,'LC,EP';...
                    'PD','Period doubling (flip)',1,(1),'LC,LC2,PD';...
                    'PDNS','Flip-Neimark-Sacker',2,(1),'NS,PD';...
                    'R1','1:1 Resonance',2,(1),'LPC,NS';...
                    'R2','1:2 Resonance',2,(1),'NS,PD';...
                    'R3','1:3 Resonance',2,(1),'NS';...
                    'R4','1:4 Resonance',2,(1),'NS';...
                    'SH','Shilnikov-Hopf',2,(1),'';...
                    'T','Transcritical (branch point) bifurcation',1,1,'EP,T';...
                    'TLS','Three Leading eigenvalues (Stable)',2,3,'';...
                    'TLU','Three Leading eigenvalues (Unstable)',2,3,'';...
                    'ZH','Zero-Hopf',2,3,'H,F,NS';...
                    'Het','Heteroclinic bif.',1,2,'Het,Het';...
                    'Hom','Homoclinic bif.n',1,2,'HSN,Hom,Hom';...
                    'HTHet','Homotopy heteroclinic',1,(1),'HTHet,Het';...
                    'HTHom','Homotopy homoclinic',1,(1),'HTHom,Hom';...
                    'HTHSN','Homotopy saddle node',1,(1),'HSN,HTHSN'};
            else
                %label, descr,   codims,   mindim, possible branches
                %(mindim between brackets = not yet checked
                %MAPS: NS, PD, LP, BP,CH, PDNS, LPNS, R1-4, NSNS,LPPD, CP.
                %              'LP','Limit point (map)',1,(1),'LPm';...
                %                    'BP','Branch Point of Fold curves',2,(1),'BP';...
                
                pntprops={...
                    'CH','Chenciner (generalized Neimark-Sacker) bifurcation',2,2,'';...
                    'CP','Cusp bifurcation',2,1,'';...
                    'F','Fold bifurcation (limit point)',1,1,'F';...
                    'EP','Equilibrium',0,1,'FP';...
                    'FP','Fixed point (map)',0,1,'FP';...
                    'GPD','Generalized Period Doubling',2,(1),'LP2';...
                    'Het','Heteroclinic bifurcation',1,(1),'Het,Het';...
                    'HHS','Homoclinic to Hyperbolic Saddle',1,(1),'';...
                    'Hom','Homoclinic bifurcation',1,(1),'Hom';...
                    'IFS','Inclination-Flip with respect to the Stable manifold',2,(1),'';...
                    'IFU','Inclination-Flip with respect to the Unstable manifold',2,(1),'';...
                    'LPNS','Fold-Neimark-Sacker bifurcation',2,3,'';...
                    'LPPD','Fold-flip',2,2,'';...
                    'NDS','Neutrally-Divergent saddle-focus (Stable)',2,(1),'';...
                    'NDU','Neutrally-Divergent saddle-focus (Unstable)',2,(1),'';...
                    'NFF','Neutral Bi-Focus',2,(1),'';...
                    'NS','Neimark-Sacker(m)',1,2,'NS,NS_Other,NS_Same';...
                    'NSF','Neutral saddle-focus',2,(1),'';...
                    'NSNS','Double Neimark-Sacker',2,4,'';...
                    'NSS','Neutral saddle',2,(1),'';...
                    'P','Any initial conditions',0,1,'FP';...
                    'PD','Period Doubling (m)',1,1,'FP2,PD';...
                    'PDNS','Flip-Neimark-Sacker bifurcation',2,3,'NS2';...
                    'R1','1:1 Resonance',2,(1),'';...
                    'R2','1:2 Resonance',2,(1),'NS2';...
                    'R3','1:3 Resonance',2,(1),'NS3';...
                    'R4','1:4 Resonance',2,(1),'NS4,LP4m1,LP4m2';...
                    'SH','Shilnikov-Hopf',2,(1),'';...
                    'T','Transcritical (branch point) bifurcation',1,1,'FP,F'};
            end
            obj.pointprops=struct('ptype',regexprep(pntprops(:,1),'([0-9]$)','$1_'),'label',pntprops(:,1),'descr',pntprops(:,2),...
                'codim',pntprops(:,3),'mindim',pntprops(:,4),'ctypes',pntprops(:,5));
            for i=1:length(obj.pointprops)
                obj.pointprops(i).ctypes=regexp(obj.pointprops(i).ctypes,'[\,]','split');
            end
            ndx=obj.getndx('pointprops','ptype','TLS');
            if any(ndx)
                obj.pointprops(ndx).label='3LS';
            end
            ndx=obj.getndx('pointprops','ptype','TLU');
            if any(ndx)
                obj.pointprops(ndx).label='3LU';
            end
            ndx=obj.getndx('pointprops','ptype','F');
            obj.pointprops(ndx).label='LP';
            ndx=obj.getndx('pointprops','ptype','T');
            obj.pointprops(ndx).label='BP';
            if obj.settings.derived.ismap
                for i=1:length(obj.pointprops)
                    obj.pointprops(i).ctypes=regexprep(obj.pointprops(i).ctypes,'^EP$','FP');
                end
                %      ndx=obj.getndx('pointprops','ptype','EP');
                %      obj.pointprops(ndx).label={'FP'};
                options_not_in_map={'PRC','dPRC','Input'};
                obj.settings.matcont=rmfield(obj.settings.matcont,options_not_in_map);
                obj.settings.grind.niter=1;
                obj.settings.matcont.AutDerivative=true;
                obj.settings.matcont.AutDerivativeIte=24;
                obj.matcont_opt=obj.matcont_opt(~ismember(obj.matcont_opt(:,1),options_not_in_map),:);
                obj.grind_opt(end+1,:)={'niter','i>0','Number of iterations of the map (default=1)',1};
                obj.matcont_opt(end+1,:)= {'AutDerivative','l','Automatic differentiation of normal form coefficients(default=1)',true};
                obj.matcont_opt(end+1,:)= {'AutDerivativeIte','i','Automatic differentiation above this iter number (default=24)',24};
                %Test if adtayl is running, if not set AutDerivative to
                %false, the model runs slower
                han=i_getodehandle('normal');
                x0=i_initvar;
                try
                    obj.open;
                    x0=x0+adtayl(0,3).*ones(size(x0));
                    res=han(0,x0);
                    canrun=isa(res,'adtayl');
                catch err
                    disp('Warning cannot run AutoDerivative (using slower methods instead)');
                    disp(err.message);
                    canrun=false;
                end
                obj.close;
                obj.settings.matcont.AutDerivative=canrun;
                obj.matcont_opt{strcmp(obj.matcont_opt(:,1),'AutDerivative'),4}=canrun;
            end
            
            N0=i_initvar;
            s2 = sprintf('%g ', round(N0 * 10000) / 10000);
            obj.add_points(struct('id','P1','label','P','msg',sprintf('initial point: %s',s2),'data',[],'N0',N0));
            for i=1:numel(obj.curveprops)
                if numel(obj.curveprops(i).singularity)==1&&isempty(obj.curveprops(i).singularity{1})
                    obj.curveprops(i).singularity={};
                end
            end
            for i=1:numel(obj.pointprops)
                if numel(obj.pointprops(i).ctypes)==1&&isempty(obj.pointprops(i).ctypes{1})
                    obj.pointprops(i).ctypes={};
                end
            end
            obj.settings.derived.IgnoredSings={obj.pointprops([obj.pointprops.mindim]>obj.settings.derived.ndim).ptype};
        end
        
        
        function help(obj)
            open(obj);
            open('doc_matcont.html')
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
                obj.get('IgnoreSingularity');
                descr=[obj.matcont_opt;obj.grind_opt];
                allopts=[[fields(obj.settings.matcont);fields(obj.settings.grind)] [struct2cell(obj.settings.matcont);struct2cell(obj.settings.grind)]];
                %%% NOTE Here I assume that the options structure and grind_opt have the same order!!
                if nargin==2&&strcmp('-nondefault',varargin{1})
                    args=[descr(:,1),allopts(:,end)];
                    for i=size(args,1):-1:1
                        if strcmp(args{i,1},'IgnoreSingularity')
                            defaultsings={obj.pointprops([obj.pointprops.mindim]>obj.settings.derived.ndim).ptype};
                            setdiffer=setxor(obj.settings.derived.IgnoredSings,defaultsings);
                            if isempty(setdiffer)
                                args(i,:)=[];
                            else
                                ndx=find(ismember(setdiffer,defaultsings));
                                for i1=1:length(ndx)
                                    setdiffer{ndx(i1)}=['-' setdiffer{ndx(i1)}];
                                end
                                args{i,2}=setdiffer;
                            end
                        elseif struccmp(descr{i,4},args{i,end})
                            args(i,:)=[];
                        elseif strcmp(args{i,1},'activepars')&&(islogical(args{i,end})&&all(args{i,end}))
                            args(i,:)=[];
                        end
                    end
                    for i=1:size(args,1)
                        if iscell(args{i,2})&&~isempty(args{i,2})
                            args(i,2)={args(i,2)};
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
                        aval=obj.get(varargin{2});
                        if iscell(aval)
                            aval={aval};
                        end
                        args = struct('name',varargin{2},'type',descr{ndx,2},'descr',descr{ndx,3},'default',descr{ndx,4},'value',aval);
                        return;
                    end
                else
                    args=[descr(:,[1 3]),allopts(:,end)];
                    for i=1:size(args,1)
                        args{i,end}=any2str(args{i,end});
                    end
                end
                return;
            end
            if nargin==2
                opt=varargin{1};
                if isfield(obj.settings.matcont,opt)
                    if strcmpi(opt,'IgnoreSingularity')&&~isempty(obj.settings.derived.frompoint)
                        if isempty(obj.settings.derived.ctype)
                            ndx=obj.points(obj.getndx('points','id',obj.settings.derived.frompoint)).propndx;
                            obj.settings.derived.ctype=obj.pointprops(ndx).ctypes{1};
                        end
                        ndx=obj.getndx('curveprops','ctype',obj.settings.derived.ctype);
                        singlist=obj.curveprops(ndx).singularity;
                        args=find(ismember(singlist,obj.settings.derived.IgnoredSings));
                        obj.settings.matcont.IgnoreSingularity=args;
                    else
                        args=obj.settings.matcont.(opt);
                    end
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
                ndx=find(obj.getndx('curveprops','ctype',varargin{3}),1);
                codim=obj.curveprops(ndx).npars;
                if ~isempty(varargin{2})
                    ndx2=find(obj.getndx('points','id',varargin{2}),1);
                    if ~any(ndx2)
                        ndx2=obj.getndx('pointprops','ptype',varargin{2});
                    else
                        ndx2=obj.points(ndx2).propndx;
                    end
                    if ~any(ndx2)
                        error('grind:matcont','Point "%s" not found',varargin{2})
                    end
                    [~,funargs]=obj.get_init_fun(obj.pointprops(ndx2).ptype,obj.curveprops(ndx).ctype);
                    funargs=regexp(funargs,'[^(^)]*','match','once');
                else
                    funargs={};
                end
                if isempty(funargs)
                    found=false;
                end
                %  ndx=find(obj.getndx('curveprops','ctype',varargin{3}),1);
                if obj.curveprops(ndx).npars==3
                    funargs=[funargs {'par3'}];
                end
                f=[fieldnames(obj.settings.matcont);fieldnames(obj.settings.grind)];
                ndx=ismember(funargs,f);
                %   keyopts=[funargs(ndx) {'relh'}];
                keyopts=[funargs(ndx) {'StepAmp','MaxNumPoints'}];    
                args=obj.get;
                if strcmp(varargin{1},'alloptions')
                    ndx1=true(size(args,1),1);
                    f=find(strcmp(obj.first_argument_par,args(:,1)));
                    ndx1(f:end)=false;
                    if codim==1
                        ndx1(ismember(args(:,1),{'par2','par3','parranges2'}))=false;
                    elseif codim==2
                        ndx1(ismember(args(:,1),{'par3'}))=false;
                    end
                    ndx1(ismember(args(:,1),keyopts))=true;
                    ndx=ndx1;
                else
                    ndx=ismember(args(:,1),keyopts);
                end
                args=args(ndx,:);
            end
        end
        function args1=set(obj,varargin)
            %% set options
            if nargin==1
                %without arguments a list of current settings is given
                disp('Current settings:');
                fprintf('g_cont.set(...\n');
                obj.get('IgnoreSingularity');
                fields1=fieldnames(obj.settings.grind);
                fields2=fieldnames(obj.settings.matcont);
                f=[fields1; fields2];
                f1=[repmat({'grind'},size(fields1)); repmat({'matcont'},size(fields2))];
                [~,ndx]=sort(lower(f));
                f=f(ndx);
                f1=f1(ndx);
                for j=1:length(f)
                    fprintf('    ''%s'',%s',f{j},any2str(obj.settings.(f1{j}).(f{j})));
                    if j==length(f1)
                        fprintf(');\n');
                    else
                        fprintf(',...\n')
                    end
                end
                return;
            end
%             pairs1=obj.matcont_opt(:,1);
%             for i=1:length(pairs1)
%                 if ~isempty(obj.matcont_opt(i,2))
%                     pairs1{i}=sprintf('%s(%s)',pairs1{i},obj.matcont_opt{i,2});
%                 end
%             end
%             pairs2=obj.grind_opt(:,1);
%             for i=1:length(pairs2)
%                 if ~isempty(obj.grind_opt(i,2))
%                     pairs2{i}=sprintf('%s(%s)',pairs2{i},obj.grind_opt{i,2});
%                 end
%             end

            ext_opts={'frompoint','c#s','code of the point from which the continuation starts (for instance ''EP1'') (default: empty)','';...
                        'ctype','c#s','type of curve for continuations (for instance EP)','';...
                         'file','s','open a previously saved conteq session',''};

            %pairs3={'frompoint(c#s)';'ctype(c#s)';'file(s)'};
            args= i_parseargs([obj.matcont_opt;obj.grind_opt;ext_opts].','if nargs==1,deffields=''file'';else,deffields=''par1,frompoint,ctype'';end;','-c,-out,-l,-defaults,-coco,-matcont',varargin);
            if any(strcmp(args.opts,'-defaults'))
                obj.setdefaults;
            end
            if isfield(args,'frompoint')  %if we only define frompoint ctype should be emptied
                obj.settings.derived.ctype='';
            end
            if isfield(args,'IgnoreSingularity')
                if isfield(args,'ctype')
                    obj.settings.derived.ctype=args.ctype;
                end
                ndx=obj.getndx('curveprops','ctype',obj.settings.derived.ctype);
                singlist=obj.curveprops(ndx).singularity;
                
                if ischar(args.IgnoreSingularity)
                    if isempty(args.IgnoreSingularity)
                        args.IgnoreSingularity={};
                    else
                        args.IgnoreSingularity=regexp(args.IgnoreSingularity,',','split');
                    end
                end
                if iscell(args.IgnoreSingularity)
                    %To remove a singularity from the default list add "-"
                    %in front of the label
                    %
                    defaultsings={obj.pointprops([obj.pointprops.mindim]>obj.settings.derived.ndim).ptype};
                    removelist=regexp(args.IgnoreSingularity,'(?<=^-).*','match','once');
                    ndx=cellfun('isempty',removelist);
                    removelist=removelist(~ndx);
                    args.IgnoreSingularity(~ndx)=removelist;
                    for i1=1:length(args.IgnoreSingularity)
                        if ~any(obj.getndx('pointprops','ptype',args.IgnoreSingularity{i1}))
                            error('grind:matcont:set','Unknown singularity: %s',args.IgnoreSingularity{i1});
                        end
                    end
                    i1=ismember(defaultsings,removelist);
                    defaultsings(i1)=[];
                    obj.settings.derived.IgnoredSings=unique([defaultsings args.IgnoreSingularity(ndx)]);
                elseif isnumeric(args.IgnoreSingularity)
                    if max(args.IgnoreSingularity)>length(singlist)
                        error('grind:matcont','Singularity index too large')
                    end
                    ndx=false(size(singlist));
                    ndx(args.IgnoreSingularity)=true;
                    ndx1= ismember(obj.settings.derived.IgnoredSings,singlist(~ndx));
                    obj.settings.derived.IgnoredSings(ndx1)=[];
                    obj.settings.derived.IgnoredSings=unique([obj.settings.derived.IgnoredSings,singlist(ndx)]);
                end
            end
            if isfield(args,'UserfunctionsInfo')&&isfield(args.UserfunctionsInfo,'label')
                ndx=[obj.pointprops(:).mindim]~=-1;
                user_labels={obj.pointprops(~ndx).label};
                all_labels={obj.pointprops(ndx).label,obj.pointprops(ndx).ptype};
                changed=false;
                for i=1:length(args.UserfunctionsInfo)
                    while any(args.UserfunctionsInfo(i).label(end)=='0123456789')||any(strcmp(args.UserfunctionsInfo(i).label,all_labels))
                        args.UserfunctionsInfo(i).label=[args.UserfunctionsInfo(i).label '_'];
                        changed=true;
                    end
                    if any(strcmp(args.UserfunctionsInfo(i).label,user_labels))
                        ndx1=obj.getndx('pointprops','label',args.UserfunctionsInfo(i).label);
                        obj.pointprops(ndx1).descr=args.UserfunctionsInfo(i).name;
                    end
                end
                if changed
                    warning('grind:set:matcont','Changed labels of userfunction(s) to make them unique: %s',sprintf('%s ',args.UserfunctionsInfo(:).label))
                end
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
            matcontfields=obj.matcont_opt(:,1);
            derivedfields={'frompoint';'ctype'};
            for i=1:length(f)
                if any(strcmp(f{i},grindfields))
                    obj.settings.grind.(f{i})=args.(f{i});
                elseif any(strcmp(f{i},derivedfields))
                    obj.settings.derived.(f{i})=args.(f{i});
                elseif any(strcmp(f{i},matcontfields))
                    obj.settings.matcont.(f{i})=args.(f{i});
                end
            end
            if ~isempty(obj.settings.grind.par2)&&strcmp(obj.settings.grind.par1,obj.settings.grind.par2)
                obj.settings.grind.par2='';
                error('grind:matcont','The second free parameter must be different from the first')
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
                warning('grindmatcont:notactive','parameter %s is not active',obj.settings.grind.par1)
                ipar1=find(obj.settings.grind.activepars,1);
                args.par1=obj.settings.derived.allpars{ipar1};
                obj.settings.grind.par1=args.par1;
            end
            ipar2=find(strcmp(obj.settings.grind.par2, obj.settings.derived.allpars),1);
            if ~isempty(ipar2)&&~obj.settings.grind.activepars(ipar2)
                warning('grindmatcont:notactive','parameter %s is not active',obj.settings.grind.par2)
                pars=obj.settings.grind.activepars;
                pars(ipar1)=false;
                ipar2=find(pars,1);
                args.par2=obj.settings.derived.allpars{ipar2};
                obj.settings.grind.par2=args.par2;
            elseif isempty(ipar2)&&~isempty(obj.settings.grind.par2)
                error('grind:matcont:nopar','%s is not a parameter',args.par2)
            end
            ipar3=find(strcmp(obj.settings.grind.par3, obj.settings.derived.allpars),1);
            
            obj.settings.derived.freepars=[ipar1,ipar2,ipar3];
            
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
            if nargout>0
                args1=args;
            end
            
        end
        
        
        function out=handles(obj,filename)
            global  g_grind;
            if nargin<2
                filename='';
            end
            if ~isempty(filename)
                fid=fopen(filename,'w');
                [~,name]=fileparts(filename);
                fprintf(fid,'function out = %s\n',name);
                fprintf(fid,'out{9}=[];\nout{1} = @init;\nout{2} = @fun_eval;\n');
                if isempty(g_grind.syms.errormsg)&&obj.settings.grind.symbolic
                    if ~isempty(g_grind.syms.Jacobian)
                        fprintf(fid,'out{3} = @jacobian;\n');
                    end
                    if ~isempty(g_grind.syms.Jacobianp)
                        fprintf(fid,'out{4} = @jacobianp;\n');
                    end
                    if ~isempty(g_grind.syms.Hessian)
                        fprintf(fid,'out{5} = @hessian;\n');
                    end
                    if ~isempty(g_grind.syms.Hessianp)
                        fprintf(fid,'out{6} = @hessianp;\n');
                    end
                    if ~isempty(g_grind.syms.der3)
                        fprintf(fid,'out{7} = @der3;\n');
                    end
                    if ~isempty(g_grind.syms.der4)
                        fprintf(fid,'out{8} = @der4;\n');
                    end
                    if ~isempty(g_grind.syms.der5)
                        fprintf(fid,'out{9} = @der5;\n');
                    end
                    %                         out{3} = i_getodehandle('Jacobian',activepars);
                    %                         out{4} = i_getodehandle('Jacobianp',activepars);
                    %                         out{5} = i_getodehandle('Hessian',activepars);
                    %                         out{6} = i_getodehandle('Hessianp',activepars);
                    %                         out{7} = i_getodehandle('der3',activepars);
                    %                         out{8} = i_getodehandle('der4',activepars);
                    %                         out{9} = i_getodehandle('der5',activepars);
                end
                separ='% --------------------------------------------------------------------------';
                fprintf(fid,'\n\n%s\n',separ);
                fprintf(fid,'function [tspan,y0,options] = init\nhandles = %s;\ny0=zeros(%d,1);\n',name,g_grind.statevars.dim);
                fprintf(fid,'options = odeset(''Jacobian'',handles(3),''JacobianP'',handles(4),''Hessians'',handles(5),''HessiansP'',handles(6));\ntspan = [0 10];\n');
                fprintf(fid,'\n\n%s\n',separ);
                fid2=fopen([g_grind.odefile '.m'],'r');
                odefile=textscan(fid2,'%s', 'delimiter','\n', 'whitespace','');
                odefile=odefile{1};
                fclose(fid2);
                funline=regexprep(odefile{2},sprintf('=%s(',g_grind.odefile),'=%s(');
                odefile{2}=sprintf(funline,'fun_eval');
                
                fprintf(fid,'%s\n',odefile{:});
                if isempty(g_grind.syms.errormsg)
                    if ~isempty(g_grind.syms.Jacobian)
                        plotjac(fid,g_grind.syms.Jacobian,separ,sprintf(funline,'jacobian'))
                    end
                    if ~isempty(g_grind.syms.Jacobianp)
                        plotjac(fid,g_grind.syms.Jacobianp,separ,sprintf(funline,'jacobianp'))
                    end
                    if ~isempty(g_grind.syms.Hessian)
                        plotjac(fid,g_grind.syms.Hessian,separ,sprintf(funline,'hessian'))
                    end
                    if ~isempty(g_grind.syms.Hessianp)
                        plotjac(fid,g_grind.syms.Hessianp,separ,sprintf(funline,'hessianp'))
                    end
                    if ~isempty(g_grind.syms.der3)
                        plotjac(fid,g_grind.syms.der3,separ,sprintf(funline,'der3'))
                    end
                    if ~isempty(g_grind.syms.der4)
                        plotjac(fid,g_grind.syms.der4,separ,sprintf(funline,'der4'))
                    end
                    if ~isempty(g_grind.syms.der5)
                        plotjac(fid,g_grind.syms.der5,separ,sprintf(funline,'der5'))
                    end
                end
                
                fclose(fid);
            else
                %gets the matcont handles to the grind model
                out{9}=[];
                out{1} = @()my_init(obj);
                out{2} = i_getodehandle('matcont',obj.settings.derived.allpars(obj.settings.grind.activepars));
                activepars=sprintf(',%s',obj.settings.derived.allpars{obj.settings.grind.activepars});
                %matcont handles vectorized model correctly
                if isempty(g_grind.syms.errormsg)&&obj.settings.grind.symbolic
                    if ~isempty(g_grind.syms.Jacobian)
                        out{3} = i_getodehandle('Jacobian',activepars);
                    end
                    if ~isempty(g_grind.syms.Jacobianp)
                        out{4} = i_getodehandle('Jacobianp',activepars);
                    end
                    if ~isempty(g_grind.syms.Hessian)
                        out{5} = i_getodehandle('Hessian',activepars);
                    end
                    if ~isempty(g_grind.syms.Hessianp)
                        out{6} = i_getodehandle('Hessianp',activepars);
                    end
                    if ~isempty(g_grind.syms.der3)
                        out{7} = i_getodehandle('der3',activepars);
                    end
                    if ~isempty(g_grind.syms.der4)
                        out{8} = i_getodehandle('der4',activepars);
                    end
                    if ~isempty(g_grind.syms.der5)
                        out{9} = i_getodehandle('der5',activepars);
                    end
                end
                if ~isempty(obj.settings.grind.Userhandles)
                    if isa(obj.settings.grind.Userhandles,'function_handle')
                        out{10}=obj.settings.grind.Userhandles;
                    else
                        out=[out obj.settings.grind.Userhandles(:)'];
                    end
                end
            end
        end
        
        function internal_consistency(obj)
            global lpds lds HTHetds HTHSNds HTHomds;
            obj.open;
            nproblems=0;
            lpds.nphase=obj.settings.derived.ndim;
            lpds.BranchParams=[];
            lds.BranchParams=[];
            HTHetds.index=1;
            HTHSNds.index=1;
            HTHomds.index = 1;
            for i=1:length(obj.curveprops)
                hh=obj.curveprops(i).handle;
                if ~isempty(hh)
                    try
                        bb=hh();
                        [~,slist]=feval(bb{9});
                        slist=str2cell(slist)';
                    catch err
                        slist={};
                        
                        if strcmp('MATLAB:TooManyOutputs',err.identifier)
                            % fprintf('%s,'''';...\n',obj.curveprops(i).ctype)
                            
                        elseif strcmp('MATLAB:UndefinedFunction',err.identifier)
                            fprintf('%s does not exist\n',func2str(hh));
                            nproblems=nproblems+1;
                        else
                            fprintf('cannot run %s\n',func2str(hh));
                            nproblems=nproblems+1;
                        end
                    end
                    slist1=slist;
                    k=0;
                    for j=1:length(slist1)
                        ndx=find(obj.getndx('pointprops','label',slist1{j}));
                        if isempty(ndx)
                            fprintf('Cannot find %s\n',slist1{j});
                            nproblems=nproblems+1;
                        else
                            slist(j+k:j+k+length(ndx)-1)={obj.pointprops(ndx).ptype};
                            k=k+length(ndx)-1;
                        end
                    end
                    curve_slist=obj.curveprops(i).singularity;
                    if ~all(ismember(curve_slist,slist))
                        disp(slist);
                        disp(obj.curveprops(i).singularity)
                        nproblems=nproblems+1;
                        fprintf('list of singularities curve "%s" does not match\n',obj.curveprops(i).ctype);
                    end
                    %                    s=sprintf('%s,',slist{:});
                    %                    fprintf('%s,''%s'';...\n',obj.curveprops(i).ctype,s(1:end-1))
                end
            end
            for i=1:numel(obj.pointprops)
                pnt=obj.pointprops(i).ptype;
                for j=1:numel(obj.pointprops(i).ctypes)
                    [initfun,argnames]=obj.get_init_fun(pnt,obj.pointprops(i).ctypes{j});
                    if isempty(initfun)
                        nproblems=nproblems+1;
                        fprintf('cannot get name for init function init_%s_%s\n',pnt,obj.pointprops(i).ctypes{j});
                    end
                    if ~isempty(initfun)&&~exist(func2str(initfun),'file')
                        nproblems=nproblems+1;
                        fprintf('cannot find init function %s\n',func2str(initfun));
                    elseif ~isempty(initfun)
                        fid=fopen([func2str(initfun) '.m'],'r');
                        lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                        fclose(fid);
                        lines=lines{1};
                        f=find(strncmp(lines,'function ',9),1);
                        fun=lines{f};
                        f=strfind(fun,'(');
                        fun=fun(f:end);
                        read_args=regexp(fun,'(?<=[,\( ])[A-Za-z0-9_]*','match');
                        ndx=strncmp(argnames,'dumm:',5);
                        argnames=argnames(~ndx);
                        if numel(read_args)~=numel(argnames)
                            nproblems=nproblems+1;
                            fprintf('number of arguments differs %s\n',func2str(initfun));
                        else
                            syns={'x','xtot';'bp','BranchingPar';'bp','BranchingPars';'amplitude','amp';'eps','amp';'h','amp';...
                                'up','UParam';'sp','SParam';'asp','ActiveSParam';'aup','ActiveUParam';'varargin','(BranchingPars)';...
                                'J','niter';'n','niter'};
                            for k=1:numel(read_args)
                                if ~strcmpi(read_args{k},argnames{k})
                                    ndx1=strcmp(read_args{k},syns(:,1));
                                    if ~any(ndx1)||~any(strcmp(argnames{k},syns(ndx1,2)))
                                        nproblems=nproblems+1;
                                        fprintf('Arg %d differs %s : %s <> %s\n',i,func2str(initfun),read_args{k},argnames{k});
                                    end
                                end
                            end
                            
                        end
                    end
                end
            end
            [~,~,init_funcs]=obj.get_init_fun(pnt,'');
            for i=1:size(init_funcs,1)
                ids=regexp(init_funcs{i,1},'[_](?![_])','split');
                ndx1=obj.getndx('pointprops','ptype',ids{1});
                ndx2=obj.getndx('curveprops','ctype',ids{2});
                if ~any(ndx1)
                    fprintf('Init_functs %s: Point %s not found\n',init_funcs{i,1},ids{1})
                    nproblems=nproblems+1;
                end
                if ~any(ndx2)
                    fprintf('Init_functs %s: Curve %s not found\n',init_funcs{i,1},ids{2})
                    nproblems=nproblems+1;
                end
            end
            if nproblems==0
                disp('No inconsistencies found');
            else
                fprintf('%d inconsistencies\n',nproblems);
            end
            obj.close;
        end
    end
    
    methods (Static, Access = public)
        function res=multipliers(curve)
            global g_grind;
            if strcmp(curve.ctype,'LC')
                %all cyclic solutions?
                if isfield(curve.data.curveds,'nphase')
                    dim=curve.data.curveds.nphase;
                end
                dim2=size(curve.results.f,2);
                res=curve.results.f(end-dim+1:end,:);
                %remove one multiplier that is closest to 1
                %(as the multipliers of Poincaremap have dim-1)
                [min1,ndx]=min(abs(res-1),[],1);
                ndx1=sub2ind([dim,dim2],ndx,1:dim2);
                res(ndx1)=[];
                res=reshape(res,[dim-1,dim2]);
                res(:,min1>1e-4)=nan;
            elseif any(strcmp(curve.ctype,{'FP','FP2'}))
                dim=g_grind.statevars.dim;
                res=curve.results.f(end-dim+1:end,:);
            else
                res=[];
            end
        end
        
        function repairBugs(ismap,engine_path) %I found a bug in MATCONTM that is crucial
            %check new versions if it is repaired
            if nargin==0
                grind_matcont.repairBugs(true);
                grind_matcont.repairBugs(false);
                return;
            end
            if nargin<2
                engine_path=grind_matcont.findupdate(false, '',ismap);
            end
            oldcd=cd;
            if ismap
                %version MatContM_2017.09.06,matcontm5p3,matcontm5p4
                cd(engine_path);
                cd AD
                cd @adtayl
                fid=fopen('subsasgn.m','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                if ~any(strcontains(lines,'isempty(a)'))
                    lines{11}= sprintf('      if isempty(a)%%grind\n         a=b;\n      end\n%s',lines{11});
                    fid=fopen('subsasgn.m','w');
                    fprintf(fid,'%s\n',lines{:});
                    fclose(fid);
                end
                
                %version MatContM_2017.09.06,matcontm5p3,matcontm5p4
                %
                cd(engine_path);
                cd Continuer
                fid=fopen('contmrg.m','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                if any(strcontains(lines,'eval(['))
                    lines{18}='  val=opt.(strtrim(allopt(i,:)));%grind' ;
                    lines{20}='     options.(strtrim(allopt(i,:)))=val;%grind';
                    fid=fopen('contmrg.m','w');
                    fprintf(fid,'%s\n',lines{:});
                    fclose(fid);
                end
                %version MatContM_2017.09.06, repaired in matcontm5p3
                %
                cd(engine_path);
                cd LimitPointMap
                fid=fopen('lpvecthessvect.m','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                ndx=strcontains(lines,'jac(:,:,h)=lpmjac(xit,p,h);');
                if any(ndx)
                    lines(ndx)={'  jac(:,:,h)=lpmjac(x1,p,h);%grind'};
                    fid=fopen('lpvecthessvect.m','w');
                    fprintf(fid,'%s\n',lines{:});
                    fclose(fid);
                end
                
                %                   val=opt.(strtrim(allopt(i,:)));%  eval(['val = opt.' allopt(i,:) ';']);
                %   if ~isempty(val)
                %       options.(strtrim(allopt(i,:)))=val; %   eval(['options.' allopt(i,:) '= val;']);
                %   end
                
            else
                %version matcont6p6,matcont6p7,matcont6p8, matcont6p9,
                %matcont6p10, matcont6p11
                %init_BT_Hom
                %replace disp(table(a,b,d,e,a1,b1));
                %with fprintf('a = %g, b = %g, d = %g,\ne = %g, a1 = %g, b1 = %g\n',a,b,d,e,a1,b1);
                cd(engine_path);
                cd Homoclinic
                fid=fopen('init_BT_Hom.m','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                ndx=strcontains(lines,'disp(table(a,b,d,e,a1,b1));');
                if any(ndx)
                    lines(ndx)={'fprintf(''a = %g, b = %g, d = %g,\ne = %g, a1 = %g, b1 = %g\n'',a,b,d,e,a1,b1); %grind'};
                    fid=fopen('init_BT_Hom.m','w');
                    fprintf(fid,'%s\n',lines{:});
                    fclose(fid);
                end
                
                %version matcont6p6,matcont6p7,matcont6p8 repaired in matcont6p9
                %nf_NSNS missing rearr functions
                rearrfunc=sprintf('function [x,p,T] = rearr(x0) %%grind\n%% [x,p] = rearr(x0)\n%% Rearranges x0 into coordinates (x) and parameters (p)\nglobal lds\n\np = lds.P0;\np(lds.ActiveParams) = x0(lds.PeriodIdx+1:lds.PeriodIdx+2);\nx = x0(lds.coords);\nT = x0(lds.PeriodIdx);\n');
                
                cd(engine_path);
                cd LimitCycleCodim2
                fid=fopen('nf_NSNS.m','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                ndx=strcontains(lines,'function [x,p,T] = rearr(x0)');
                if ~any(ndx)
                    lines(end+1)={rearrfunc};
                    fid=fopen('nf_NSNS.m','w');
                    fprintf(fid,'%s\n',lines{:});
                    fclose(fid);
                end
            end
            cd(oldcd);
        end
        
        function [engine_path,matcontfile]=findupdate(checkupdates, url,ismap)
            %try to find the path of matcont and writes to grind.cfg for
            %future use. If matcont is not found it is downloaded from
            %url.
            %Checkupdates checks if there are newer versions.
            function writecfg(engine_path)
                if exist('grind.cfg','file')
                    fid=fopen('grind.cfg','r');
                    lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                    lines=lines{1};
                    fclose(fid);
                    ndx=~cellfun('isempty',(strfind(lines,fileitem)));
                    if ~any(ndx)
                        ndx=length(lines)+1;
                    end
                    lines{ndx}=sprintf('%s%s',fileitem,engine_path);
                else
                    lines={sprintf('%s%s',fileitem,engine_path)};
                end
                fid=fopen(fullfile(grindpath,'grind.cfg'),'w');
                fprintf(fid,'%s\n',lines{:});
                fclose(fid);
            end
            global g_grind
            if nargin==0
                checkupdates=[]; %only update if engine is not found
            end
            engine_path='';
            if nargin<=1||isempty(url)
                url = 'http://content.alterra.wur.nl/webdocs/internet/aew/downloads/';
            end
            if nargin<=2||isempty(ismap)
                if isempty(g_grind)
                    ismap=false;
                else
                    ismap=g_grind.solver.isdiffer;
                end
            end
            
            if ismap
                matcontfile='matcontm.m';
                fileitem='matcontmdir=';
                statfile='matcontmversion.txt';
            else
                matcontfile='matcont.m';
                fileitem='matcontdir=';
                statfile='matcontversion.txt';
            end
            % if isempty(engine_path)||~exist(engine_path,'dir')
            if exist('grind.cfg','file')
                fid=fopen('grind.cfg','r');
                lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
                lines=lines{1};
                fclose(fid);
                ndx=~cellfun('isempty',(strfind(lines,fileitem)));
                if any(ndx)
                    line=lines{ndx};
                    engine_path=regexp(line,'(?<=[=]).*','match', 'once');
                    if ~exist(engine_path,'dir')||~exist(fullfile(engine_path,matcontfile),'file')
                        engine_path='';
                    end
                end
            else
                engine_path =fileparts(which(matcontfile));
            end
            % end
            if isempty(engine_path)||~exist(engine_path,'dir')
                engine_path = findgrindfile(matcontfile);
                if isempty(engine_path)
                    if isempty(checkupdates)
                        checkupdates=true;
                    end
                else
                    writecfg(engine_path)
                end
            end
            
            if ~isempty(checkupdates)&& checkupdates
                [newversion, status] = urlread([url statfile]);
                if newversion(end) + 0 == 10
                    newversion = newversion(1:end - 1);
                end
                
                if status == 1
                    if ~isempty(engine_path)
                        if strcontains(engine_path,newversion)
                            if ismap
                                disp('No newer version of MATCONTM available');
                            else
                                disp('No newer version of MATCONT available');
                            end
                            return;
                        end
                    end
                    
                    engine_path = grindpath;
                    f = strfind(engine_path, 'grind');
                    if ~isempty(f)
                        engine_path = engine_path(1:f(1) - 2);
                    end
                    
                    zipfile=[newversion '.zip'];
                    fprintf('Downloading %s\n', zipfile)
                    unzip([url zipfile], engine_path);
                    engine_path = fullfile(engine_path, newversion);
                    writecfg(engine_path)
                    grind_matcont.repairBugs(ismap,engine_path); %I found a bug in MATCONTM that is crucial
                    if ismap
                        disp('Updated MATCONTM');
                    else
                        disp('Updated MATCONT');
                    end
                end
                
            end
        end
        function  [init_fun,argnames,init_functs]=get_init_fun(ptype,clabel)
            global g_grind;
            if ~g_grind.solver.isdiffer
                
                
                init_functs={'BP_BP','BP_BP',{'odefile','x','p','ap','BranchingPar'};...
                    'BP_F','BP_LP',{'odefile','x','p','ap'};...
                    'BPC_BPC','BPC_BPC',{'odefile','xtot','s','ap','NTST','NCOL','BranchingPars'};...
                    'BPC_LC','BPC_LC',{'odefile','xtot','v','s','NTST','NCOL','amp'};...
                    'BPC_LPC','BPC_LPC',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'BT_F','BT_LP',{'odefile','x','p','ap'};...
                    'BT_H','BT_H',{'odefile','x','p','ap'};...
                    'BT_Hom','BT_Hom',{'odefile','x','s','p','ap','NTST','NCOL','TTolerance','amp','extravec'};...
                    'CH_NS','CHNS_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'CP_EP','EP_EP',{'odefile','x','p','ap'};...
                    'CP_F','CP_LP',{'odefile','x','p','ap'};...
                    'CPC_LPC','CPC_LPC',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'EP_EP','EP_EP',{'odefile','x','p','ap'};...
                    'F_EP','LP_EP',{'odefile','x','p','ap'};...
                    'F_F','LP_LP',{'odefile','x','p','ap','(BranchingPars)'};...
                    'GH_H','GH_H',{'odefile','x','p','ap'};...
                    'GH_LPC','GH_LPC',{'odefile','xtot','p','s','ap','NTST','NCOL','amp','(BranchingPars)'};...
                    'GPD_PD','GPD_PD',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'H_EP','H_EP',{'odefile','x','p','ap'};...
                    'H_H','H_H',{'odefile','x','p','ap'};...
                    'H_LC','H_LC',{'odefile','x','p','ap','amp','NTST','NCOL'};...
                    'Het_Het','Het_Het',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HH_H','HH_H',{'odefile','x','p','ap'};...
                    'HH_NS1','HH_NS1',{'odefile','xtot','p','s','ap','NTST','NCOL','amp'};...
                    'HH_NS2','HH_NS2',{'odefile','xtot','p','s','ap','NTST','NCOL','amp'};...
                    'Hom_Hom','Hom_Hom',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'Hom_HSN','Hom_HSN',{'odefile','xtot','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HSN_Hom','HSN_Hom',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HSN_HSN','HSN_HSN',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HTHet_Het','HTHet_Het',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HTHet_HTHet','HTHet_HTHet',{'odefile','xtot','v','s','p','ap','UParam','ActiveUParam','SParam','ActiveSParam','NTST','NCOL','T','eps1','eps1Tol'};...
                    'HTHom_Hom','HTHom_Hom',{'odefile','x','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HTHom_HTHom','HTHom_HTHom',{'odefile','xtot','v','s','p','ap','UParam','ActiveUParam','SParam','ActiveSParam','NTST','NCOL','T','eps1','eps1Tol'};...
                    'HTHSN_HSN','HTHSN_HSN',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'HTHSN_HTHSN','HTHSN_HTHSN',{'odefile','xtot','v','s','p','UParam','ActiveUParam','SParam','ActiveSParam','NTST','NCOL','T','eps1','eps1Tol'};...
                    'LC_Hom','LC_Hom',{'odefile','xtot','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'LC_HSN','LC_HSN',{'odefile','xtot','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'LC_LC2','LC_LC',{'odefile','xtot','v','s','par','ap','NTST','NCOL'};...
                    'LPC_LC','LPC_LC',{'odefile','xtot','v','s','par','ap','NTST','NCOL'};...
                    'LPC_LPC','LPC_LPC',{'odefile','xtot','s','ap','NTST','NCOL','(BranchingPars)'};...
                    'LPNS_LPC','LPNS_LPC',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'LPNS_NS','LPNS_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'LPPD_PD','LPPD_PD',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'NCH_Hom','NCH_Hom',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'NCH_HSN','NCH_HSN',{'odefile','xtot','v','s','p','ap','NTST','NCOL','extravec','T','eps0','eps1'};...
                    'NS_LC','NS_LC',{'odefile','xtot','v','s','par','ap','NTST','NCOL'};...
                    'NS_NS','NS_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'P_EP','P_EP',{'odefile','x','p','ap','ndays'};...
                    'P_LC','P_LC',{'odefile','x','p','ap','NTST','NCOL','ndays','cycletol'};...
                    'PD_LC','PD_LC',{'odefile','xtot','s','NTST','NCOL','amp'};...
                    'PD_LC2','PD_LC2',{'odefile','xtot','s','ap','NTST','NCOL','amp'};...
                    'PD_PD','PD_PD',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'PDNS_NS','PDNS_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'PDNS_PD','PDNS_PD',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'R1__LPC','R1_LPC',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'R1__NS','R1_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'R2__NS','R2_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'R2__PD','R2_PD',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'R3__NS','R3_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'R4__NS','R4_NS',{'odefile','xtot','s','ap','NTST','NCOL'};...
                    'T_EP','BP_EP',{'odefile','x','p','s','amp','dumm:copy cds from curve'};...
                    'T_T','LP_LP',{'odefile','x','p','ap','(BranchingPars)'};...
                    'ZH_F','ZH_LP',{'odefile','x','p','ap'};...
                    'ZH_H','ZH_H',{'odefile','x','p','ap'};...
                    'ZH_NS','ZH_NS',{'odefile','xtot','p','s','ap','NTST','NCOL','amp'}};
            else
                %Map functions
                init_functs={'CP_F','LPm_LPm',{'mapfile','x','p','ap','niter','(BranchingPars)'};...
                    'EP_FP','FPm_FPm',{'mapfile','x','p','ap','niter'};...
                    'F_F', 'LPm_LPm',{'mapfile','x','p','ap','niter','(BranchingPars)'};...
                    'FP_1DMan', 'FPm_1DMan', {'x','p','niter'};...
                    'FP_FP','FPm_FPm', {'mapfile','x','p','ap','niter'};...
                    'GPD_LP2', 'GPD_LP2m',{'mapfile','amp','x','p','ap','niter','(BranchingPars)'};...
                    'Het_Het','Het_Het', {'mapfile','C','p','ap','niter'};...
                    'Hom_Hom','Hom_Hom', {'mapfile','C','p','ap','niter'};...
                    'LPPD_NS2', 'LPPD_NS2m',{'mapfile','amp','x','p','ap','niter'};...
                    'NS_NS', 'NSm_NSm', {'mapfile','x','p','ap','niter'};...
                    'NS_NS_Other', 'NSm_NSm_Other', {'mapfile','x','p','ap','amp','niter'};...
                    'NS_NS_Same', 'NSm_NSm_Same',  {'mapfile','x','p','ap','amp','niter'};...
                    'P_FP','P_FPm', {'mapfile','x','p','ap','niter','ndays'};...
                    'PD_FP2','PDm_FP2m', {'mapfile','x','p','s','amp','niter'};...
                    'PD_PD','PDm_PDm',  {'mapfile','x','p','ap','niter','(BranchingPars)'};...
                    'PDNS_NS2', 'PDNS_NS2m', {'mapfile','amp','x','p','ap','niter'};...
                    'R2__NS2','R2_NS2m',  {'mapfile','amp','x','p','ap','niter'};...
                    'R3__NS3', 'R3_NS3m', {'mapfile','amp','x','p','ap','niter','(BranchingPars)'};...
                    'R4__LP4m1', 'R4_LP4m1',{'mapfile','amp','x','p','ap','niter','(BranchingPars)'};...
                    'R4__LP4m2', 'R4_LP4m2',{'mapfile','amp','x','p','ap','niter','(BranchingPars)'};...
                    'R4__NS4', 'R4_NS4m',{'mapfile','amp','x','p','ap','niter','(BranchingPars)'};...
                    'T_F', 'LPm_LPm',  {'mapfile','x','p','ap','niter','(BranchingPars)'};...
                    'T_FP','BPm_FPm',{'mapfile','x','p','s','amp','niter','dumm:copy cds from curve'}}; %last param can also be direction of amp
            end
            ndx=strcmp([ptype '_' clabel],init_functs(:,1));
            %             %debug%%%%
            %             args={};
            %             for i=1:size(init_functs,1)
            %                 args =[args init_functs{i,2}{:}];
            %             end
            %             args=unique(args);
            %             g_grind.cont.inifun_args=args;
            %             %end debug %%%%
            argnames=init_functs(ndx,3);
            if ~isempty(argnames)
                argnames=argnames{1};
                init_fun=str2func(['init_' init_functs{ndx,2}]);
            else
                init_fun=[];
            end
        end
        
    end
    methods (Access = protected)  %for testing public must become private
        function ndx=getsamecurve(obj,curve)
            %same frompoint,same ctype, same freepars same p0
            %(except freepars)
            %can be overwritten to check other essential parameters:
            %Backward, niter
            ndx=obj.getsamecurve@grind_cont(curve);
            if any(ndx)
                fndx=find(ndx);
                for i=1:length(fndx)
                    [~,fdiffer]=struccmp(curve.data.settings,obj.curves(fndx(i)).data.settings);
                    if ~isempty(intersect(fdiffer,{'Backward','niter'}))
                        ndx(fndx(i))=false;
                    end
                end
            end
        end
        function setdefaults(obj)
            global g_grind;
            sett=obj.matcont_opt(:,[1 4]).';
            obj.settings.matcont=struct(sett{:});
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
            obj.settings.grind.symbolic=~g_grind.statevars.vector&&isempty(g_grind.syms.errormsg);
            obj.settings.grind.activepars=true(size(pars));
            obj.settings.derived.allpars=pars;
            obj.settings.derived.direction=[true true];
            obj.settings.derived.freepars=[];
            obj.settings.derived.parranges=struct('par',{},'range',[NaN NaN]);
        end
        
        function run_point(obj,frompoint,ctype)
            global g_grind cds;
            if obj.settings.derived.ismap
                h='M';
            else
                h='';
            end
            fprintf(['\nRunning <a href="https://sourceforge.net/projects/matcont">MatCont%s[%s]- Numerical Bifurcation Analysis Toolbox in Matlab</a>' ...
                '\nCitation: Dhooge, A., W. Govaerts, Yu.A. Kuznetsov, H.G.E. Meijer and B. Sautois, 2008 \nNew features of the software MatCont for bifurcation analysis of dynamical systems.\nMCMDS Vol. 14, No. 2, pp 147-175'...
                ' <a href="http://dx.doi.org/10.1080/13873950701742754">doi: 10.1080/13873950701742754</a> \n\n'],h,obj.settings.derived.version);
            
            ptype=obj.pointprops(frompoint.propndx).ptype;
            [initfun,argnames]=obj.get_init_fun(ptype,ctype);
            %determine if we need to go in two directions to fill the
            %parranges
            %           freepars=frompoint.p0(obj.settings.derived.freepars);
            %             ran=obj.settings.grind.parranges1;
            %             if length(freepars)==2
            %                 ran=[ran; obj.settings.grind.parranges2];
            %             end
            obj.settings.derived.direction=[true true]; %we cannot predict which direction forward is
            %obj.settings.derived.direction=[any(isnan(ran(:,1))|ran(:,1)<freepars) any(isnan(ran(:,2))|ran(:,2)>freepars)];
            if isempty(argnames)
                error('grind:matcont:unknowncurve','Point-curve combination is not allowed');
            end
            if obj.settings.grind.symbolic&&(isempty(g_grind.syms.Jacobian)||isempty(g_grind.syms.Jacobianp))
                try
                    if i_hastoolbox('symbolic')&&~g_grind.statevars.vector
                        enterjac('-s');
                    end
                catch %#ok<CTCH>
                    obj.settings.grind.symbolic = false;
                    warning('conteq:symfail', 'Cannot determine symbolic differentials, using numeric approximations instead');
                    %  g_grind.syms.Jacobian = {};
                    %  g_grind.syms.Jacobianp = {};
                    %  g_grind.syms.Hessian={};
                end
                obj.settings.grind.symbolic = isempty(g_grind.syms.errormsg);
            end
            propndx=obj.getndx('curveprops','ctype',ctype);
            codim=obj.curveprops(propndx).npars;
            if codim<0
                codim=2;
            end
            if codim>numel(obj.settings.derived.freepars)
                error('grind:matcont:freepars','Not enough free parameters selected')
            end
            args=cell(size(argnames));
            for i=1:length(argnames)
                switch argnames{i}
                    case {'odefile','mapfile'}
                        args{i}=@obj.handles;%matcontmodel;
                    case 'x'
                        if obj.settings.grind.sdjitter>0
                            args{i}=frompoint.x0+randn(size(frompoint.x0))*obj.settings.grind.sdjitter;
                        else
                            args{i}=frompoint.x0;
                        end
                    case 'T'
                        args{i}=frompoint.data.T/2;
                    case 'par'
                        for j=1:length(obj.curves)
                            for k=1:length(obj.curves(j).results.s)
                                if strcmp(frompoint.id,obj.curves(j).results.s(k).data.id)
                                    args{i}=obj.curves(j).results.s(k).data.parametervalues;
                                    break;
                                end
                            end
                        end
                    case 'dumm:copy cds from curve'
                        %dummy parameter for BP just copy cds
                        for j=1:length(obj.curves)
                            for k=1:length(obj.curves(j).results.s)
                                if strcmp(frompoint.id,obj.curves(j).results.s(k).data.id)
                                    cds= obj.curves(j).data.cds;
                                    break;
                                end
                            end
                        end
                        args(i)=[];
                    case 'xtot'
                        for j=1:length(obj.curves)
                            for k=1:length(obj.curves(j).results.s)
                                if strcmp(frompoint.id,obj.curves(j).results.s(k).data.id)
                                    args{i}=obj.curves(j).results.x;
                                    break;
                                end
                            end
                        end
                    case 'v'
                        for j=1:length(obj.curves)
                            for k=1:length(obj.curves(j).results.s)
                                if strcmp(frompoint.id,obj.curves(j).results.s(k).data.id)
                                    args{i}=obj.curves(j).results.v;
                                    break;
                                end
                            end
                        end
                    case 's'
                        % args{i}=frompoint; %if Matcont needs other fields than s.data this line should be changed
                        for j=1:length(obj.curves)
                            for k=1:length(obj.curves(j).results.s)
                                if strcmp(frompoint.id,obj.curves(j).results.s(k).data.id)
                                    args{i}=obj.curves(j).results.s(k);
                                    break;
                                end
                            end
                        end
                    case 'niterx2'
                        %this is with the initial function PDm_FP2m here
                        %the niters is automatically multiplied with 2
                        if isfield(frompoint.data,'niter')
                            obj.settings.derived.realiters(2)=obj.settings.grind.niter;
                            obj.settings.grind.niter=frompoint.data.niter;
                        end
                        obj.settings.derived.realiters(1)=obj.settings.grind.niter*2;
                        args{i}=obj.settings.grind.niter;
                    case 'p'
                        args{i}=frompoint.p0(obj.settings.grind.activepars);
                    case 'ap'
                        args{i}=obj.all2active(obj.settings.derived.freepars(1:codim));
                    case 'BranchingPar'
                        args{i}=frompoint.data.param;
                        ndx1=strcmp(argnames,'ap');
                        if numel(args{ndx1})==2
                            args(ndx1)={[args{ndx1} frompoint.data.param]};
                        end
                    case 'BranchingPars'
                        bp=obj.get('BranchingPars');
                        if isempty(bp)
                            error('grind:matcont','No branching pars defined, use set(''BranchingPars'',PAR) to set a branching parameter');
                        elseif ischar(bp)
                            args{i}=find(strcmp(obj.settings.derived.allpars,bp));
                        else
                            args{i}= bp;
                        end
                    case '(BranchingPars)' %optional
                        args{i}=obj.get('BranchingPars');
                    otherwise
                        %other arguments should be added as options
                        [args{i},found]=obj.get(argnames{i});
                        if ~found
                            error('grind:matcont:unsupported','Continuation not yet supported (argument: %s)',argnames{i});
                        end
                end
            end
            
            % initfun=str2func(sprintf('init_%s_%s',label,clabel));

            try
                [x0,v0]=initfun(args{:});
                if isempty(x0)
                    error('grind:matcont:initfun_fail','Error initialization function %s',func2str(initfun))
                end
            catch err1
                err2=MException(sprintf('matcont:%s',func2str(initfun)),'Error in matcont initialization function "%s"',func2str(initfun));
                disp(getReport(err1))
                i_waitbar([]);
                throw(addCause(err2,err1));
            end
            run_cont(obj,obj.curveprops(propndx).handle,x0,v0)
            if isfield(obj.settings.derived,'realiters')
                obj.settings.grind.niter=obj.settings.derived.realiters(2);
                obj.settings.derived=rmfield(obj.settings.derived,'realiters');
            end
        end
        function expand_curve(obj,frompoint,ctype,dirs)
            global cds;
            if nargin<4
                %dirs must be a vector of 2 logicals [true true] for both
                %directions [false true] for right only etc.
                dirs=obj.settings.derived.direction;
            end
            ndx=obj.getndx('curves','ctype',ctype)&obj.getndx('curves','frompoint',frompoint.id);
            if sum(ndx)>1
                curve.data.settings=obj.get('-nondefault');
                curve.freepars=obj.settings.derived.freepars(1:codim);
                curve.frompoint=frompoint.id;
                ndx=getndx(obj,'curves','samecurve',curve);
            end
            if sum(ndx)==1
                curveold=obj.curves(ndx);
                
                if exist('cpl.m','file')<=0
                    obj.open;
                end
                if isfield(curveold.data,'expanded')
                    curveold.data.expanded=curveold.data.expanded+1;
                else
                    curveold.data.expanded=1;
                end
                %expand to the right:
                curve1=curveold;
                if dirs(2)
                    cds=curveold.data.cds;
                    cds.options.Backward=0;
                    [curve1.results.x,curve1.results.v,curve1.results.s,curve1.results.h,curve1.results.f]=...
                        cont(curveold.results.x,curveold.results.v,curveold.results.s,...
                        curveold.results.h,curveold.results.f,cds);
                else
                    curve1=curveold;
                end
                %expand to the left:
                % curve1=flipcurve(curve1);
                curve2=curve1;
                if dirs(1)
                    cds.options.Backward=1;
                    [curve2.results.x,curve2.results.v,curve2.results.s,curve2.results.h,curve2.results.f]=...
                        cont(curve1.results.x,curve1.results.v,curve1.results.s,...
                        curve1.results.h,curve1.results.f,cds);
                else
                    curve2=curve1;
                end
                add_curve(obj,curve2);
                % obj.close;
            end
        end
        function run_cont(obj,varargin)
            
            %       global cds;
            %low level cont procedure (direct call to cont after settings)
            curve=struct('freepars',obj.settings.derived.freepars,'ctype',obj.settings.derived.ctype,'frompoint',obj.settings.derived.frompoint,...
                'data',[],'color','b','propndx',find(obj.getndx('curveprops','ctype',obj.settings.derived.ctype)),'results',struct('startndx',1,'x',[],'v',[],'s',[],'h',[],'f',[],'stabil',[]));
            %set current settings in matcont
            curve.color=obj.curveprops(curve.propndx).color;
            curve.freepars=curve.freepars(1:obj.curveprops(curve.propndx).npars);
            sett1=eval_settings(obj,curve);
            f=fieldnames(sett1);
            sett=contset;
            %             if obj.settings.derived.ismap&&obj.settings.derived.ndim<2&&strcmp(obj.settings.derived.ctype,'FP')
            %                 %we need to ignore Niemark Sacker in one dimensional model (else we get error)
            %                 sett1.IgnoreSingularity=[1 sett1.IgnoreSingularity];
            %             end
            for i=1:length(f)
                sett=contset(sett,f{i},sett1.(f{i}));
            end
            if obj.settings.derived.direction(1)
                curve1=docont(obj,curve,varargin,contset(sett,'Backward',1));
            else
                curve1=[];
            end
            try
                if obj.settings.derived.direction(2)
                    %                  sett = contset;
                    % sett = contset(sett,'Multipliers',1);
                    % %disp('>> opt = contset(opt,''Backward'',0); ');
                    % sett = contset(sett,'Backward',0);
                    % sett = contset(sett,'MaxNumPoints',300);
                    % %disp('>>opt = contset(opt,''Singularities'',1); ');
                    % sett = contset(sett,'Singularities',1);
                    curve2=docont(obj,curve,varargin,contset(sett,'Backward',0));
                else
                    curve2=[];
                end
            catch err
                %in case of error anyway add the backwards results
                if ~isempty(curve1)
                    add_curve(obj,curve1)
                end
                rethrow(err);
            end
            add_curve(obj,obj.merge_curves(curve1,curve2))
        end
        
        function [tspan,y0,options] = my_init(obj)
            global g_grind;
            han = obj.handles;
            y0=zeros(1,obj.settings.derived.ndim);
            options = odeset('Jacobian',han(3),'JacobianP',han(4),'Hessians',han(5),'HessiansP',han(6),'Vectorized',g_grind.solver.opt.Vectorized);
            tspan = [0 10];
        end
        
        function settings1 = eval_settings(obj, curve)
            %           relh = obj.settings.grind.relh; %#ok<NASGU>
            %             if isempty(obj.settings.grind.parranges1)
            %                 range = NaN;
            %             else
            %                 range = obj.settings.grind.parranges1;
            %                 range = abs(range(2) - range(1));
            %             end
            %
            %             if isnan(range)&&~isempty(obj.settings.derived.freepars)
            %                 range = abs(evalin('base', obj.settings.derived.allpars{obj.settings.derived.freepars(1)}));
            %             end
            %
            %             if isnan(range)||range==0
            %                 range = 1; %#ok<NASGU>
            %             end
            
            if ~isempty(obj.settings.matcont.Backward)
                if obj.settings.matcont.Backward
                    obj.settings.derived.direction=obj.settings.derived.direction & [true false];
                else
                    obj.settings.derived.direction=obj.settings.derived.direction & [false, true];
                end
            end
            %Automatic removal of singularities that are impossible
            settings1.IgnoreSingularity=obj.get('IgnoreSingularity');
            if ~isempty(settings1.IgnoreSingularity)
                ndx=curve.propndx;
                singlist=obj.curveprops(ndx).singularity;
                singlist=singlist(settings1.IgnoreSingularity);
                s=sprintf('%s,',singlist{:});
                fprintf('Ignored singlarities: %s\n',s(1:end-1));
            end
            
            settings1=obj.settings.matcont;
            f=fieldnames(settings1);
            for i=1:length(f)
                if ischar(settings1.(f{i}))
                    disp(settings1.(f{i}))
                    v=eval(settings1.(f{i}));
                    if ~isnan(v)
                        settings1.(f{i})=v;
                    end
                end
            end
            

            settings1.MinStepsize=settings1.MinStepsize*obj.settings.grind.StepAmp;
            settings1.MaxStepsize=settings1.MaxStepsize*obj.settings.grind.StepAmp;
            settings1.InitStepsize=settings1.InitStepsize*obj.settings.grind.StepAmp;
            if ~isempty(obj.settings.grind.Userhandles)&&numel(obj.settings.grind.Userhandles)==numel(obj.settings.matcont.UserfunctionsInfo)
                settings1.Userfunctions=numel(obj.settings.grind.Userhandles);
                for i=1:length(obj.settings.matcont.UserfunctionsInfo)
                    ndx=obj.getndx('pointprops','label',obj.settings.matcont.UserfunctionsInfo(i).label);
                    if ~any(ndx)
                        %mindim=-1 indicates that it is a user function
                        obj.pointprops(end+1)=struct('ptype',obj.settings.matcont.UserfunctionsInfo(i).label,'label',obj.settings.matcont.UserfunctionsInfo(i).label,...
                            'descr',obj.settings.matcont.UserfunctionsInfo(i).name,'codim',1,'mindim',-1,'ctypes',{{}});
                    end
                end
            else
                settings1.Userfunctions=0;
            end
        end
        
        
        function curve=merge_curves(obj,curve1,curve2)
            %             curve1=struct('freepars',obj.settings.derived.freepars,'ctype',obj.settings.derived.ctype,...
            %                 'p0',obj.settings.derived.p0,'x0',obj.settings.derived.x0,'settings',eval_settings(obj,curve),'results',struct('startndx',1,'x',[],'v',[],'s',[],'h',[],'f',[]));
            if isempty(curve1)||isempty(curve1.results.x)
                curve=curve2;
                return;
            elseif isempty(curve2)||isempty(curve2.results.x)
                curve=curve1;
                return;
            end
            if ~strcmp(curve1.ctype,curve2.ctype)||any(curve1.freepars~=curve2.freepars)
                error('grind:matcont:merge','Cannot merge incompatible types')
            end
            s1=obj.curveprops(curve1.propndx).resfun(obj,curve1);
            pars1=s1.parvalues(:,curve1.freepars(end));
            s1=obj.curveprops(curve2.propndx).resfun(obj,curve2);
            pars2=s1.parvalues(:,curve2.freepars(end));
            if curve1.results.startndx==1&&curve2.results.startndx==1
                if length(pars1)>1&&length(pars2)>1
                    diffpars1=pars1(2)-pars1(1);
                    diffpars2=pars2(2)-pars1(1);
                    if diffpars1<diffpars2
                        %diffpars1 = backwards
                        curve=curve1;
                        curvef=curve2;
                    else
                        curve=curve2;
                        curvef=curve1;
                    end
                else
                    curve=curve2;
                    curvef=curve1;
                end
                %Remove first and last singularities that will be in the middle
                %and adapt the indices
                curve=flipcurve(curve);
                curve.results.s=curve.results.s(1:end-1);
                n=size(curve.results.x,2);
                for i=1:length(curvef.results.s)
                    curvef.results.s(i).index=n+curvef.results.s(i).index;
                end
                curvef.results.s=curvef.results.s(2:end);
                
                %reverse curve
                curve.results.startndx=n;
                curve.results.x=[curve.results.x,curvef.results.x];
                curve.results.v=[-curve.results.v,curvef.results.v];
                curve.results.s=[curve.results.s;curvef.results.s];
                curve.results.h=[curve.results.h,curvef.results.h];
                curve.results.f=[curve.results.f,curvef.results.f];
            end
        end
        
        
        function add_curve(obj,curve)
            if isempty(curve.results.x)
                disp('Continuation failed');
                return;
            end
            % if strcmp(curve.ctype,'EP')
            p0=obj.points(obj.getndx('points','id',curve.frompoint)).p0;
            curve=add_stability(curve,obj,p0);
            % end
            %  for i=1:length(curve.results.s)
            %      curve.results.s(i).label=strtrim(curve.results.s(i).label);
            %  end
            ids=obj.add_points(curve);
            for i=1:length(ids)
                curve.results.s(i).data.id=ids{i};
            end
            ndx=find(obj.getndx('curves','samecurve',curve),1);
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
                point.propndx=find(obj.getndx('pointprops','label',regexp(spoint.id,'([A-Za-z+0-9]*_)|([A-Za-z+]*)','match','once')));
                lab=obj.pointprops(point.propndx).label;
                pntndx=[];
            else
                ndx=find(obj.getndx('points','x0+p0',{x0;p0}));
                pntndx=[];
                if ~isempty(ndx)
                    labs=obj.pointprops([obj.points(ndx).propndx]).label;
                    pntndx=find(strcmp(labs,strtrim(spoint.label)),1);
                    if ~isempty(pntndx)
                        pntndx=ndx(pntndx);
                        point=obj.points(pntndx);
                        return;
                    end
                end
                if ~isempty(regexp(spoint.label,'BP[0-9]', 'once'))
                    point.data.param=round(str2double(spoint.label(3:end)));
                    point.data.label=spoint.label;
                    apar=obj.settings.derived.allpars{point.data.param};
                    spoint.label='BP';
                    spoint.msg=sprintf('branch point (parameter "%s")',apar);
                    point.propndx=find(obj.getndx('pointprops','ptype','BP'));
                elseif strcmp(strtrim(spoint.label),'BP')
                    point.propndx=find(obj.getndx('pointprops','ptype','T'));
                else
                    point.propndx=find(obj.getndx('pointprops','label',strtrim(spoint.label)));
                end
                oldndx=find(obj.getndx('points','label',strtrim(spoint.label)));
                %find translation for label For instance LP->F
                if isempty(point.propndx)
                    disp(spoint);
                    error('grind:matcont','Unknown type of point, probably your MatCont version does not match (use updategrind to update)')
                end
                lab=obj.pointprops(point.propndx).ptype;
                point.id=[lab int2str(length(oldndx)+1)];
            end
            if isfield(point,'propndx')
                point.data.msg=obj.pointprops(point.propndx).descr;
            end
            switch lab
                case {'EP','P'}
                    point.data.msg=spoint.msg;
                case 'H'
                    if isfield(spoint.data,'lyapunov')&&spoint.data.lyapunov>0
                        point.id=sprintf('H+%d',length(oldndx)+1);
                    else
                        point.id=sprintf('H%d',length(oldndx)+1);
                    end
                case 'NS'
                    if obj.settings.derived.ismap
                        c=spoint.data.c;
                    else
                        c=spoint.data.nscoefficient;
                    end
                    if real(c)>0
                        point.id=sprintf('NS+%d',length(oldndx)+1);
                    else
                        point.id=sprintf('NS%d',length(oldndx)+1);
                    end
                case 'PD'
                    if obj.settings.derived.ismap
                        if isfield(obj.settings.derived,'realiters')
                            point.data.niter=obj.settings.derived.realiters(1);
                        elseif isfield(obj.settings.grind,'niter')
                            point.data.niter=obj.settings.grind.niter;
                        end
                        point.data.msg=sprintf('Period doubling from %d to %d',point.data.niter,point.data.niter*2);
                        b=spoint.data.b;
                    else
                        b=spoint.data.pdcoefficient;
                    end
                    if real(b)<0
                        point.id=sprintf('PD+%d',length(oldndx)+1);
                    else
                        point.id=sprintf('PD%d',length(oldndx)+1);
                    end
                    
            end
        end
        
        
        
        function init_engine(obj,remove)
            if nargin==1
                remove=false;
            end
            if obj.settings.derived.ismap
                matcontfile='matcontm.m';
            else
                matcontfile='matcont.m';
            end
            if xor(obj.settings.derived.opened,remove)
                return;
            end
            matcont_in_path=exist(matcontfile,'file')>0;
            if ~matcont_in_path&&~isempty(obj.settings.derived.engine_path)||~exist(obj.settings.derived.engine_path,'dir')
                obj.settings.derived.engine_path=obj.findupdate([],'',obj.settings.derived.ismap);
            end
            
            if obj.settings.derived.ismap
                fullpath={ 'Continuer', 'FixedPointMap','Systems','LimitPointMap',...
                    'PeriodDoublingMap', 'NeimarkSackerMap','MultilinearForms','AD',...
                    'Heteroclinic', 'HeteroclinicT', 'Homoclinic', 'HomoclinicT',...
                    'InvManifolds',fullfile('Testruns','CodStock'), ...
                    fullfile('Testruns','LeslieGower'), fullfile('Testruns','Tnfmap'),...
                    fullfile('Testruns','Connections'),fullfile('Testruns','InvManifolds'),...
                    fullfile('GUI','Imported'), ''};
            else
                fullpath={'','Continuer','Equilibrium', 'LimitCycle','LimitCycleCodim2',...
                    'PeriodDoubling','Systems','LimitPoint','Hopf', 'LimitPointCycle',...
                    'NeimarkSacker','BranchPoint', 'BranchPointCycle', 'Homoclinic',...
                    'HomoclinicSaddleNode', 'HomotopySaddle', 'HomotopySaddleNode', 'HomotopyHet',...
                    'Heteroclinic','MultilinearForms','Help','SBML','GUI'};
            end
            if ~remove
                for i=1:length(fullpath)
                    addpath(fullfile(obj.settings.derived.engine_path,fullpath{i}));
                end
            elseif remove
                for i=1:length(fullpath)
                    rmpath(fullfile(obj.settings.derived.engine_path,fullpath{i}));
                end
            end
        end
    end
    
    properties (Hidden= true,Access=public)
        %Options:
        %name,validation,description,matcont default
        
        %          MoorePenrose: []
        %         SymDerivative: []
        %        SymDerivativeP: []
        %         TestFunctions: []
        %             WorkSpace: []
        %              Locators: []
        %          ActiveParams: []
        %         ActiveUParams: []
        %         ActiveSParams: []
        %          ActiveSParam: []
        matcont_opt={'InitStepsize','n>0#s','the initial stepsize (default: 0.01)',0.01;...
            'MinStepsize','n>0#s','the minimum stepsize to compute the next point on the curve (default: 1e-5)',1e-5;...
            'MaxStepsize','n>0#s','the maximum stepsize (default: 0.1)',0.1;...
            'MaxCorrIters','i>0','maximum number of correction iterations (default: 10)',10;...
            'MaxNewtonIters','i>0','maximum number of Newton-Raphson iterations before switching to Newton-Chords in the corrector iterations (default: 3)',3;...
            'MaxTestIters','i>0','maximum number of iterations to locate a zero of a testfunction (default: 10)',10;...
            'Increment','n>0','the increment to compute the derivatives numerically (default: 10-5)',1e-5;...
            'FunTolerance','n>0','tolerance of function values: ||F(x)|| £ FunTolerance is the first convergence criterium of the Newton iteration (default: 10-6)',1e-6;...
            'VarTolerance','n>0','tolerance of coordinates: ||dx|| £ VarTolerance is the second convergence criterium of the Newton iteration (default: 10-6)',1e-6;...
            'TestTolerance','n>0','tolerance of test functions (default: 10-5)',1e-5;...
            'Singularities','l','boolean indicating the presence of a singularity matrix (default: 1)',true;...
            'MaxNumPoints','i>0','maximum number of points on the curve (default: 300)',300;...
            'Backward','l#E','boolean indicating the direction of the continuation (sign of the initial tangent vector) v0, []=two directions (default: [])',[];...
            'CheckClosed','i>=0','number of points indicating when to start to check if the curve is closed (0 = do not check) (default: 50)',50;...
            'Adapt','i','number of points after which to call the adapt-function while computing the curve (default: 1=adapt always)',1;...
            'IgnoreSingularity','i#c#s','vector containing indices of singularities which are to be ignored (default: empty)','';...
            'Multipliers','l','boolean indicating the computation of the multipliers (default: 1)',true;...
            'Eigenvalues','l','boolean indicating the computation of the eigenvalues (default: 0)',true;...
            'Userfunctions','l','boolean indicating the presence of user functions (default: false)',false;...
            'UserfunctionsInfo','r#E','is an array with structures containing information about the userfunctions (default: empty)',[];...
            'PRC','s','variable indicating the computation of the phase response curve (default: empty)',[];...
            'dPRC','s','variable indicating the computation of the derivative of the phase response curve (default: empty)',[];...
            'Input','n','vector representing the input given to the system for the computation of the phase response curve (default: 0)',0};
        %            'TSearchOrder','l','Search order of the tangent vector increasing=0, decreasing=1 (default=0)',false;...

        first_argument_par='ActiveSParam'; %this is the first of the options that are context specifif

        grind_opt={'StepAmp','n>0','Amplifier for stepsize settings (default: 1)',1;...
            'par1','p','first free parameter (default: '''')','';...
            'parranges1','n&length(n)==2','range for the first parameter (default: [NaN NaN])',[NaN NaN];...
            'par2','p#E','second free parameter (default: '''')','';...
            'par3','p#E','third free parameter (default: '''')','';...
            'parranges2','n&length(n)==2','range for the second parameter (default: [NaN NaN])',[NaN NaN];...
            'mindist','n>0','minimum distance between different special points (default: 1E-5)',1E-5;...
            'stateranges','n','min/max value for each of the state variables, or one row if all are the same (default: [Nan NaN])',[NaN NaN];...
            'symbolic','l','if available use the symbolic toolbox for Jacobians (default 1)',true;...
            'sdjitter','n>=0','standard deviation of jitter added to equilibrium points (default: 0)',0;...
            'Userhandles','f#E','Handles to user functions (default: [])',[];...
            'activepars','l#c','indices of the active parameters (default: all 1)','true(size(allpars))';...
            ...%the next options are only needed for some curves (see get_init_fun)
            'ActiveSParam','i#s','Active stable parameters for HTHSN curves(comma delimited) (default: '')','';...
            'ActiveUParam','i#s','Active unstable parameters for HTHSN curves(comma delimited) (default: '')','';...
            'BranchingPars','p#E','List of indexes of branching parameters (default: [])',[];...
            'NTST','i>0','number of mesh intervals (Default: 40)',40;...
            'NCOL','i>0','number of collocation nodes (default 4)',4';...
            'SParam','i#s','Free stable parameter for HTHSN curves(comma delimited) (default: '')','';...
            'TTolerance','n','Tolerance of period for homoclinic bifurcation (BT->Hom) (default: 1e-5)',1e-5;...
            'UParam','i#s','Free unstable parameter for HTHSN curves(comma delimited) (default: '')','';...
            'cycletol','n>0','tolerance for finding a limit cycle from an initial point (default: 1e-2)',1E-2;...
            'amp','n','Initial amplitude used in sevaral initial points (BP->EP/H->LC) (default: 1E-6)',1E-6;...
            'eps0','n','First tolerance parameter for homoclinc bifurcation (Hom) (default: 0.01)',0.01;...
            'eps1','n','Second tolerance parameter for homoclinc bifurcation (Hom) (default: 0.01)',0.01;...
            'eps1Tol','n','Tolerance for eps1 - HTHSN curves (default: 1e-2)',1e-2;...
            'extravec','l&length(l)==3','Vector indicating which parameters should be adapted for homoclinic bifurcation (Hom) (default [0 1 1])',[false true true];...
            'ndays','n>0','number of days to stabilize an initial point (backwards or forwards) (default: 1000)',1000};
        
        %ctype npars args handle descr
        curveprops=struct('ctype',{},'npars',[],'color',[],'handle',[],'descr',[])
        % id descr codim ctypes
        pointprops=struct('ptype',{},'label',[],'descr',[],'codim',[],'mindim',[],'ctypes',[]);
    end
end

function res = get_curve_var(obj,curve,varname,index)
if isnumeric(curve)
    curve=obj.curves(curve);
end
if isfield(curve.data.cds,'nUserf')
    nUserf=curve.data.cds.nUserf;
else
    nUserf=0;
end
switch varname
    case 'timestep'
        res=curve.results.h(1,:);
    case 'niters'
        res=curve.results.h(2,:);
    case 'userfunct'
        if nUserf==0
            res=[];
        else
            res=curve.results.h(3:3+nUserf-1,:);
            if index<=nUserf
                res=res(index,:);
            end
        end
    case 'testfunct'
        res=curve.results.h(3+nUserf:end,:);
        if ischar(index)
            index=find(strcmp(index,obj.curveprops(curve.propndx).singularity));
        end
        if ~isempty(index)
            index=curve.data.cds.ActTest==index;
            if any(index)
                res=res(index,:);
            end
        end
end
end

% what is in results.h?
% cds.h=current stepsize
% it= num of Newton iterations
% testvals are the testvals of the curve
% uservals are the user functions
%   if Singularities & ~Userfunctions
%       hout(:,i) = [cds.h;it;cds.testvals(3-cds.atv,:)'];
%   elseif Userfunctions & ~Singularities
%       hout(:,i) = [cds.h;it;cds.uservals(3-cds.utv,:)'];
%   elseif Userfunctions & Singularities
%       hout(:,i) = [cds.h;it;cds.uservals(3-cds.utv,:)';cds.testvals(3-cds.atv,:)'];
%   else
%       hout(:,i) = [cds.h;it];
%   end


function [s,xtra]=extract_curve_eq(obj,curve)
ndim=obj.settings.derived.ndim;
xtra.stabil=[];
%all parameters inculding the free pars
p0new=obj.points(obj.getndx('points','id',curve.frompoint)).p0;
p0=p0new(:,ones(1,size(curve.results.x,2)+1));
x=curve.results.x;
x(:,end+1)=NaN;
xtra.stabil=curve.results.stabil;
xtra.stabil(1,end+1)=-1;
p0(curve.freepars,:)=x(ndim+1:ndim+numel(curve.freepars),:);
xtra.stabil=transpose(xtra.stabil);
s=struct('parvalues',transpose(p0),'pars',{obj.settings.derived.allpars},'Y',permute(transpose(x(1:ndim,:)), [3 2 1]),'t',zeros(1,1,size(x,2)),'perm',[]);
end

function [s,xtra]=extract_curve_cycle(obj,curve)
if isfield(curve.data.curveds,'ncoords')
    ndim=curve.data.curveds.ncoords;
end
x=curve.results.x;
p0=repmat(obj.points(obj.getndx('points','id',curve.frompoint)).p0,[1,size(x,2)]);
p0(curve.freepars,:)=x(end-numel(curve.freepars)+1:end,:);
if size(x,1)-numel(curve.freepars)>ndim
    periods=x(ndim+1:ndim+1,:);
else
    periods=zeros(1,size(x,2))+curve.results.s(1).data.T;
end
nstatevar=obj.settings.derived.ndim;
x(ndim+1:ndim+nstatevar,:)=NaN;
YY=permute(reshape(x(1:ndim+nstatevar,:),nstatevar,ndim/nstatevar+1,size(x,2)),[2,1,3]);
%YY=YY(:,:,1:end-1);
%p0=p0(:,1:end-1);
%xtra.stabil=curve.results.stabil(1:end-1);
xtra.stabil=repmat(curve.results.stabil,[size(YY,1),1]);
xtra.stabil=permute(xtra.stabil,[1 3 2]);
tt=zeros(size(YY,1),1,size(YY,3));
for i=1:size(YY,3)
    tt(:,1,i)=linspace(0,periods(i),size(YY,1));
end
s=struct('parvalues',transpose(p0),'pars',{obj.settings.derived.allpars},'Y',YY,'t',tt,'perm',[]);
end

function curve=add_stability(curve,obj,p0)
%recalculate eigenvalues to have all information
    function stabil=getstability(eigenvals,isdiffer)
        stabil=zeros(1,size(eigenvals,2),'int16');
        iscomplex=any(imag(eigenvals)>0,1);
        if isdiffer
            isstable = max(abs(eigenvals), [], 1) < 1;
            issaddle = (min(abs(eigenvals), [], 1) < 1) & (max(abs(eigenvals), [], 1) > 1); %No convert
        else
            isstable=max(real(eigenvals),[],1)<=0;
            issaddle=~isstable&min(real(eigenvals),[],1)<0;
        end
        stabil(isstable&~iscomplex)=1; %stable point
        stabil(isstable&iscomplex)=2; %stable spiral
        stabil(issaddle)=3;      %saddle point
        stabil(~isstable&~issaddle&~iscomplex)=4; %unstable node
        stabil(~isstable&~issaddle&iscomplex)=5; %unstable spiral
    end
global g_grind;
if isempty(curve.results.x)
    curve=[];
    return;
end
ndim=obj.settings.derived.ndim;
if any(strcmp(curve.ctype,{'EP','FP'}))
    x=curve.results.x(1:ndim,:);
    jac_handle=i_getodehandle('coco_jac',[],iif(isempty(g_grind.syms.errormsg),'','numonly'));
    p0=p0(:,ones(size(x,2),1));
    for i=1:size(curve.freepars)
        p0(curve.freepars(i),:)=curve.results.x(ndim+i,:);
    end
    eigenvals = zeros(size(x));
    if strcmp(g_grind.solver.opt.Vectorized, 'on')
        Js = jac_handle(x,p0);
        for i = 1:size(x, 2)
            eigenvals(:, i) = eig(Js(:, :, i));
        end
    else
        for i = 1:size(x, 2)
            J = jac_handle(x(:,i),p0(:,i));
            eigenvals(:, i) = eig(J);
        end
    end
    curve.results.stabil=getstability(eigenvals,g_grind.solver.isdiffer);
elseif any(strcmp(curve.ctype,{'NS','NS1','NS2','NS3','NS4','NS_Same','NS_Other'}))
    %first test function is positive for subcritical
    %    c=curve.results.h(3+obj.settings.matcont.Userfunctions,:);
    c=get_curve_var(obj,curve,'testfunct','CH');
    curve.results.stabil=zeros(size(c),'int16')+6;
    curve.results.stabil(c<0)=7;
    curve.results.stabil(c==111|c==500)=-1; %neutral saddle curve
elseif any(strcmp(curve.ctype,{'PD','PD2'}))
    %last test function is negative for subcritical
    % b=curve.results.h(end,:);
    b=get_curve_var(obj,curve,'testfunct','GPD');
    curve.results.stabil=zeros(size(b),'int16')+7;
    curve.results.stabil(b<0)=6;
elseif strcmp(curve.ctype,'H')
    %first lyapunov coefficient (determining sub/supercritical
    %hopf) negative is supercritical
    %6=sub and 7=supercit
    lyap=get_curve_var(obj,curve,'testfunct','GH');
    %    lyap=curve.results.h(end,:);
    curve.results.stabil=zeros(size(lyap),'int16')+6;
    curve.results.stabil(lyap<0)=7;
    curve.results.stabil(lyap==500)=-1; %neutral saddle
else
    mult=grind_matcont.multipliers(curve);
    if ~isempty(mult)
        curve.results.stabil=getstability(mult,true);
    else
        curve.results.stabil=[];
    end
end
end

function curve=flipcurve(curve)
n=size(curve.results.x,2);
for i=1:length(curve.results.s)
    curve.results.s(i).index=n+1-curve.results.s(i).index;
end

%reverse curve
curve.results.startndx=n;
curve.results.x=fliplr(curve.results.x);
curve.results.v=fliplr(-curve.results.v);
curve.results.s=flipud(curve.results.s);
curve.results.s(1).label= '00';
curve.results.s(1).msg= 'This is the first point of the curve';
curve.results.s(end).label= '99';
curve.results.s(1).msg= 'This is the last point of the curve';
curve.results.h=fliplr(curve.results.h);
curve.results.f=fliplr(curve.results.f);
end



function plotjac(fid,jac,separ,func)
fprintf(fid,'\n\n%s\n%s\n',separ,func);
if iscell(jac)
    siz=size(jac);
    jac= regexprep(jac,'g_X1[\(]([0-9]*)[\,][\:][\)]','g_X1($1)');
    if numel(siz)==2
        fprintf(fid,'g_X2=zeros(%d,%d);\n',siz);
        for i=1:siz(1)
            for j=1:siz(2)
                if ~any(strcmp(jac(i,j),{'0','0.0'}))
                    fprintf(fid,'g_X2(%d,%d) = %s;\n',i,j,jac{i,j});
                end
            end
        end
    else
        fprintf(fid,'g_X2=zeros(%d,%d,%d);\n',siz);
        for i=1:siz(1)
            for j=1:siz(2)
                for k=1:siz(3)
                    if ~any(strcmp(jac(i,j,k),{'0','0.0'}))
                        fprintf(fid,'g_X2(%d,%d,%d) = %s;\n',i,j,k,jac{i,j,k});
                    end
                end
            end
        end
    end
elseif isstruct(jac)
    ssiz=sprintf('%d,',jac.size);
    fprintf(fid,'g_X2=zeros(%s);\n',ssiz(1:end-1));
    for i=1:numel(jac.unique)
        sjac= regexprep(jac.unique(i).equation,'g_X1[\(]([0-9]*)[\,][\:][\)]','g_X1($1)');
        if numel(jac.unique(i).indices)==1
            fprintf(fid,'g_X2(%d) = %s;\n',jac.unique(i).indices,sjac);
        else
            sind=sprintf('%d,',jac.unique(i).indices);
            fprintf(fid,'g_X2([%s]) = %s;\n',sind(1:end-1),sjac);
        end
    end
end


end
function curve1=docont(obj,curve1,args,sett1)
global cds;
try
    [curve1.results.x,curve1.results.v,curve1.results.s,curve1.results.h,curve1.results.f]=cont(args{:},sett1);
    curve1.data.cds=cds;
    curve1.data.settings=obj.get('-nondefault');
    if ~isempty(obj.curveprops(curve1.propndx).curveds)
        curveds=obj.curveprops(curve1.propndx).curveds;
        eval(sprintf('global %s',curveds));
        curve1.data.curveds=eval(curveds);
    end
    if obj.curveprops(curve1.propndx).iscycle
        if strcmp(curve1.ctype,'LC')&&numel(curve1.results.s)>2
            curve1.results.s(end).label='LC';
        end
        %  curve1.data.lds=lds;
    end
    if size(curve1.results.f,2)>size(curve1.results.x,2)
        %sometimes f has one point more
        curve1.results.f=curve1.results.f(:,size(curve1.results.x,2));
    end
catch err1
    %as matcont does not use identifiers we glue an
    %identifier to the error
    disp(getReport(err1));
    i_waitbar([]);
    err2=MException('matcont:cont','Continuation error in matcont');
    throw(addCause(err2,err1))
end
i_waitbar(1);
end