function success=i_split_sde
%split stochastic model in drift and diffusion part (dwiener)
global g_grind
success=isfield(g_grind.solver,'dwiener')||isfield(g_grind.solver,'djump');
if success && (~isfield(g_grind.syms,'drift')||isempty(g_grind.syms.drift))
    %extracts the lines with state variables from the model\
    %remove comments
    model = regexprep(g_grind.model,'%.*','');
    vars=regexp(model,'[A-Za-z_][A-Za-z0-9_]*(?=[\(:0-9,\)]*''[ =])','match','once');
    ndx=~cellfun('isempty',vars);
    vars=vars(ndx);
    model1=regexp(model(ndx),'(?<=[=]).*','match','once');
    for i=1:length(model1);
        mod1= model1{i};
        while mod1(end)==';';
            mod1=mod1(1:end-1);
        end
        f=strfind(mod1,'%');
        if ~isempty(f)
            mod1=mod1(1:f(1)-1);
        end
        model1{i}=mod1;
    end
    model=parsed_equation(model1);
    funcs=str2cell(g_grind.funcs);
    if any(strcmp(funcs,'end;')|strcmp(funcs,'end'))
        error('grind:dwiener:split','Cannot split this model (for stochast_heun) as there are compound statements in the auxiliary variables (e.g. if else for case etc)')
    end
    for i=length(funcs):-1:1
        model=model.subs(funcs{i});
    end
    maxeq=0;
    for i=1:length(model)
        f=find(strcmp('dwiener',model(i).fs));
        if length(f)>maxeq
            maxeq=length(f);
        end
    end
    dwien=cell(length(model),maxeq);
    dwien(:,:)={'0'};
    corrs=cell(length(model),maxeq);
    jcorrs=corrs;
    dwien_dx=cell(length(model),maxeq);
    dwien_dx(:,:)={'0'};
    drift1=model;
    for i=1:length(model)
        s=model(i).structure;
        ndx=find(strcmp({s(:).oper},'dwiener'));
        for j=1:length(ndx)
            left=s(ndx(j)).leftvar;
            for k=1:length(s)
                found=any(strcmp(s(k).args,left));
                if found&&~any(strcmp(s(k).oper,{'+','-'}))
                    error('grind:dwiener','%s\nCannot split this model (for stochast_heun): dwiener should be in an additive term',char(model(i)));
                end
            end
        end
    end
    if isfield(g_grind.solver,'djump')
        for i=1:length(model)
            f=find(strcmp('djump',drift1(i).fs));
            for j=length(f):-1:1
                [pars,jj]=drift1(i).getfunctionpars(f(j));
                pars=regexprep(pars,'^''|''$','');
                args=djump(pars{:});
                if isfield(args,'corr')
                    corrs{i,j}=args.corr;
                else
                    corrs{i,j}=0;
                end
                fs=[drift1(i).fs(1:f(j)-1) {'0'} drift1(i).fs(jj(end,2)+2:end)];
                drift1(i).equation=sprintf('%s',fs{:});
            end
        end
        jumps=cell(length(model),max([g_grind.solver.djump.args(:).num]));
        for i=1:length(model)
            if ~g_grind.statevars.vector
                siz=[1 1];
            else
                siz=[g_grind.statevars.dims{i}.dim1 g_grind.statevars.dims{i}.dim2];
            end
            jumps(i,:)={sprintf('zeros(%d,%d)',siz(1),siz(2))};
            for j=1:length(g_grind.solver.djump.args)
                if g_grind.solver.djump.args(j).var==i
                    jumppars=[g_grind.solver.djump.args(j).timing_pars,g_grind.solver.djump.args(j).size_pars];
                    if g_grind.statevars.vector
                    else
                        for k1=1:length(g_grind.statevars.names)
                            f1=strcmp(jumppars,g_grind.statevars.names{k1});
                            if any(f1)
                                jumppars(f1)={sprintf('g_X1(%d,:)',k1)};
                            end
                        end
                    end
                    parslist=sprintf('%s,',jumppars{:});
                    if ~isempty(parslist)
                        parslist=parslist(1:end-1);
                    end
                    jumps{i,g_grind.solver.djump.args(j).num}=sprintf('i_djump(%d,t,%s)',j,parslist);
                end
            end
        end
        g_grind.solver.djump.jumps=jumps;
    for i=1:numel(corrs)
        siz=size(corrs);
        if length(siz)==2
            siz(3)=1;
        end
        if all(siz==[g_grind.statevars.dim,g_grind.statevars.dim,size(jumps,2)])
            g_grind.solver.djump.corr=corrs;
        end
    end
    if ~isfield(g_grind.solver.djump,'corr')||isempty(g_grind.solver.djump.corr)
        totcorrs=zeros(g_grind.statevars.dim,g_grind.statevars.dim,size(jumps,2));
        L=totcorrs;
        if g_grind.statevars.vector
            for j=1:size(jumps,2)
                for i=1:size(corrs,1)
                    d=g_grind.statevars.dims{i};
                    ndx=d.from:d.to;
                    if ischar(corrs{i,j})
                        cor=evalin('base',corrs{i,j});
                    else
                        cor=corrs{i,j};
                    end
                    if isempty(cor)
                        cor=0;
                    elseif cor==1
                        cor=1-1E-10;
                    end
                    totcorrs(ndx,ndx,j)=cor;
                end
                totcorrs(:,:,j)=setdiagon(totcorrs(:,:,j),1);
            end
            for j=1:size(jumps,2)
                L(:,:,j)=chol(totcorrs(:,:,j));
            end
            g_grind.solver.djump.corr=totcorrs;
            g_grind.solver.djump.L=L;
        end
    end

    end
    
    for i=1:length(model)
        f=find(strcmp('dwiener',drift1(i).fs));
        for j=length(f):-1:1
            [pars,jj]=drift1(i).getfunctionpars(f(j));
            pars=regexprep(pars,'^''|''$','');
            args=dwiener(pars{:});
            dwien{i,j}=args.fun;
            if isfield(args,'dfundx')
                dwien_dx{i,j}=args.dfundx;
            else
                dwien_dx{i,j}='#';
            end
            if isfield(args,'corr')
                corrs{i,j}=args.corr;
            else
                corrs{i,j}=0;
            end
            fs=[drift1(i).fs(1:f(j)-1) {'0'} drift1(i).fs(jj(end,2)+2:end)];
            drift1(i).equation=sprintf('%s',fs{:});
        end
    end
    if any(any(strcmp(dwien_dx,'#')))
        dwien_dx=[];
    end
    drift1=cell(drift1);
    dim=length(vars);
    twodim=false;
    for i = 1:dim
        %replace statevar symbols with g_X1(i)
        if g_grind.statevars.vector
            h=sprintf('zeros(%d,%d)+',g_grind.statevars.dims{i}.dim1,g_grind.statevars.dims{i}.dim2);
            if g_grind.statevars.dims{i}.dim2==1
                newvar=sprintf('g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to)',i,i);
            else
                twodim=true;
                newvar=sprintf('reshape(g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to),g_grind.statevars.dims{%d}.dim1,g_grind.statevars.dims{%d}.dim2)',i,i,i,i);
            end
            drift1{i}=[h drift1{i}];
            for l=1:size(dwien,2)
                dwien{i,l}=[h dwien{i,l}];
            end
            if ~isempty(dwien_dx)
                for l=1:size(dwien_dx,2)
                    dwien_dx{i,l}=[h dwien_dx{i,l}];
                end
            end
        else
            newvar=sprintf('g_X1(%d,:)', i);
        end
        drift1=regexprep(drift1,sprintf('\\<%s\\>',vars{i}),newvar);
        dwien=regexprep(dwien,sprintf('\\<%s\\>',vars{i}),newvar);
        if ~isempty(dwien_dx)
            dwien_dx=regexprep(dwien_dx,sprintf('\\<%s\\>',vars{i}),newvar);
        end
    end
    if twodim
        for i=1:size(drift1,1)
            drift1{i}=sprintf('reshape(%s,[%d,1])',drift1{i},g_grind.statevars.dims{i}.dim1*g_grind.statevars.dims{i}.dim2);
        end
        for i=1:size(dwien,1)
            for j=1:size(dwien,2)
                dwien{i,j}=sprintf('reshape(%s,[%d,1])',dwien{i,j},g_grind.statevars.dims{i}.dim1*g_grind.statevars.dims{i}.dim2);
            end
        end
        for i=1:size(dwien_dx,1)
            for j=1:size(dwien_dx,2)
                dwien_dx{i,j}=sprintf('reshape(%s,[%d,1])',dwien_dx{i,j},g_grind.statevars.dims{i}.dim1*g_grind.statevars.dims{i}.dim2);
            end
        end
    end
    g_grind.syms.drift=drift1;
    g_grind.syms.diffusion=dwien;
    g_grind.syms.diffusion_dx=dwien_dx;
    g_grind.solver.dwiener.corr=[];
    for i=1:numel(corrs)
        siz=size(corrs);
        if length(siz)==2
            siz(3)=1;
        end
        if all(siz==[g_grind.statevars.dim,g_grind.statevars.dim,maxeq])
            g_grind.solver.dwiener.corr=corrs;
        end
    end
    if isempty(g_grind.solver.dwiener.corr)
        totcorrs=zeros(g_grind.statevars.dim,g_grind.statevars.dim,maxeq);
        L=totcorrs;
        if g_grind.statevars.vector
            for j=1:maxeq
                for i=1:size(corrs,1)
                    d=g_grind.statevars.dims{i};
                    ndx=d.from:d.to;
                    if ischar(corrs{i,j})
                        cor=evalin('base',corrs{i,j});
                    else
                        cor=corrs{i,j};
                    end
                    if isempty(cor)
                        cor=0;
                    elseif cor==1
                        cor=1-1E-10;
                    end
                    totcorrs(ndx,ndx,j)=cor;
                end
                totcorrs(:,:,j)=setdiagon(totcorrs(:,:,j),1);
            end
            for j=1:maxeq
                L(:,:,j)=chol(totcorrs(:,:,j));
            end
            g_grind.solver.dwiener.corr=totcorrs;
            g_grind.solver.dwiener.L=L;
        end
    end
    
end
