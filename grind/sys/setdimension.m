%SETDIMENSION   Change the dimension of a (state) variable
%   This function changes the dimension of a vector/matrix state variable or parameter. It changes
%   state variables dynamically but it doesn’t automatically check for necessary changes in vector
%   parameters. However, you can first indicate how the dimensions of the parameters are linked
%   to the dimensions of the state variables. Only these linked parameters are updated automatically 
%   if the size of a state variable changes. If a parameter is a scalar variable (=one value) it 
%   is never changed.
%   Note that if dimensions are changed during a run, model results can get lost if the size of a 
%   state variable is shrinking.
%
%
% Usage:
%    SETDIMENSION('VAR',[NROWS,NCOLS]) - set the dimension of variable/parameter VAR to NROWS and NCOLS
%    SETDIMENSION('VAR',NROWS) - assumes NCOLS=1, or NCOLS=NROWS if the variable is already square
%    SETDIMENSION PAR [VAR1(1) VAR2(2)] CREATEFUN DIAGON  - links the dimensions
%       of parameter PAR to the state variable VAR, VAR(1) is the first dimension for variable VAR etc.
%       CREATEFUN is a function to create a new element, DIAGON is a special function/value for the diagonal.
%    SETDIMENSION PAR1 PAR2 PAR3 VAR  - the list of parameters PAR1..PAR3 have the same dimension as statevariable VAR.
%    SETDIMENSION('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'afterscript' [string] - script that will be run after changing the dimensions.
%     'createfun' [string or function_handle] - script that is run for resizeing a parameter (default mean of the previous values)
%     'diagon' [number or empty] - Set diagional with a value (default = NaN)
%     'pars' [parameter+] - list of parameters that are associated with a state variable
%     'parsize' [string] - state variable used set parameter value, can be defined as [X(1) Y(1)] (first dim of X, first dim of Y)
%     'size' [integer>0 and length(integer)<=2] - new size of the state variable
%     'var' [state variable] - state variable to change
%     'vector' [logical] - set vector notation on or off (used internally)
%   SETDIMENSION('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' AFTERSCRIPT - defines a script that will be run each time after changing the dimensions.
%     '-l' - Lists the defined dimensions of the parameters.
%     '-v' VECTOR - set vector notation on or off
%
%
%  See also model, definespace
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setdimension')">commands setdimension</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function res = setdimension(varargin)
%(avar, dims, createfun ,diagon)
global g_Y g_grind;
fieldnams={'afterscript', 's', 'script that will be run after changing the dimensions.','';...
   'diagon', 'n#E', 'Set diagional with a value (default = NaN)',NaN;...
   'pars', 'p+', 'list of parameters that are associated with a state variable','';...
   'parsize', 's', 'state variable used set parameter value, can be defined as [X(1) Y(1)] (first dim of X, first dim of Y)',[];...
   'createfun', 's#f', 'script that is run for resizeing a parameter (default mean of the previous values)',[];...
   'size', 'i>0&length(i)<=2', 'new size of the state variable',[10 1];...
   'var', 'v', 'state variable to change',i_statevars_names(1);...
   'vector', 'l', 'set vector notation on or off (used internally)',true}';
args=i_parseargs(fieldnams,...
'if(hasoption(''-a'')),deffields=''afterscript'';elseif(hasoption(''-v'')),deffields=''vector'';elseif(argtype(1,''v'')),deffields=''var,size'';else,deffields=''pars(+),parsize,createfun,diagon'';end','-a,-l,-v',varargin);
%   'par', 'p#E', 'name of parameter to set size','';...

if strcmp(g_grind.model{1},'%external odefile')
    error('grind:setdimension:externalode','Cannot use this command when an external odefile is used');
end
if isfield(args,'afterscript')
    g_grind.onsetdimension=args.afterscript;
end
if any(strcmp(args.opts,'-v'))||isfield(args,'vector') %-vector on or -vector off - set vector modus of GRIND on or off
    if ~isfield(args,'vector')
        args.vector=true;
    end
    [ggrind,fs]=i_analysemodel(g_grind.model);
    beforeeq=strcmp(fs,'=');
    for i=1:size(beforeeq,1)
        f1=find(beforeeq(i,:),1);
        if ~isempty(f1)
            beforeeq(i,1:f1(1))=true;
        end
    end
    if args.vector
    %vectorize (if necessary)
    if ~ggrind.statevars.vector
       if ~any(any(strcmp(fs,'.*')|strcmp(fs,'./')|strcmp(fs,'.^')))
           if ~isempty(g_grind.scheme)
               ndx=find(strncmp(g_grind.scheme,'%exp=',4));
               for k=1:length(ndx)
                   s=g_grind.scheme{ndx(k)};
                   s=vectorize(s(6:end));
                   g_grind.scheme{ndx(k)}=sprintf('%%exp=%s',s);
               end
           end
           fs(strcmp(fs,'*'))={'.*'};
           fs(strcmp(fs,'/'))={'./'};
           fs(strcmp(fs,'^'))={'.^'};
           fs(strcmp(fs,'||'))={'|'};
           fs(strcmp(fs,'&&'))={'&'};
           %replace || and && with | and &?
       end
       %add (1:2) after statevars 
       v=strcmp(fs,':')&beforeeq;
       if any(any(v))
            for i=1:size(fs,1)
              if any(v(i,:))
                  f1=find(strcmp(fs(i,:),'(')&beforeeq(i,:));
                  if ~isempty(f1)
                  f2=find(strcmp(fs(i,:),')')&beforeeq(i,:),1);
                  fs(i,f1+1:f2)={''};
                  fs{i,f1}='(1:2,1:2)';
                  end
              end
            end
            g_grind.model=transpose(ggrind.locfunc);
            h=length(g_grind.model);
            for i=1:size(fs,1)
              g_grind.model{i+h}= sprintf('%s',fs{i,:});
           end
       else
       if ggrind.solver.isdiffer
           v=strcmp(fs,'t')&beforeeq;
           for i=1:size(fs,1)
              if any(v(i,:))
                  f1=find(strcmp(fs(i,:),'('),1);
                  fs{i,f1}='(1:2,1:2)(';
              end
              g_grind.model{i}= sprintf('%s',fs{i,:});
           end
       else
           v=strcmp(fs,'''')&beforeeq;
           fs(v)={'''(1:2,1:2)'};
           g_grind.model=transpose(ggrind.locfunc);
           h=length(g_grind.model);
           for i=1:size(fs,1)
              g_grind.model{i+h}= sprintf('%s',fs{i,:});
           end
       end
       end
       i_makemodel;
    end
    else
        %remove vector definitions if necessary
        f=find(strcmp('%vec=1',g_grind.scheme));
        if ~isempty(g_grind.scheme)&& ~isempty(f)
           g_grind.scheme(f)=[];
        end
        if g_grind.statevars.vector
           v=strcmp(fs,':')&beforeeq;
           for i=1:size(fs,1)
              if any(v(i,:))
                  f1=find(strcmp(fs(i,:),'(')&beforeeq(i,:),1);
                  f2=find(strcmp(fs(i,:),')')&beforeeq(i,:),1);
                  fs(i,f1:f2)={''};
              end
              g_grind.model{i}= sprintf('%s',fs{i,:});
           end
           i_makemodel;
        end
    end
    return; 
end
if any(strcmp(args.opts,'-l'))
    if isfield(g_grind,'pardim')
        par('modelonly')
        disp('Parameter dimension       createfun      diagonal value');
        for i=1:length(g_grind.pardim)
        fprintf('%5s     [%s(%d),%s(%d)]     %s    %g\n',g_grind.pardim(i).name,...
            g_grind.pardim(i).dim1.var,g_grind.pardim(i).dim1.dim,g_grind.pardim(i).dim2.var,...
            g_grind.pardim(i).dim2.dim,g_grind.pardim(i).createfun,g_grind.pardim(i).diagon);
        end
    else
        disp('no parameter dimensions specified');
    end
    if isfield(g_grind,'onsetdimension')
        fprintf('onsetdimension = "%s"\n',g_grind.onsetdimension);
    end
    return;
end
    
if (~isfield(args,'var')||~isfield(args,'size'))&&~isfield(args,'pars')&&~any(strcmp(args.opts,'-a'))
   error('GRIND:setdimension:ArgError','Not enough arguments for setdimension');
end

if isfield(args,'pars')&&isfield(args,'parsize')
    %define how the size of parameter(s) should be scaled
    if ~isfield(g_grind,'pardim')
        g_grind.pardim=struct('name',{},'dim1',struct('var','','dim',0),'dim2',struct('var','','dim',0),'createfun',[],'diagon',NaN);
    end
    vars=regexp(args.parsize,'[A-Za-z_][A-Za-z0-9_]*','match');
    vardim=regexp(args.parsize,'[0-9]*','match');
    if length(vars)==1&&isempty(vardim)
        dim1=struct('var',vars{1},'dim',1);
        dim2=struct('var',vars{1},'dim',2);
    elseif length(vars)==2&&length(vardim)==2
        dim1=struct('var',vars{1},'dim',str2double(vardim{1}));
        dim2=struct('var',vars{2},'dim',str2double(vardim{2}));
    else
        error('grind:setdimension:parsize','Error in the size of the parameter');
    end
    if ~isfield(args,'diagon')
        args.diagon=NaN;
    end
    for i=1:length(args.pars)
        if ~isfield(args,'createfun')
           createfun=sprintf('mean(%s(:))',args.pars{i});
        else
           createfun=args.createfun;
        end
        f=find(strcmp({g_grind.pardim(:).name},args.pars{i}),1);
        if isempty(f)
            f=numel(g_grind.pardim)+1;
        end
        g_grind.pardim(f)=struct('name',args.pars{i},'dim1',dim1,'dim2',dim2,'createfun',createfun,'diagon',args.diagon);
    end
elseif isfield(args,'var')&&~isempty(args.var)&&isfield(args,'size')
   if numel(args.size)==1
       args.size=[args.size 1];
   end
   siz=evalin('base',sprintf('size(%s)',args.var));
   dims = args.size;
   if prod(dims)>1&&~g_grind.statevars.vector
       setdimension('-vector','on');
   end
   if numel(dims)==1
      if siz(1) == siz(2)
         dims(2) = dims(1);
      else
         if siz(1) ~= 1
             dims(2) = 1;
         else
             dims(2) = dims(1);
             dims(1) = 1;
         end
     end
   end
   dim1 = dims(1);
   dim2 = dims(2);
   if isempty(g_grind)||~g_grind.statevars.vector
      setdimvar(args.var, dim1, dim2, siz);
      if nargout > 0
         res = evalin('base', args.var);
      end
      return;
   end
   ivect = 0;
   for i = 1:size(g_grind.statevars.vectnames, 2)
      if strcmp(args.var, g_grind.statevars.vectnames{i})
         ivect = i;
         break;
      end
   end
   if ivect == 0
      setdimvar(args.var, dim1, dim2, siz);
      if nargout > 0
         res = evalin('base', args.var);
      end
      return;
   end
   ggrind = g_grind;
   try
      nams = g_grind.funcnames.names;
      g_grind.funcnames = [];
      g_grind.funcnames.names = nams;
      differ = dim1 * dim2 - g_grind.statevars.dims{ivect}.dim1 * g_grind.statevars.dims{ivect}.dim2;
      g_grind.statevars.dims{ivect}.dim1 = dim1;
      g_grind.statevars.dims{ivect}.dim2 = dim2;
      g_grind.statevars.dims{ivect}.to = g_grind.statevars.dims{ivect}.to + differ;
      g_grind.statevars.dim = g_grind.statevars.dim + differ;
      for i = ivect + 1:size(g_grind.statevars.vectnames, 2)
         g_grind.statevars.dims{i}.from = g_grind.statevars.dims{i}.from + differ;
         g_grind.statevars.dims{i}.to = g_grind.statevars.dims{i}.to + differ;
      end
      g_grind.statevars.names = {};
      if ~isempty(g_Y)  && (size(g_Y, 2) < g_grind.statevars.dim)
         g_Y(size(g_Y, 1), g_grind.statevars.dim) = 0;
      end
      setdimvar(args.var, dim1, dim2, siz);
      if nargout > 0
         res = evalin('base', args.var);
      end
      if isfield(g_grind,'pardim')
          for i=1:length(g_grind.pardim)
            if strcmp(args.var,g_grind.pardim(i).dim1.var)||strcmp(args.var,g_grind.pardim(i).dim2.var)
                s=evalin('base',sprintf('size(%s)',g_grind.pardim(i).name));
                if prod(s)>1
                   setdimvar(g_grind.pardim(i).name,getdim(g_grind.pardim(i).dim1.var,g_grind.pardim(i).dim1.dim),getdim(g_grind.pardim(i).dim2.var,g_grind.pardim(i).dim2.dim),...
                     s, g_grind.pardim(i).createfun, g_grind.pardim(i).diagon);
                end
            end
          end
      end
      if isfield(g_grind,'onsetdimension')&&~isempty(g_grind.onsetdimension)
          evalin('base',g_grind.onsetdimension);
      end
   catch err
      g_grind = ggrind; 
      rethrow(err);
   end
         %update g_grind.model with the new dimensions
      [ggrind,fs]=i_analysemodel(g_grind.model);
      beforeeq=strcmp(fs,'=');
      for i=1:size(beforeeq,1)
         f1=find(beforeeq(i,:),1);
         if ~isempty(f1)
            beforeeq(i,1:f1(1))=true;
         end
      end
      if ~isempty(g_grind.scheme)&& ~any(strcmp('%vec=1',g_grind.scheme))
          g_grind.scheme=[{'%vec=1'} g_grind.scheme];
      end
      for j=1:length(g_grind.statevars.vectnames)
          ndx=strcmp(fs,g_grind.statevars.vectnames{j})&beforeeq;
          dim=g_grind.statevars.dims{j};
          for i=1:size(fs,1)
              if any(ndx(i,:))
                  f1=find(strcmp(fs(i,:),'('),1);
                  if ~isempty(f1)
                     f2=find(strcmp(fs(i,:),')'),1);
                     fs(i,f1+1:f2)={''};
                     if dim.dim2==1
                         fs{i,f1}=sprintf('(1:%d)',dim.dim1);
                     else
                         fs{i,f1}=sprintf('(1:%d,1:%d)',dim.dim1,dim.dim2);
                     end
                  end
              end
          end
          hasrow=changescheme(g_grind.statevars.vectnames{j}, '%row=', int2str(dim.dim1));
          if ~hasrow
             addscheme(g_grind.statevars.vectnames{j},'%row=', int2str(dim.dim1));
          end
          hascol=changescheme(g_grind.statevars.vectnames{j}, '%col=', int2str(dim.dim2));
          if ~hascol
             addscheme(g_grind.statevars.vectnames{j},'%col=', int2str(dim.dim2));
          end
      end
           g_grind.model=transpose(ggrind.locfunc);
           h=length(g_grind.model);
           for i=1:size(fs,1)
              g_grind.model{i+h}= sprintf('%s',fs{i,:});
           end
          fid=fopen(which(g_grind.odefile),'r');
          odefile=textscan(fid,'%s','delimiter','\n','whitespace','');
          odefile=odefile{1};
          fclose(fid);
          if isfield(ggrind.solver,'dwiener')
              g_grind.solver.dwiener=ggrind.solver.dwiener;
          end
        if all(strcontains(odefile,'reshape('))
             i_makemodel; 
       end
end


%out('-defaults');
function res=getdim(avar,dim)
s=evalin('base',sprintf('size(%s)',avar));
res=s(dim);
function setdimvar(avar, dim1, dim2, siz, createfun, diagon)
if nargin < 5
   createfun = [];sprintf('mean(%s(:))', avar);
end
if nargin < 6
   diagon = nan;
end
if (dim1 - siz(1)==0)&&(dim2 - siz(2)==0)
   return;
end
differ = dim1 * dim2 - prod(siz);
if differ < 0
   if (dim1<=siz(1))&&(dim2<=siz(2))
      evalin('base',sprintf('%s=%s(1:%d,1:%d);',avar,avar,dim1,dim2));
   else
      evalin('base',sprintf('%s=reshape(%s(1:%d),%d,%d);',avar,avar,dim1*dim2,dim1,dim2));
   end
elseif differ == 0
   evalin('base',sprintf('%s=reshape(%s,%d,%d);',avar,avar,dim1,dim2));
   else %if (dim1<=siz(1))&&(dim2<=siz(2))
   oldvar =  evalin('base', avar);
   newvar = zeros(dim1, dim2)+mean(oldvar(:));
   if ~isempty(createfun)
      for i = siz(1):dim1
         for j = 1:dim2
            newvar(i, j) = evalin('base', createfun);
         end
      end      
      for i = 1:dim1
         for j= size(2):dim2
            newvar(i, j) = evalin('base', createfun);
         end
      end
   end
   if (dim1>=siz(1))&&(dim2>=siz(2))
      newvar(1:siz(1), 1:siz(2)) = mean(oldvar(:));
   else
      newvar(1:prod(siz)) = oldvar(:);
   end
   if ~isnan(diagon)
      newvar = setdiagon(newvar, diagon);
   end
   assignin('base', avar, newvar)
   %   evalin('base',sprintf('g_h321=%s;%s=zeros(%d,%d)+mean(mean(%s));%s(1:%d,1:%d)=g_h321;clear(''g_h321'');',avar,avar,dim1,dim2,avar,avar,siz(1),siz(2)));
   %else
   %   evalin('base',sprintf('g_h321=%s;%s=zeros(%d,%d)+mean(mean(%s));%s(1:%d)=g_h321;clear(''g_h321'');',avar,avar,dim1,dim2,avar,avar,siz(1)*siz(2)));
end

function addscheme(pname, plabel, pnew)
global g_grind;
if isfield(g_grind, 'scheme')&&~strcmp(pnew,'[]')
   k=find(strcmp(g_grind.scheme,sprintf('%%sym=%s',pname)));
  % k=1;
   while k < length(g_grind.scheme)
      s = g_grind.scheme{k};
      if strncmp(s, '%sym=', 5) && strcmp(s(6:end), pname)
         if (k < length(g_grind.scheme)) &&~(strncmp(g_grind.scheme{k}, '%[', 2))
            g_grind.scheme = {g_grind.scheme{1:k}, sprintf('%s%s',plabel,pnew),g_grind.scheme{k + 1:end}};
         end
         return;
      end
      k = k + 1;
   end
end

function found = changescheme(pname, plabel, pnew)
global g_grind;
found = 0;
if isfield(g_grind, 'scheme')&&~strcmp(pnew,'[]')
    k=find(strcmp(g_grind.scheme,sprintf('%%sym=%s',pname)));
    if ~isempty(k)
        while (k < length(g_grind.scheme)) &&~(strncmp(g_grind.scheme{k}, '%[', 2) || strncmp(g_grind.scheme{k}, plabel, 5))
            k = k + 1;
        end
        if (k < length(g_grind.scheme))  && strncmp(g_grind.scheme{k}, plabel, 5)
            s = g_grind.scheme{k};
            found = 1;
            g_grind.scheme{k} = [s(1:5), pnew];
        end
    end
end
