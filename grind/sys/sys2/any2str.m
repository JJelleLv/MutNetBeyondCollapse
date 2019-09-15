% res=any2str(val,precision)
% any class is converted to a string that is executable 
% supported classes: numeric classes,logical,cell,char,struct,function_handle
% precision is the number of digits (see mat2str) default=15
%
function res = any2str(aval,precision)
if nargin<=1
    precision=15;
end
if isnumeric(aval)||islogical(aval)
    if isa(aval,'double')||isa(aval,'logical')
       res=mat2str(aval,precision);
       if strcontains(res,'+i*')
           res=regexprep(res,'+i[*]([0-9.eE-])*','+$1i');
       end
    else
       res=mat2str(aval,precision,'class');
    end
elseif  ischar(aval)||isa(aval,'string')
     if isempty(aval)&&sum(size(aval))>0
         siz=sprintf('%d,',size(aval));
         res=sprintf('char(zeros(%s))',siz(1:end-1));
     else
         if ischar(aval)
             delim='''';
         else
             delim='"';
         end
         aval=strrep(aval,delim,[delim delim]);
         new_line=char(10); %#ok<CHARTEN>
         if ~strcontains(char(aval),new_line)
             res=[delim aval delim];
         else
            res=sprintf('sprintf(%s%s%s)',delim,strrep(strrep(aval,new_line,'\n'),'%','%%'),delim);
         end
     end 
elseif isstruct(aval)
    f=fieldnames(aval);
    vals=cell(2,numel(f));
    if numel(aval)==1
        for i=1:length(f)
            vals{2,i}=any2str(aval.(f{i}),precision);
            if iscell(aval.(f{i}))&&~isempty(aval.(f{i}))
                vals{2,i}=['{' vals{2,i} '}'];
            end
            vals{1,i}=['''' f{i} ''''];
        end
    else
        for i=1:length(f)
            s={aval(:).(f{i})};
            vals{2,i}=any2str(s,precision);
            vals{1,i}=['''' f{i} ''''];
        end
    end
    res=sprintf('%s,',vals{:});
    res=sprintf('struct(%s)',res(1:end-1));
elseif  isa(aval,'function_handle')
     res=func2str(aval);
     if ~strcontains(res,'@')
         res=['@' res];
     end
elseif  iscell(aval)
    vals=cell(size(aval));
    for i=1:numel(vals)
        vals{i}=any2str(aval{i},precision);
    end
    doreshape=false;
    if isempty(aval)
        if sum(size(aval))==0 
            res='{} ';
        else
            siz=sprintf('%d,',size(vals));
            res=sprintf('cell(%s) ',siz(1:end-1));
        end
    elseif size(vals,1)==1
       res=sprintf('%s,',vals{:});
    elseif size(vals,2)==1
       res=sprintf('%s;',vals{:});
    else
       res=sprintf('%s,',vals{:});
       siz=sprintf('%d,',size(vals));
       doreshape=true;
    end
    if doreshape
        res=sprintf('reshape({%s},%s) ',res(1:end-1),siz(1:end-1));
    else
        res=['{' res(1:end-1) '}'];
    end
else
    res=char(aval);
end



end

