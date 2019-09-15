function m=i_runsymstruc(s,siz2,varargin)
%global g_grind;
%make sure all elements have the correct size
if siz2>1
    flen1=find(cellfun('size',varargin,2)==1);
    if ~isempty(flen1)
        for i=1:length(flen1)
            varargin{flen1(i)}=varargin{flen1(i)}(ones(1,siz2));
        end
    end
end
%s=g_grind.syms.(symfield);
if siz2==1
    m=zeros(s.size);
    if isfield(s,'unique')
        for i=1:length(s.unique)
            m(s.unique(i).indices)=varargin{i};
        end
    end
else
    m=zeros([s.size,siz2]);
    if isfield(s,'unique')
        for i=1:length(s.unique)
            for j=1:length(s.unique(i).indices)
                c=index2sub(s.size,s.unique(i).indices(j));
                    m(c{:},:)=varargin{i};
            end
        end
    end
end
function a=index2sub(siz,index)
a=cell(size(siz));
[a{1:end}]=ind2sub(siz,index);


