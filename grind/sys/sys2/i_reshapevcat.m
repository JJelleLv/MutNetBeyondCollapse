function m=i_reshapevcat(siz,varargin)
if siz(end)>1
    flen1=find(cellfun('size',varargin,2)==1);
    if ~isempty(flen1)
        for i=1:length(flen1)
            varargin{flen1(i)}=varargin{flen1(i)}(ones(1,siz(end)));
        end
    end
end
m=reshape(vertcat(varargin{:}),siz);
end
