%fast generation of list of elements of a matrix
%either one argument (which is the variable)
%or allelems(varname, siz);
function elems=allelems(p,siz)
if nargin==1
    siz=size(p);
    p=inputname(1);
end

if length(siz)==1
    siz=[siz siz];
end

if prod(siz)==1
    elems={p};
else
    [r,c]=meshgrid(1:siz(2),1:siz(1));
    if siz(2)==1
        ss=sprintf([p,'(%d);'],c(:));
    else
        ss=sprintf([p,'(%d,%d);'],transpose([c(:),r(:)]));
    end

    elems=regexp(ss(1:end-1),';','split');
    elems=reshape(elems,siz);
end

