function A=rep_elem(A,elem,ndx)
    if numel(ndx)==numel(elem)
        A(ndx)=elem;
    else
       A(ndx(1),ndx(2))=elem;
    end

