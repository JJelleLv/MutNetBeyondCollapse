function s=i_statevars_names(i)
global g_grind
if nargin==0
   if g_grind.statevars.vector
        s=cell(1,g_grind.statevars.dim);
        k=1;
        for j=1:length(g_grind.statevars.vectnames)
            %very fast:
            nams=allelems(g_grind.statevars.vectnames{j},[g_grind.statevars.dims{j}.dim1 g_grind.statevars.dims{j}.dim2]);
            nams=nams(:);
            s(k:k+length(nams)-1)=nams;
            k=k+length(nams);
        end
    else
        s=g_grind.statevars.names;
   end
else
    if length(i)>1
        s={};
        for k=1:length(i)
            s=[s {i_statevars_names(i(k))}]; %#ok<AGROW>
        end

        return;
    end

    if i<=0||i>g_grind.statevars.dim
        s='';
        return;
    end

    if g_grind.statevars.vector
        j=1;
        while (j<=length(g_grind.statevars.dims))&&(g_grind.statevars.dims{j}.to<i)
            j=j+1;
        end

        i=i-g_grind.statevars.dims{j}.from+1;
        [i1,i2]=ind2sub([g_grind.statevars.dims{j}.dim1,g_grind.statevars.dims{j}.dim2],i);
        if g_grind.statevars.dims{j}.dim2>1
            s=sprintf('%s(%d,%d)',g_grind.statevars.vectnames{j},i1,i2);
        else
            s=sprintf('%s(%d)',g_grind.statevars.vectnames{j},i1);
        end

    else
        if length(i)==1
            s=g_grind.statevars.names{i};
        else
            s=sprintf('%s;',g_grind.statevars.names{i});
        end
    end

end

