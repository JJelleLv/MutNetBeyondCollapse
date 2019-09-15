function ok=check_conn(M)
%CHECK_CONN - check connections of network
%check if unipartite network is fully connected. Returns true if full
%connections exist.
%Alogrithm: just try all connections

%set diagonal to 0
M(1:(size(M,1)+1):end)=0;
degree=full(sum(M,2));
if any(degree==0)
    %there cannot be any without links
    ok=false;
    return;
end

issymmetric=all(all(triu(M)==triu(transpose(M))));
conn=transpose(false(size(degree)));
if issymmetric
    conn(1)=true;
    ok=checknode(1,M,conn);
else
    for i=1:length(conn)
        conn=transpose(false(size(degree)));
        conn(i)=true;
        ok=checknode(i,M,conn);
        if ~ok
            return;
        end

    end

end

function [ok,conn]=checknode(i,M,conn)
%Check connections of node i (recursive)
%stop if all conn are true (or if all links are visited)
if all(conn)
    ok=true;
else
    ok=false;
    connectsto=find(M(i,:)>0&~conn);
    conn(connectsto)=true;
    for j=1:length(connectsto)
        [ok,conn]=checknode(connectsto(j),M,conn);
        if ok
            return;
        end

    end

end

