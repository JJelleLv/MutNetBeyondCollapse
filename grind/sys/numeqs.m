function eqs=numeqs(eqs)
%add a random number to each column to see the difference between [0 10 0] and [0 0 10]
eqs1=bsxfun(@times,eqs,rand(1,size(eqs,2)));
%determine the distance with the first element
dist=sum((bsxfun(@minus,eqs1,rand(1,size(eqs,2)))).^2,2);
[~,ndx]=sort(dist,'descend');
eqs1=eqs(ndx,:);
%distance between subsequent points
d=sum((eqs1(1:end-1,:)-eqs1(2:end,:)).^2,2);
d1=sort(d,'descend');
i=2;
while cv(d1(i:end))>1E-10
    i=i+1;
end
crit=d1(i)*1.5;
eqs=eqs1([true;d>crit],:);
size(eqs)
function res=cv(x)
res=var(x)/mean(x);



        
