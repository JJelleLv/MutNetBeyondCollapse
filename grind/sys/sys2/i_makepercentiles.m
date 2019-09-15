function Result=i_makepercentiles(A,percentiles)                   
if max(percentiles)>2
   percentiles=percentiles/100;
end

Result=zeros(size(A,1),length(percentiles));
for i=1:size(A,1)
   Result(i,:)=percents(A(i,:),percentiles);
end

%function p=i_makepercentiles(A1,percents);
function p=percents(A1,percents)
A1=transpose(sort(A1));
m=length(A1);
while (m>1)&&isnan(A1(m))
   m=m-1;
end

q = (0.5:m - 0.5)./m;
A1 = [min(A1); A1(1:m); max(A1)];
q = [0 q 100];
p = interp1(q,A1,percents);

   
