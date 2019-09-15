function [tr, datar] = i_concatdata(t1, data1, t2, data2)
if isempty(t1)
   tr = t2;
   datar = data2;
   return;
end

if isempty(t2)
   tr = t1;
   datar = data1;
   return;
end


ts = sort([t1; t2]);
ndx=~isnan(ts);
if ~all(ndx)
    ts=ts(ndx);
    data1=data1(~isnan(t1),:);
    data2=data2(~isnan(t2),:);
    t1=t1(~isnan(t1));
    t2=t2(~isnan(t2));
end
l1 = size(data1, 2);
l2 = size(data2, 2);
datar = zeros(length(ts), l1 + l2) + NaN;
if any(diff(t1)<0)||any(diff(t2)<0) %x may decrease
    tr=[t1;t2];
    datar(1:length(t1),1:l1)=data1;
    datar(length(t1)+1:end,l1 + 1:l1 + l2)=data2;
    return;
end
tr=zeros(length(ts),1)+NaN;
k = 1;
i = 1;

while i <= length(ts)
   tt = ts(i);
   i1=find(t1 == tt);
   n1=length(i1);
   if n1>0
      datar(k:k+n1-1, 1:l1) = data1(i1, :);
   end
   i2=find(t2 == tt);
   n2=length(i2);
   if n2>0
      datar(k:k+n2-1, l1 + 1:l1 + l2) = data2(i2, :);
   end
   tr(k:k+max(n1,n2)-1)=ts(i);
   i=i+n1+n2;
   k=k+max(n1,n2);
end

tr = tr(1:k-1);
datar = datar(1:k-1, :);
