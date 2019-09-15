function nyears=yrcycles(yearlength)
global g_Y g_t;
ndays=round(g_t(length(g_t))-g_t(1));
data=interp1(g_t,g_Y(:,1),1:ndays);
nyrs=floor(ndays/yearlength);
yrdiff=zeros(nyrs,1);
for i=2:nyrs
   k=(i-1)*yearlength+1;
   yrdiff(i)=sum(data(1:yearlength)-data(k:k+yearlength-1));
end
yrdiff=yrdiff/sum(data(1:yearlength));
nyears=1;
for i=2:nyrs
   sumdiff=0;
   n=0;
   for j=nyears+1:nyears:nyrs
      sumdiff=sumdiff+yrdiff(j);
      n=n+1;
   end
   if sumdiff/n<0.001
      break;
   end
   nyears=nyears+1;
end

