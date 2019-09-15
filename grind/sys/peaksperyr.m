function npeaks=peaksperyr(yearlength,ivar,rel_height)
global g_Y g_t;
if nargin<2
   ivar=1;
end
if nargin<3
   rel_height=1;
end
ndays=round(g_t(length(g_t))-g_t(1));
data=interp1(g_t,g_Y(:,ivar),1:ndays);
nyrs=floor(ndays/yearlength);
%peaks=zeros(nyrs,1);
peakk=zeros(100,1);
mink=zeros(100,1);
k=2;
imin=1;
imax=1;
for j=1:nyrs
   for i=3:yearlength-2
      k=k+1;
      if (data(k-2)<data(k))&&(data(k-1)<data(k))&&(data(k+1)<data(k))&&(data(k+2)<data(k))
         peakk(imax)=k;
         imax=imax+1;    
      %   disp(sprintf('%g  : %g %g %g %g %g',[k,data(k-2),data(k-1),data(k),data(k+1),data(k+2)]));
      end
      if (data(k-2)>data(k))&&(data(k-1)>data(k))&&(data(k+1)>data(k))&&(data(k+2)>data(k))
         mink(imin)=k;
         imin=imin+1;
      %   disp(sprintf('%g  : %g %g %g %g %g',[k,data(k-2),data(k-1),data(k),data(k+1),data(k+2)]));
      end
   end
   k=k+4;
end
if rel_height>0.999
   npeaks=imax/nyrs;
else
   npeaks=0;
   for i=1:imax-1
      mindiff=9999;
      for j=i-1:i+1
         if (j>0)&&(j<imin)&&(abs(peakk(i)-mink(j))<mindiff)
            mindiff=abs(peakk(i)-mink(j));
            reldiff=data(mink(j))/data(peakk(i));
         end
      end
      if (mindiff<9999)&&(reldiff<rel_height)
         npeaks=npeaks+1;
      end   
   end
   npeaks=npeaks/nyrs;
end


