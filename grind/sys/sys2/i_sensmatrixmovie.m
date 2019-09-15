function i_sensmatrixmovie(s_matrix, vlabels, plabels)
   times=cell2num(vlabels(:,2));
   utimes=unique(times);
   if isempty(utimes)||size(vlabels,2)<2
       return;
   end

   vtexts=vlabels(strcmp(num2str(utimes(1)),vlabels(:,2)),1);
   amatrix=zeros(length(utimes),length(plabels),length(vtexts));
   drawmovie=1;
for i=1:length(utimes)
   ndx=times==utimes(i);
   if sum(ndx)==length(vtexts)
        amatrix(i,:,:)=s_matrix(:,ndx);
   else
       drawmovie=0;
   end

end

if drawmovie
    i_matrixmovie(amatrix,vtexts,plabels);
else
    disp('cannot make matrix movie as there are missing values');
end

