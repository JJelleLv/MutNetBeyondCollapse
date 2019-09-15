function [DirCHmat]=func_CHdirVect(Vect_mat)

DATAlength=length(Vect_mat(1,:));

DirCHmat=nan(1,DATAlength-1);
for DATANR=1:(DATAlength-1)
    
    Vect1=Vect_mat(:,DATANR);
    Vect2=Vect_mat(:,(DATANR+1));
    CosTheta=dot(Vect1,Vect2)/(norm(Vect1)*norm(Vect2));
    ThetaInDegrees=acosd(CosTheta);
    
    DirCHmat(1,DATANR)=ThetaInDegrees;
end