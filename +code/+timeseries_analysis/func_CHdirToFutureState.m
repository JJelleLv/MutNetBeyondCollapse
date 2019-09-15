function [DirCHtoFutureStateMat,DirSIMILARITYtoFutureStateMat]=func_CHdirToFutureState(Vect_mat,Vect_comp)

DATAlength=length(Vect_mat(1,:));
DIMs=length(Vect_mat(:,1));

DirCHtoFutureStateMat=nan(1,DATAlength);
DirSIMILARITYtoFutureStateMat=nan(1,DATAlength);
for DATANR=1:DATAlength
    
    %% determine angle (in degrees)
    Vect1=Vect_mat(:,DATANR);
    
    CosTheta=dot(Vect1,Vect_comp)/(norm(Vect1)*norm(Vect_comp));
    ThetaInDegrees=acosd(CosTheta);
    
    DirCHtoFutureStateMat(1,DATANR)=ThetaInDegrees;
    
    %% similarity measure
    [~,CDFvalue]=code.timeseries_analysis.func_VectAngleProb(ThetaInDegrees,DIMs);
    DirSIMILARITYtoFutureStateMat(1,DATANR)=1-CDFvalue;
    
end