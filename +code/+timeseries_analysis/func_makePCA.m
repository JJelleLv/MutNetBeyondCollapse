function [VEC_cov,EIG_cov,PC1_ALL_vect,PC1_EXPL_var,PC1_skewness,NONCOLLAPSEDstate,Xaxis_makePCA]=func_makePCA(noiseseries,window,stepsize)

%% make start and end position of moving window
DATAlength=length(noiseseries(:,1));
NRdata=length(noiseseries(1,:));

STARTlist=[1:stepsize:(DATAlength-window+1)]';
ENDlist=[window:stepsize:DATAlength]';

NRsteps=length(STARTlist);

%% empty datasets
COV=cell(1,NRsteps);
VEC_cov=cell(1,NRsteps);
EIG_cov=cell(1,NRsteps);

PC1_NONEXT_vect=cell(1,NRsteps);
PC1_EXPL_var=nan(NRsteps,1);

PC1_ALL_vect=cell(1,NRsteps);

PC1_skewness=nan(NRsteps,1);

NONEXT_SpecNRs=cell(1,NRsteps);
NONCOLLAPSEDstate=nan(NRsteps,1);

%% PCA analysis
for StepNR=1:NRsteps
    
    %% empty datasets
    COV{StepNR}=nan(NRdata,NRdata);
    PC1_ALL_vect{StepNR}=nan(NRdata,1);
    
    %% get the data belonging to this step
    StepSERIES=noiseseries([STARTlist(StepNR,1):ENDlist(StepNR,1)],:);
    
    %% non-extinct species
    NONEXT_SpecNRs{StepNR}=find(sum(StepSERIES<=0.01)==0);
    NONEXT_StepSERIES=StepSERIES(:,NONEXT_SpecNRs{StepNR});
    
    %% save when non-collapsed
    if length(NONEXT_SpecNRs{StepNR})==NRdata
        NONCOLLAPSEDstate(StepNR,1)=1;
    end
    
    %% if at least one species non-extinct
    if length(NONEXT_SpecNRs{StepNR})>=1
        
        %% PCA over last steps
        COV{StepNR}=cov(NONEXT_StepSERIES);
        [VEC_cov{StepNR},EIG_cov{StepNR}]=eig(COV{StepNR});
        
        %% store info of PC1
        DOM_EIG_NR=find((diag(EIG_cov{StepNR})==max(max(real(EIG_cov{StepNR})))));
        PC1_NONEXT_vect{StepNR}=VEC_cov{StepNR}(:,DOM_EIG_NR);
        PC1_EXPL_var(StepNR,1)=(EIG_cov{StepNR}(DOM_EIG_NR,DOM_EIG_NR))/sum(sum((real(EIG_cov{StepNR}))));
        
        %% add PC vectors with nan for extinct species
        PC1_ALL_vect{StepNR}(NONEXT_SpecNRs{StepNR},:)=PC1_NONEXT_vect{StepNR};

        %% data in direcetion of PC1
        PC1_NONEXT_StepSERIES=(PC1_NONEXT_vect{StepNR}'*NONEXT_StepSERIES');
        PC1_NONEXT_StepSERIES=PC1_NONEXT_StepSERIES-mean(PC1_NONEXT_StepSERIES);
        
        %% determine skewness - in direction of PC1
        PC1_skewness(StepNR,1)=skewness(PC1_NONEXT_StepSERIES);

    end
    
end

%% Xaxis_meanrange - used for plotting
%Xaxis_makePCA=([0.5*window:stepsize:DATAlength-0.5*window+1]')./DATAlength; %% point in the middle of window
Xaxis_makePCA=([window:stepsize:DATAlength+1]')./DATAlength; %% point at the end of window