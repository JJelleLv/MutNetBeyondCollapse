function analyseCriticalRangeFile = analyseCriticalRange(NETnr,SETnr,ExplVartrend_window,ExplVartrend_stepsize,analyseTppFile,analyseMakePCAFile,replaceAnalyseCriticalRangeFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% function that analyses timeseries and makes figures %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%%%% Output File %%%%
%%%%%%%%%%%%%%%%%%%%%

CHintMUTnumFolder = sprintf('data%sCHintMUTnumerical%s%d_INIT_NET',filesep,filesep,NETnr);
analyseCriticalRangeFile = sprintf('%s%s%d_CriticalRangeAnalysis',CHintMUTnumFolder,filesep,SETnr);

%% skip if this set already exists
if (~replaceAnalyseCriticalRangeFile)
    if exist(sprintf('%s.mat',analyseCriticalRangeFile),'file')==2
        return
    end
end

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

analyseTppDATA=load(analyseTppFile);
analyseMakePCADATA=load(analyseMakePCAFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tipping point info
TPP_ANALYSIS_DONE=analyseTppDATA.TPP_ANALYSIS_DONE;
P_DELTAab_tpp=analyseTppDATA.P_DELTAab_tpp;
A_DELTAab_tpp=analyseTppDATA.A_DELTAab_tpp;
Mtpp=analyseTppDATA.Mtpp;

%% PCA series
P_PC1_EXPL_var=analyseMakePCADATA.P_PC1_EXPL_var;
A_PC1_EXPL_var=analyseMakePCADATA.A_PC1_EXPL_var;
P_PC1_DirCHmat=analyseMakePCADATA.P_PC1_DirCHmat;
A_PC1_DirCHmat=analyseMakePCADATA.A_PC1_DirCHmat;
P_Xaxis_makePCA=analyseMakePCADATA.P_Xaxis_makePCA;
A_Xaxis_makePCA=analyseMakePCADATA.A_Xaxis_makePCA;
P_DirSIMILARITYtoFutureStateMat=analyseMakePCADATA.P_DirSIMILARITYtoFutureStateMat;
A_DirSIMILARITYtoFutureStateMat=analyseMakePCADATA.A_DirSIMILARITYtoFutureStateMat;

%% pValue threshold
TRSH_pValue=0.05;

%% dMtpp
dMtpp=0.01;

%% similarity threshold
similarityTRSH=0.99;

%%%%%%%%%%%%%%%%%%%%
%%%% empty data %%%%
%%%%%%%%%%%%%%%%%%%%

CritRange_Mmax=[];

P_CritRange_FOUND=0;
P_CritRange_Mmin=[];
P_STARTstep_CritRange=[];
P_ExplVar_CritRange=[];
P_kendTau_ExplVar_CritRange=[];
P_pValue_ExplVar_CritRange=[];

A_CritRange_FOUND=0;
A_CritRange_Mmin=[];
A_STARTstep_CritRange=[];
A_ExplVar_CritRange=[];
A_kendTau_ExplVar_CritRange=[];
A_pValue_ExplVar_CritRange=[];

%%%%%%%%%%%%%%%%%%
%%%% Analysis %%%%
%%%%%%%%%%%%%%%%%%

CRITRANGE_ANALYSIS_DONE=0;
if TPP_ANALYSIS_DONE==1
    
    CRITRANGE_ANALYSIS_DONE=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Determine Critical Range %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CritRange_Mmax=Mtpp-dMtpp;
    
    %% steps of PCA ouput to analyse
    NRsteps_makePCA=length(P_PC1_EXPL_var);
    ENDstep_CritRange=floor(CritRange_Mmax.*(NRsteps_makePCA));
    if isnan(P_PC1_EXPL_var(ENDstep_CritRange,1))==1;
        fprintf('P END STEP IS NaN!')
        ENDstep_CritRange=max(find(isnan(P_PC1_EXPL_var)==0));
    end
    if isnan(A_PC1_EXPL_var(ENDstep_CritRange,1))==1;
        fprintf('A END STEP IS NaN!')
        ENDstep_CritRange=max(find(isnan(A_PC1_EXPL_var)==0));
    end
    
    %% moving window
    STARTlist=[1:ExplVartrend_stepsize:(ENDstep_CritRange-ExplVartrend_window+1)]';
    ENDlist=[ExplVartrend_window:ExplVartrend_stepsize:ENDstep_CritRange]';
    
    NRsteps=length(STARTlist);

    P_kendTrend_ExplVar=nan(NRsteps,2);
    A_kendTrend_ExplVar=nan(NRsteps,2);
    P_Trend_ExplVar=zeros(NRsteps,1);
    A_Trend_ExplVar=zeros(NRsteps,1);
    for StepNR=1:NRsteps
        
        %% get the data belonging to this step
        P_ExplVarSERIES=P_PC1_EXPL_var([STARTlist(StepNR,1):ENDlist(StepNR,1)],:);
        A_ExplVarSERIES=A_PC1_EXPL_var([STARTlist(StepNR,1):ENDlist(StepNR,1)],:);
        
        %% check if there is a positive trend
        [P_kendTau_ExplVarSERIES,P_pValue_ExplVarSERIES]=corr([STARTlist(StepNR,1):ENDlist(StepNR,1)]',P_ExplVarSERIES,'type','Kendall');
        [A_kendTau_ExplVarSERIES,A_pValue_ExplVarSERIES]=corr([STARTlist(StepNR,1):ENDlist(StepNR,1)]',A_ExplVarSERIES,'type','Kendall');
        
        %% save Tau and pValue to list
        P_kendTrend_ExplVar(StepNR,1)=P_kendTau_ExplVarSERIES;
        P_kendTrend_ExplVar(StepNR,2)=P_pValue_ExplVarSERIES;
        A_kendTrend_ExplVar(StepNR,1)=A_kendTau_ExplVarSERIES;
        A_kendTrend_ExplVar(StepNR,2)=A_pValue_ExplVarSERIES;
        
        %% if Tau positive and p lower than 0.05, there is a positive trend
        if P_pValue_ExplVarSERIES<=TRSH_pValue
            if P_kendTau_ExplVarSERIES>=0
                P_Trend_ExplVar(StepNR,1)=1;
            end
            if P_kendTau_ExplVarSERIES<=0
                P_Trend_ExplVar(StepNR,1)=-1;
            end
        end
        if A_pValue_ExplVarSERIES<=TRSH_pValue
            if A_kendTau_ExplVarSERIES>=0
                A_Trend_ExplVar(StepNR,1)=1;
            end
            if A_kendTau_ExplVarSERIES<=0
                A_Trend_ExplVar(StepNR,1)=-1;
            end
        end
        
        %% if positive trend, check if positive trend continues till end of CritRange
        if P_Trend_ExplVar(StepNR,1)==1 && P_CritRange_FOUND==0
            P_CritRange_FOUND=1;

            for P_ENDstepNR_CritRange=ENDlist(StepNR,1):ENDlist(NRsteps,1) %%NRsteps - recently changed mistake!?
                
                %% check positive trend over CritRange
                P_ExplVar_CritRange=P_PC1_EXPL_var([STARTlist(StepNR,1):P_ENDstepNR_CritRange],:);
                [P_kendTau_ExplVar_CritRange,P_pValue_ExplVar_CritRange]=corr([STARTlist(StepNR,1):P_ENDstepNR_CritRange]',P_ExplVar_CritRange,'type','Kendall');
                if P_kendTau_ExplVar_CritRange<=0 || P_pValue_ExplVar_CritRange>TRSH_pValue
                    P_CritRange_FOUND=0;
                end
                
            end
            if P_CritRange_FOUND==1
                P_STARTstep_CritRange=STARTlist(StepNR,1);
                P_CritRange_Mmin=P_Xaxis_makePCA(P_STARTstep_CritRange,1);%+0.5.*(ExplVartrend_window./NRsteps_makePCA);
            else
                P_ExplVar_CritRange=[];
                P_kendTau_ExplVar_CritRange=[];
                P_pValue_ExplVar_CritRange=[];
                P_STARTstep_CritRange=[];
                P_CritRange_Mmin=[];
            end

        end
        
        if A_Trend_ExplVar(StepNR,1)==1 && A_CritRange_FOUND==0
            A_CritRange_FOUND=1;
            for A_ENDstepNR_CritRange=ENDlist(StepNR,1):ENDlist(NRsteps,1)
                
                %% check positive trend over CritRange
                A_ExplVar_CritRange=A_PC1_EXPL_var([STARTlist(StepNR,1):A_ENDstepNR_CritRange],:);
                [A_kendTau_ExplVar_CritRange,A_pValue_ExplVar_CritRange]=corr([STARTlist(StepNR,1):A_ENDstepNR_CritRange]',A_ExplVar_CritRange,'type','Kendall');
                if A_kendTau_ExplVar_CritRange<=0 || A_pValue_ExplVar_CritRange>TRSH_pValue
                    A_CritRange_FOUND=0;
                end
            
            end
            if A_CritRange_FOUND==1
                A_STARTstep_CritRange=STARTlist(StepNR,1);
                A_CritRange_Mmin=A_Xaxis_makePCA(A_STARTstep_CritRange,1);%+0.5.*(ExplVartrend_window./NRsteps_makePCA);   
            else
                A_ExplVar_CritRange=[];
                A_kendTau_ExplVar_CritRange=[];
                A_pValue_ExplVar_CritRange=[];
                A_STARTstep_CritRange=[];
                A_CritRange_Mmin=[];
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Kendall rank corr M/ExplVar - trend %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    P_Mrange=([P_STARTstep_CritRange:ENDstep_CritRange]./NRsteps_makePCA)';
    A_Mrange=([A_STARTstep_CritRange:ENDstep_CritRange]./NRsteps_makePCA)';
    
    P_PC1_EXPL_var_range=P_PC1_EXPL_var([P_STARTstep_CritRange:ENDstep_CritRange],1);
    A_PC1_EXPL_var_range=A_PC1_EXPL_var([A_STARTstep_CritRange:ENDstep_CritRange],1);
    
    %% data for correlation
    if P_CritRange_FOUND==1
        [P_kendTau_MExplVar,P_pValue_MExplVar]=corr(P_Mrange,P_PC1_EXPL_var_range,'type','Kendall');
    end
    
    %% data for correlation
    if A_CritRange_FOUND==1
        [A_kendTau_MExplVar,A_pValue_MExplVar]=corr(A_Mrange,A_PC1_EXPL_var_range,'type','Kendall');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Kendall rank corr M/DirCH - trend %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ENDstep_DirCHmat_CritRange=ENDstep_CritRange-1;    
    if isnan(P_PC1_DirCHmat(1,ENDstep_DirCHmat_CritRange))==1
        ENDstep_DirCHmat_CritRange=max(find(isnan(P_PC1_DirCHmat)==0));
    end
    if isnan(A_PC1_DirCHmat(1,ENDstep_DirCHmat_CritRange))==1
        ENDstep_DirCHmat_CritRange=max(find(isnan(A_PC1_DirCHmat)==0));
    end
    
    P_PC1_DirCHmat_range=P_PC1_DirCHmat(1,[P_STARTstep_CritRange:ENDstep_DirCHmat_CritRange])';
    A_PC1_DirCHmat_range=A_PC1_DirCHmat(1,[A_STARTstep_CritRange:ENDstep_DirCHmat_CritRange])';
    
    P_Mrange_minus1=([P_STARTstep_CritRange:ENDstep_DirCHmat_CritRange]./NRsteps_makePCA)';
    A_Mrange_minus1=([A_STARTstep_CritRange:ENDstep_DirCHmat_CritRange]./NRsteps_makePCA)';

    %% data for correlation
    if P_CritRange_FOUND==1
        [P_kendTau_MDirCH,P_pValue_MDirCH]=corr(P_Mrange_minus1,P_PC1_DirCHmat_range,'type','Kendall');
    end
    
    %% data for correlation
    if A_CritRange_FOUND==1
        [A_kendTau_MDirCH,A_pValue_MDirCH]=corr(A_Mrange_minus1,A_PC1_DirCHmat_range,'type','Kendall');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Similarity in critical range %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if P_CritRange_FOUND==1
        P_FRAC_Similar_CritRange=sum(P_DirSIMILARITYtoFutureStateMat(:,[P_STARTstep_CritRange:ENDstep_CritRange])>similarityTRSH)/(ENDstep_CritRange-P_STARTstep_CritRange+1);
    end
    
    if A_CritRange_FOUND==1
        A_FRAC_Similar_CritRange=sum(A_DirSIMILARITYtoFutureStateMat(:,[A_STARTstep_CritRange:ENDstep_CritRange])>similarityTRSH)/(ENDstep_CritRange-A_STARTstep_CritRange+1);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%
%%%% Store output %%%%
%%%%%%%%%%%%%%%%%%%%%%

save(analyseCriticalRangeFile,'P_CritRange_Mmin','P_CritRange_FOUND','A_CritRange_Mmin','CritRange_Mmax', ...
    'A_CritRange_FOUND','CRITRANGE_ANALYSIS_DONE')

% save(analyseCriticalRangeFile,'CritRange_Mmax','P_CritRange_FOUND','A_CritRange_FOUND','P_PC1_noiseseries_CritRange','P_PC2_noiseseries_CritRange','A_PC1_noiseseries_CritRange','A_PC2_noiseseries_CritRange', ...
%     'P_CritRange_FOUND','P_kendTau_ExplVar_CritRange','P_pValue_ExplVar_CritRange','P_STARTstep_CritRange','P_CritRange_Mmin', ...
%     'A_CritRange_FOUND','A_kendTau_ExplVar_CritRange','A_pValue_ExplVar_CritRange','A_STARTstep_CritRange','A_CritRange_Mmin', ...
%     'P_VEC_cov','P_EIG_cov','P_PC1_NONEXT_vect','P_PC1_ALL_vect','P_PC1_EXPL_var','P_PC2_NONEXT_vect','P_PC2_ALL_vect','P_PC2_EXPL_var','P_PC1_skewness','P_PC2_skewness','P_PC1_PEAKSrate','P_PC2_PEAKSrate','P_NONEXT_SpecNRs','P_NONCOLLAPSEDstate','P_PC1_vect_eucl', ...
%     'A_VEC_cov','A_EIG_cov','A_PC1_NONEXT_vect','A_PC1_ALL_vect','A_PC1_EXPL_var','A_PC2_NONEXT_vect','A_PC2_ALL_vect','A_PC2_EXPL_var','A_PC1_skewness','A_PC2_skewness','A_PC1_PEAKSrate','A_PC2_PEAKSrate','A_NONEXT_SpecNRs','A_NONCOLLAPSEDstate','A_PC1_vect_eucl', ...
%     'P_mdlOrigin_PCACh','P_regParOrigin_PCACh','P_regParPvalsOrigin_PCACh','P_regParConfOrigin_PCACh','P_regRsqOrigin_PCACh','P_regResidOrigin_PCACh','P_regPvalOrigin_PCACh', ...
%     'P_mdlLin_PCACh','P_regParLin_PCACh','P_regParPvalsLin_PCACh','P_regParConfLin_PCACh','P_regRsqLin_PCACh','P_regResidLin_PCACh','P_regPvalLin_PCACh', ...
%     'P_mdlQuadr_PCACh','P_regParQuadr_PCACh','P_regParPvalsQuadr_PCACh','P_regParConfQuadr_PCACh','P_regRsqQuadr_PCACh','P_regResidQuadr_PCACh','P_regPvalQuadr_PCACh', ...
%     'A_mdlOrigin_PCACh','A_regParOrigin_PCACh','A_regParPvalsOrigin_PCACh','A_regParConfOrigin_PCACh','A_regRsqOrigin_PCACh','A_regResidOrigin_PCACh','A_regPvalOrigin_PCACh', ...
%     'A_mdlLin_PCACh','A_regParLin_PCACh','A_regParPvalsLin_PCACh','A_regParConfLin_PCACh','A_regRsqLin_PCACh','A_regResidLin_PCACh','A_regPvalLin_PCACh', ...
%     'A_mdlQuadr_PCACh','A_regParQuadr_PCACh','A_regParPvalsQuadr_PCACh','A_regParConfQuadr_PCACh','A_regRsqQuadr_PCACh','A_regResidQuadr_PCACh','A_regPvalQuadr_PCACh', ...
%     'P_kendTau_PCACh','P_kendPval_PCACh','A_kendTau_PCACh','A_kendPval_PCACh','P_pearrho_PCACh','P_pearPval_PCACh','A_pearrho_PCACh','A_pearPval_PCACh', ...
%     'P_Mrange','P_PC1_EXPL_var_range','P_PC1_DirCHmat_range','A_Mrange','A_PC1_EXPL_var_range','A_PC1_DirCHmat_range', ...
%     'P_kendTau_MExplVar','P_pValue_MExplVar','A_kendTau_MExplVar','A_pValue_MExplVar', ...
%     'P_kendTau_MDirCH','P_pValue_MDirCH','A_kendTau_MDirCH','A_pValue_MDirCH', 'P_FRAC_Similar_CritRange', 'A_FRAC_Similar_CritRange', ...
%     'CRITRANGE_ANALYSIS_DONE');

fprintf('> Critical Range Analysis\n');
