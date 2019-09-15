function analyseTppFile = analyseTPP(NETnr,SETnr,analyseTPP_window_dAb,analyseTPP_window_PCA,analyseTPP_stepsize,CHintMUTnumericalFile,replaceAnalyseTppFile)

%%%%%%%%%%%%%%%%%%%%%
%%%% Output File %%%%
%%%%%%%%%%%%%%%%%%%%%

CHintMUTnumFolder = sprintf('data%sCHintMUTnumerical%s%d_INIT_NET',filesep,filesep,NETnr);
analyseTppFile = sprintf('%s%s%d_TppAnalysis',CHintMUTnumFolder,filesep,SETnr);

%% skip if this set already exists
if (~replaceAnalyseTppFile)
    if exist(sprintf('%s.mat',analyseTppFile),'file')==2
        return
    end
end

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

%% load numerical run
CHintMUTnumericalDATA = load(CHintMUTnumericalFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ptime=CHintMUTnumericalDATA.Ptime;
Atime=CHintMUTnumericalDATA.Atime;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Empty datasets %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

beforeTPP_PC1_ALL_vect_P=[];
beforeTPP_PC1_ALL_vect_P_eucl=[];
beforeTPP_PC1_skewness_P=[];
beforeTPP_VEC_cov_P=[];
beforeTPP_EIG_cov_P=[];

P_DELTAab_tpp=[];
PlNRs_collapsed=[];
regParOrigin_PCACh_P=[];

beforeTPP_PC1_ALL_vect_A=[];
beforeTPP_PC1_ALL_vect_A_eucl=[];
beforeTPP_PC1_skewness_A=[];
beforeTPP_VEC_cov_P=[];
beforeTPP_EIG_cov_P=[];

A_DELTAab_tpp=[];
AnNRs_collapsed=[];
regParOrigin_PCACh_A=[];

%%%%%%%%%%%%%%%%%%
%%%% Analysis %%%%
%%%%%%%%%%%%%%%%%%

%% make start and end position of moving window
DATAlength=length(Ptime(:,1));

STARTlist=[1:analyseTPP_stepsize:(DATAlength-analyseTPP_window_dAb+1)]';
ENDlist=[analyseTPP_window_dAb:analyseTPP_stepsize:DATAlength]';

NRsteps=length(STARTlist);

%% special start list for PCA analysis
STARTlist_PCA=ENDlist-analyseTPP_window_PCA+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIND TIPPINGPOINT %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% euclidean distance above which it is considered a tipping point
TRSH_EUCL_dist=0.5;

%% delta StepNR (to determine difference in Ab)
dSteps_TPP=3;

%% set NaN or 0
TPP_FOUND=0;
TPP_StepNR=NaN;
Mtpp=NaN;
EUCL_dists=nan(1,NRsteps);

%% check if not collapsed at step 1
PlNRs_collapsed_step1=find(sum(Ptime([STARTlist(1,1):ENDlist(1,1)],:)<=0.01)>=1)';
AnNRs_collapsed_step1=find(sum(Atime([STARTlist(1,1):ENDlist(1,1)],:)<=0.01)>=1)';

NRPl_collapsed_step1=length(PlNRs_collapsed_step1);
NRAn_collapsed_step1=length(AnNRs_collapsed_step1);

if NRPl_collapsed_step1==0 && NRAn_collapsed_step1==0
    for StepNR=2:(NRsteps-1)
        
        %% get the data belonging to the previous step
        StepSERIES_P_prev=Ptime([STARTlist((StepNR-1),1):ENDlist((StepNR-1),1)],:);
        StepSERIES_A_prev=Atime([STARTlist((StepNR-1),1):ENDlist((StepNR-1),1)],:);
        
        %% get the data belonging to the next step
        StepSERIES_P_next=Ptime([STARTlist((StepNR+1),1):ENDlist((StepNR+1),1)],:);
        StepSERIES_A_next=Atime([STARTlist((StepNR+1),1):ENDlist((StepNR+1),1)],:);
        
        %% mean of series
        meanStepSERIES_P_prev=mean(StepSERIES_P_prev);
        meanStepSERIES_A_prev=mean(StepSERIES_A_prev);
        
        meanStepSERIES_P_next=mean(StepSERIES_P_next);
        meanStepSERIES_A_next=mean(StepSERIES_A_next);
        
        %% euclidean distance
        EUCL_dists(1,StepNR)=(sum((meanStepSERIES_P_prev-meanStepSERIES_P_next).^2)+sum((meanStepSERIES_A_prev-meanStepSERIES_A_next).^2)).^0.5;
        
        %% number below extinction trsh
        NRspecNEAREXT=sum(meanStepSERIES_P_next<0.1)+sum(meanStepSERIES_A_next<0.1);
        
        %% determine if this is a TPP
        if TPP_FOUND==0 && EUCL_dists(1,StepNR)>=TRSH_EUCL_dist && NRspecNEAREXT>=1
            TPP_FOUND=1;
            TPP_StepNR=StepNR;
            Mtpp=TPP_StepNR./NRsteps;
        end
        
    end
else
    disp('Extinctions at step 1!!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Tipping point analysis %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beforeTPP_StepNR=TPP_StepNR-dSteps_TPP;
afterTPP_StepNR=TPP_StepNR+dSteps_TPP;

TPP_ANALYSIS_DONE=0;
if TPP_FOUND==1 && beforeTPP_StepNR>=1 && afterTPP_StepNR<=NRsteps
    
    TPP_ANALYSIS_DONE=1;
    
    %% timeseries before TPP - for dAb
    beforeTPP_StepSERIES_P=Ptime([STARTlist(beforeTPP_StepNR,1):ENDlist(beforeTPP_StepNR,1)],:);
    beforeTPP_StepSERIES_A=Atime([STARTlist(beforeTPP_StepNR,1):ENDlist(beforeTPP_StepNR,1)],:);
    
    %% timeseries before TPP - long, for PCA
    if STARTlist_PCA(beforeTPP_StepNR,1)>=1
        beforeTPP_StepSERIES_PCA_P=Ptime([STARTlist_PCA(beforeTPP_StepNR,1):ENDlist(beforeTPP_StepNR,1)],:);
        beforeTPP_StepSERIES_PCA_A=Atime([STARTlist_PCA(beforeTPP_StepNR,1):ENDlist(beforeTPP_StepNR,1)],:);
    else
        beforeTPP_StepSERIES_PCA_P=Ptime([1:ENDlist(beforeTPP_StepNR,1)],:);
        beforeTPP_StepSERIES_PCA_A=Atime([1:ENDlist(beforeTPP_StepNR,1)],:);
        analyseTPP_window_PCA=ENDlist(beforeTPP_StepNR,1);
    end
    
    %% timeseries after TPP
    afterTPP_StepSERIES_P=Ptime([STARTlist(afterTPP_StepNR,1):ENDlist(afterTPP_StepNR,1)],:);
    afterTPP_StepSERIES_A=Atime([STARTlist(afterTPP_StepNR,1):ENDlist(afterTPP_StepNR,1)],:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% CHANGE IN ABUNDANCE %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% collapsed species after TPP
    PlNRs_collapsed=find(sum(afterTPP_StepSERIES_P<=0.01)>=1)';
    AnNRs_collapsed=find(sum(afterTPP_StepSERIES_A<=0.01)>=1)';
    
    NRPl_collapsed=length(PlNRs_collapsed);
    NRAn_collapsed=length(AnNRs_collapsed);
    
    %% mean abundances before and after
    meanBeforeTPP_P=mean(beforeTPP_StepSERIES_P)';
    meanBeforeTPP_A=mean(beforeTPP_StepSERIES_A)';
    
    meanAfterTPP_P=mean(afterTPP_StepSERIES_P)';
    meanAfterTPP_A=mean(afterTPP_StepSERIES_A)';
    
    %% change in abundance - before after
    P_DELTAab_tpp=meanAfterTPP_P-meanBeforeTPP_P;
    A_DELTAab_tpp=meanAfterTPP_A-meanBeforeTPP_A;
    
    %% check no extinct species
    if (sum(meanBeforeTPP_P<=0.5))>=1 || (sum(meanBeforeTPP_A<=0.5))>=1
        error('extinctions before tipping point')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PCA BEFORE COLLAPSE %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %% PCA before TPP
    [beforeTPP_VEC_cov_P,beforeTPP_EIG_cov_P,beforeTPP_PC1_ALL_vect_P,beforeTPP_PC1_EXPL_var_P,beforeTPP_PC1_skewness_P,~]=code.timeseries_analysis.func_makePCA(beforeTPP_StepSERIES_PCA_P,analyseTPP_window_PCA,analyseTPP_window_PCA);
    
    beforeTPP_VEC_cov_P=beforeTPP_VEC_cov_P{1};
    beforeTPP_EIG_cov_P=beforeTPP_EIG_cov_P{1};
    beforeTPP_PC1_ALL_vect_P=beforeTPP_PC1_ALL_vect_P{1};

    [beforeTPP_VEC_cov_A,beforeTPP_EIG_cov_A,beforeTPP_PC1_ALL_vect_A,beforeTPP_PC1_EXPL_var_A,beforeTPP_PC1_skewness_A,~]=code.timeseries_analysis.func_makePCA(beforeTPP_StepSERIES_PCA_A,analyseTPP_window_PCA,analyseTPP_window_PCA);
    
    beforeTPP_VEC_cov_A=beforeTPP_VEC_cov_A{1};
    beforeTPP_EIG_cov_A=beforeTPP_EIG_cov_A{1};
    beforeTPP_PC1_ALL_vect_A=beforeTPP_PC1_ALL_vect_A{1};
    
    %% convert such that length of PC1 is equal to explained variance - and put NaN when there are collapsed species
    beforeTPP_EUCL_PC1_EXPL_var_P=beforeTPP_PC1_EXPL_var_P./(sum(beforeTPP_PC1_ALL_vect_P.^2).^0.5);
    beforeTPP_PC1_ALL_vect_P_eucl=beforeTPP_PC1_ALL_vect_P.*beforeTPP_EUCL_PC1_EXPL_var_P;
    
    beforeTPP_EUCL_PC1_EXPL_var_A=beforeTPP_PC1_EXPL_var_A./(sum(beforeTPP_PC1_ALL_vect_A.^2).^0.5);
    beforeTPP_PC1_ALL_vect_A_eucl=beforeTPP_PC1_ALL_vect_A.*beforeTPP_EUCL_PC1_EXPL_var_A;

    %% data in the direction of PC1
    beforeTPP_PC1dirSERIES_P=(beforeTPP_PC1_ALL_vect_P'*beforeTPP_StepSERIES_P');
    beforeTPP_PC1dirSERIES_P=beforeTPP_PC1dirSERIES_P-mean(beforeTPP_PC1dirSERIES_P);
    beforeTPP_PC1dirSERIES_A=(beforeTPP_PC1_ALL_vect_A'*beforeTPP_StepSERIES_A');
    beforeTPP_PC1dirSERIES_A=beforeTPP_PC1dirSERIES_A-mean(beforeTPP_PC1dirSERIES_A);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% LINEAR REGRESSION PCA/CHAb %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% regression through origin
    [~,regParOrigin_PCACh_P,~,~,~,~,~] = code.timeseries_analysis.func_regression_parOrigin(beforeTPP_PC1_ALL_vect_P_eucl,P_DELTAab_tpp);
    [~,regParOrigin_PCACh_A,~,~,~,~,~] = code.timeseries_analysis.func_regression_parOrigin(beforeTPP_PC1_ALL_vect_A_eucl,A_DELTAab_tpp);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% set direction of PC1 such that regresion is positive %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if regParOrigin_PCACh_P<0
        regParOrigin_PCACh_P=-regParOrigin_PCACh_P;

        beforeTPP_PC1_ALL_vect_P=-beforeTPP_PC1_ALL_vect_P;
        beforeTPP_PC1_ALL_vect_P_eucl=-beforeTPP_PC1_ALL_vect_P_eucl;
        beforeTPP_PC1_skewness_P=-beforeTPP_PC1_skewness_P;
        beforeTPP_PC1dirSERIES_P=-beforeTPP_PC1dirSERIES_P;
    end

    if regParOrigin_PCACh_A<0
        regParOrigin_PCACh_A=-regParOrigin_PCACh_A;
        
        beforeTPP_PC1_ALL_vect_A=-beforeTPP_PC1_ALL_vect_A;
        beforeTPP_PC1_ALL_vect_A_eucl=-beforeTPP_PC1_ALL_vect_A_eucl;
        beforeTPP_PC1_skewness_A=-beforeTPP_PC1_skewness_A;
        beforeTPP_PC1dirSERIES_A=-beforeTPP_PC1dirSERIES_A;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Similarity PCA/CHAb %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [DirCHtoFutureState_PCACh_P,DirSIMILARITYtoFutureState_PCACh_P]=code.timeseries_analysis.func_CHdirToFutureState(beforeTPP_PC1_ALL_vect_P,P_DELTAab_tpp);
    [DirCHtoFutureState_PCACh_A,DirSIMILARITYtoFutureState_PCACh_A]=code.timeseries_analysis.func_CHdirToFutureState(beforeTPP_PC1_ALL_vect_A,A_DELTAab_tpp);
        
end

%%%%%%%%%%%%%%%%%%%%%%
%%%% Store output %%%%
%%%%%%%%%%%%%%%%%%%%%%

save(analyseTppFile,'analyseTPP_window_dAb','analyseTPP_window_PCA','analyseTPP_stepsize', ...
    'beforeTPP_PC1_ALL_vect_P','beforeTPP_PC1_ALL_vect_P_eucl','beforeTPP_PC1_skewness_P', ...
    'beforeTPP_VEC_cov_P','beforeTPP_EIG_cov_P','P_DELTAab_tpp','PlNRs_collapsed','regParOrigin_PCACh_P', ...
    'beforeTPP_PC1_ALL_vect_A','beforeTPP_PC1_ALL_vect_A_eucl','beforeTPP_PC1_skewness_A', ...
    'beforeTPP_VEC_cov_P','beforeTPP_EIG_cov_P','A_DELTAab_tpp','AnNRs_collapsed','regParOrigin_PCACh_A', ...
    'Mtpp','TPP_ANALYSIS_DONE')

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Display output %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('> TPP Analysis\n');
fprintf('TPP found at: %.2f\n', Mtpp);