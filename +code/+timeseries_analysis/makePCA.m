function analyseMakePCAFile = makePCA(NETnr,SETnr,makePCAskew_window,makePCAskew_stepsize,CHintMUTnumericalFile,analyseTppFile,replaceAnalyseMakePCAFile)

%%%%%%%%%%%%%%%%%%%%%
%%%% Output File %%%%
%%%%%%%%%%%%%%%%%%%%%

CHintMUTnumFolder = sprintf('data%sCHintMUTnumerical%s%d_INIT_NET',filesep,filesep,NETnr);
analyseMakePCAFile = sprintf('%s%s%d_PCAAnalysis',CHintMUTnumFolder,filesep,SETnr);

%% skip if this set already exists
if (~replaceAnalyseMakePCAFile)
    if exist(sprintf('%s.mat',analyseMakePCAFile),'file')==2
        return
    end
end

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

%% load numerical run
CHintMUTnumericalDATA = load(CHintMUTnumericalFile);
analyseTppDATA=load(analyseTppFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ptime=CHintMUTnumericalDATA.Ptime;
Atime=CHintMUTnumericalDATA.Atime;

P_DELTAab_tpp=analyseTppDATA.P_DELTAab_tpp;
A_DELTAab_tpp=analyseTppDATA.A_DELTAab_tpp;
TPP_ANALYSIS_DONE=analyseTppDATA.TPP_ANALYSIS_DONE;

%%%%%%%%%%%%%%%%%%%%
%%%% Empty data %%%%
%%%%%%%%%%%%%%%%%%%%

P_PC1_DirCHtoFutureStateMat=[];
P_DirSIMILARITYtoFutureStateMat=[];
A_PC1_DirCHtoFutureStateMat=[];
A_DirSIMILARITYtoFutureStateMat=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine window size as fraction of time series %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windowSizeFrac=makePCAskew_window/(length(Ptime));

%%%%%%%%%%%%%%%%%%
%%%% Analysis %%%%
%%%%%%%%%%%%%%%%%%

%% PCA analysis
%[P_VEC_cov,P_EIG_cov,~,P_PC1_ALL_vect,P_PC1_EXPL_var,~,~,~,P_PC1_skewness,~,~,~,~,P_NONCOLLAPSEDstate,~,~,P_Xaxis_makePCA,~,~]=code.timeseries_analysis.func_makePCA(Ptime,makePCAskew_window,makePCAskew_stepsize);
%[A_VEC_cov,A_EIG_cov,~,A_PC1_ALL_vect,A_PC1_EXPL_var,~,~,~,A_PC1_skewness,~,~,~,~,A_NONCOLLAPSEDstate,~,~,A_Xaxis_makePCA,~,~]=code.timeseries_analysis.func_makePCA(Atime,makePCAskew_window,makePCAskew_stepsize);
[P_VEC_cov,P_EIG_cov,P_PC1_ALL_vect,P_PC1_EXPL_var,P_PC1_skewness,P_NONCOLLAPSEDstate,P_Xaxis_makePCA]=code.timeseries_analysis.func_makePCA(Ptime,makePCAskew_window,makePCAskew_stepsize);
[A_VEC_cov,A_EIG_cov,A_PC1_ALL_vect,A_PC1_EXPL_var,A_PC1_skewness,A_NONCOLLAPSEDstate,A_Xaxis_makePCA]=code.timeseries_analysis.func_makePCA(Atime,makePCAskew_window,makePCAskew_stepsize);

%% convert to matrix
P_PC1_ALL_vect_mat=cell2mat(P_PC1_ALL_vect);
A_PC1_ALL_vect_mat=cell2mat(A_PC1_ALL_vect);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Scale with explained variance %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert such that length of PC1 is equal to explained variance - and put NaN when there are collapsed species
P_EUCL_Expl_var=P_PC1_EXPL_var./(sum(P_PC1_ALL_vect_mat.^2).^0.5)'; %% this conversion is not needed because matlab output has EUCL distance 1
P_PC1_vect_mat_eucl=P_PC1_ALL_vect_mat.*(P_EUCL_Expl_var*ones(1,length(P_PC1_ALL_vect_mat(:,1))))';
A_EUCL_Expl_var=A_PC1_EXPL_var./(sum(A_PC1_ALL_vect_mat.^2).^0.5)';
A_PC1_vect_mat_eucl=A_PC1_ALL_vect_mat.*(A_EUCL_Expl_var*ones(1,length(A_PC1_ALL_vect_mat(:,1))))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% switch direction of PCs - per Step %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% determine switch direction based on previous PC
P_PC1_switch=code.timeseries_analysis.func_make_PC_Switch_eucl(P_PC1_ALL_vect_mat);
A_PC1_switch=code.timeseries_analysis.func_make_PC_Switch_eucl(A_PC1_ALL_vect_mat);

%% switch PC1 based on switch direction
P_PC1_vect_mat_eucl_switch=P_PC1_vect_mat_eucl.*(P_PC1_switch*ones(1,length(P_PC1_vect_mat_eucl(:,1))))';
A_PC1_vect_mat_eucl_switch=A_PC1_vect_mat_eucl.*(A_PC1_switch*ones(1,length(A_PC1_vect_mat_eucl(:,1))))';

%% switch direction of skewness based on switch direction
P_PC1_skewness_switch=P_PC1_skewness.*P_PC1_switch;
A_PC1_skewness_switch=A_PC1_skewness.*A_PC1_switch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% switch direction of PCs - full series %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TPP_ANALYSIS_DONE==1
    
    %% set in direction of Ab Change
    maxSTEPNR_NONC_P=max(find(P_NONCOLLAPSEDstate==1));
    maxSTEPNR_NONC_A=max(find(A_NONCOLLAPSEDstate==1));
    
    %% regression
    regParOrigin_PCACh_P=[];
    regParOrigin_PCACh_A=[];
    if maxSTEPNR_NONC_P>=1
        [~,regParOrigin_PCACh_P,~,~,~,~,~]=code.timeseries_analysis.func_regression_parOrigin(P_PC1_vect_mat_eucl_switch(:,maxSTEPNR_NONC_P),P_DELTAab_tpp);
    end
    if maxSTEPNR_NONC_A>=1
        [~,regParOrigin_PCACh_A,~,~,~,~,~]=code.timeseries_analysis.func_regression_parOrigin(A_PC1_vect_mat_eucl_switch(:,maxSTEPNR_NONC_A),A_DELTAab_tpp);
    end
    
    %% switch direction if negative regression coeff
    if regParOrigin_PCACh_P<0
        P_PC1_vect_mat_eucl_switch=-P_PC1_vect_mat_eucl_switch;
        P_PC1_skewness_switch=-P_PC1_skewness_switch;
    end
    
    if regParOrigin_PCACh_A<0
        A_PC1_vect_mat_eucl_switch=-A_PC1_vect_mat_eucl_switch;
        A_PC1_skewness_switch=-A_PC1_skewness_switch;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% determine change in direction %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P_PC1_DirCHmat]=code.timeseries_analysis.func_CHdirVect(P_PC1_vect_mat_eucl_switch);
[A_PC1_DirCHmat]=code.timeseries_analysis.func_CHdirVect(A_PC1_vect_mat_eucl_switch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% determine difference in direction between indicated and observed shift %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TPP_ANALYSIS_DONE==1
    [P_PC1_DirCHtoFutureStateMat,P_DirSIMILARITYtoFutureStateMat]=code.timeseries_analysis.func_CHdirToFutureState(P_PC1_vect_mat_eucl_switch,P_DELTAab_tpp);
    [A_PC1_DirCHtoFutureStateMat,A_DirSIMILARITYtoFutureStateMat]=code.timeseries_analysis.func_CHdirToFutureState(A_PC1_vect_mat_eucl_switch,A_DELTAab_tpp);
end

%%%%%%%%%%%%%%%%%%%%%%
%%%% Store output %%%%
%%%%%%%%%%%%%%%%%%%%%%

save(analyseMakePCAFile,'makePCAskew_window','makePCAskew_stepsize','windowSizeFrac', ...
    'P_PC1_vect_mat_eucl_switch','P_Xaxis_makePCA','P_PC1_EXPL_var','P_NONCOLLAPSEDstate', ...
    'P_DirSIMILARITYtoFutureStateMat','P_PC1_DirCHmat', ...
    'A_PC1_vect_mat_eucl_switch','A_Xaxis_makePCA','A_PC1_EXPL_var','A_NONCOLLAPSEDstate', ...
    'A_DirSIMILARITYtoFutureStateMat','A_PC1_DirCHmat')

fprintf('> PCA Analysis\n');



