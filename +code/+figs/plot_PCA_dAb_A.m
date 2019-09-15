function plot_PCA_dAb_A(analyseTppFile)

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

analyseTppDATA=load(analyseTppFile);

%%%%%%%%%%%%%%%%%%%%%
%%%% assign data %%%%
%%%%%%%%%%%%%%%%%%%%%

beforeTPP_PC1_ALL_vect_A_eucl=analyseTppDATA.beforeTPP_PC1_ALL_vect_A_eucl;
A_DELTAab_tpp=analyseTppDATA.A_DELTAab_tpp;
AnNRs_collapsed=analyseTppDATA.AnNRs_collapsed;
regParOrigin_PCACh_A=analyseTppDATA.regParOrigin_PCACh_A;

TPP_ANALYSIS_DONE=analyseTppDATA.TPP_ANALYSIS_DONE;

%%%%%%%%%%%%%%%%%%%%%%
%%%% Make figures %%%%
%%%%%%%%%%%%%%%%%%%%%%

if TPP_ANALYSIS_DONE==1
    
    % main colormap (10spec)
    COL_all=[0 0.4 1; ...
        1 0.4 0; ...
        0.2 0.8 0.2; ...
        0.8 0.2 0.8; ...
        0.1 0.7 0.9; ...
        0.9 0.7 0.1; ...
        0.2 0.4 0.2; ...
        0.6 0.1 0.2; ...
        0.2 0.1 0.6; ...
        0.1 0.1 0.1];
    
    %% plot timeseries
    figure(2), clf('reset')
    code.figs.func_scatterpatch_PCA_dAb(beforeTPP_PC1_ALL_vect_A_eucl,A_DELTAab_tpp,AnNRs_collapsed,regParOrigin_PCACh_A,COL_all);
    
end