function plot_PCA_A(analyseTppFile,analyseMakePCAskewFile,analyseCriticalRangeFile)

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

analyseTppDATA=load(analyseTppFile);
analyseMakePCAskewDATA=load(analyseMakePCAskewFile);
analyseCriticalRangeDATA=load(analyseCriticalRangeFile);

%%%%%%%%%%%%%%%%%%%%%
%%%% assign data %%%%
%%%%%%%%%%%%%%%%%%%%%

Mtpp=analyseTppDATA.Mtpp;

A_PC1_vect_mat_eucl_switch=analyseMakePCAskewDATA.A_PC1_vect_mat_eucl_switch;
A_Xaxis_makePCA=analyseMakePCAskewDATA.A_Xaxis_makePCA;

A_CritRange_Mmin=analyseCriticalRangeDATA.A_CritRange_Mmin;
A_CritRange_FOUND=analyseCriticalRangeDATA.A_CritRange_FOUND;

%%%%%%%%%%%%%%%%%%%%%%
%%%% Make figures %%%%
%%%%%%%%%%%%%%%%%%%%%%

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
figure(3), clf('reset')
code.figs.func_fig_scat_polyfit(A_PC1_vect_mat_eucl_switch,A_Xaxis_makePCA,A_CritRange_Mmin,Mtpp,A_CritRange_FOUND,COL_all)
