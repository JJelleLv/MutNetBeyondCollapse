function plot_ExplVar_Sim_A(analyseTppFile,analyseMakePCAskewFile,analyseCriticalRangeFile)

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
TPP_ANALYSIS_DONE=analyseTppDATA.TPP_ANALYSIS_DONE;

A_PC1_EXPL_var=analyseMakePCAskewDATA.A_PC1_EXPL_var;
A_NONCOLLAPSEDstate=analyseMakePCAskewDATA.A_NONCOLLAPSEDstate;
A_DirSIMILARITYtoFutureStateMat=analyseMakePCAskewDATA.A_DirSIMILARITYtoFutureStateMat;
A_Xaxis_makePCA=analyseMakePCAskewDATA.A_Xaxis_makePCA;

A_CritRange_Mmin=analyseCriticalRangeDATA.A_CritRange_Mmin;
A_CritRange_FOUND=analyseCriticalRangeDATA.A_CritRange_FOUND;

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
    figure(4), clf('reset')
    code.figs.func_fig_scat_polyfit_2yaxis((A_PC1_EXPL_var.*A_NONCOLLAPSEDstate)',A_DirSIMILARITYtoFutureStateMat,A_Xaxis_makePCA,A_Xaxis_makePCA,A_CritRange_Mmin,Mtpp,A_CritRange_FOUND,[0 0 0],[0.9 0 0])
    
end
