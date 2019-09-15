function plot_timeSeries_A(analyseMeanRangeFile,analyseTppFile,analyseMakePCAskewFile,analyseCriticalRangeFile)

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

analyseMeanRangeDATA=load(analyseMeanRangeFile);
analyseTppDATA=load(analyseTppFile);
analyseMakePCAskewDATA=load(analyseMakePCAskewFile);
analyseCriticalRangeDATA=load(analyseCriticalRangeFile);

%%%%%%%%%%%%%%%%%%%%%
%%%% assign data %%%%
%%%%%%%%%%%%%%%%%%%%%

A_noiseSeries_MEAN=analyseMeanRangeDATA.A_noiseSeries_MEAN;
A_noiseSeries_LOW=analyseMeanRangeDATA.A_noiseSeries_LOW;
A_noiseSeries_HIGH=analyseMeanRangeDATA.A_noiseSeries_HIGH;
A_Xaxis_meanRange=analyseMeanRangeDATA.A_Xaxis_meanRange;

Mtpp=analyseTppDATA.Mtpp;

windowSizeFrac=analyseMakePCAskewDATA.windowSizeFrac;

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
figure(1), clf('reset')
code.figs.func_plotSeries(A_noiseSeries_MEAN,A_noiseSeries_LOW,A_noiseSeries_HIGH,A_Xaxis_meanRange,windowSizeFrac,A_CritRange_Mmin,Mtpp,A_CritRange_FOUND,COL_all)