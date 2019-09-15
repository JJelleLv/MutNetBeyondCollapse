function analyseMeanRangeFile = meanRange(NETnr,SETnr,meanRange_window,meanRange_stepsize,CHintMUTnumericalFile,replaceanalyseMeanRangeFile)

%%%%%%%%%%%%%%%%%%%%%
%%%% Output File %%%%
%%%%%%%%%%%%%%%%%%%%%

CHintMUTnumFolder = sprintf('data%sCHintMUTnumerical%s%d_INIT_NET',filesep,filesep,NETnr);
analyseMeanRangeFile = sprintf('%s%s%d_meanRangeAnalysis',CHintMUTnumFolder,filesep,SETnr);

%% skip if this set already exists
if (~replaceanalyseMeanRangeFile)
    if exist(sprintf('%s.mat',analyseMeanRangeFile),'file')==2
        return
    end
end

%%%%%%%%%%%%%%%%%%%
%%%% load data %%%%
%%%%%%%%%%%%%%%%%%%

%% load numerical run
CHintMUTnumericalDATA=load(CHintMUTnumericalFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ptime=CHintMUTnumericalDATA.Ptime;
Atime=CHintMUTnumericalDATA.Atime;

%%%%%%%%%%%%%%%%%%
%%%% Analysis %%%%
%%%%%%%%%%%%%%%%%%

%% determines upper and lower boundary for CI
RANGE=0.9;

%% mean and confidence interval
[P_noiseSeries_MEAN,P_noiseSeries_LOW,P_noiseSeries_HIGH,P_Xaxis_meanRange]=code.timeseries_analysis.func_meanRange(Ptime,RANGE,meanRange_window,meanRange_stepsize);
[A_noiseSeries_MEAN,A_noiseSeries_LOW,A_noiseSeries_HIGH,A_Xaxis_meanRange]=code.timeseries_analysis.func_meanRange(Atime,RANGE,meanRange_window,meanRange_stepsize);

%%%%%%%%%%%%%%%%%%%%%%
%%%% Store output %%%%
%%%%%%%%%%%%%%%%%%%%%%

save(analyseMeanRangeFile,'meanRange_window','meanRange_stepsize','RANGE', ...
    'P_noiseSeries_MEAN','P_noiseSeries_LOW','P_noiseSeries_HIGH','P_Xaxis_meanRange', ...
    'A_noiseSeries_MEAN','A_noiseSeries_LOW','A_noiseSeries_HIGH','A_Xaxis_meanRange')

fprintf('> Mean/range analysis\n');