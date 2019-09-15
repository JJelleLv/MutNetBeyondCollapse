function MAIN_MutNetBeyondCollapse(NETnr,SETnr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% studying the direction of critical slowing down %%%%%%%%%%%%%%%%%%%%
%%%% prior to critical transitions in bipartite mutualistic networks %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% close all
close all

%% randomize
rng('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1: Make time series (using GRIND for MATLAB) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load time series (false) or make new simulation (true)
replaceCHintMUTnumericalFile = true;
CHintMUTnumericalFile = code.gradual_change.CHintMUTnumerical(NETnr,SETnr,replaceCHintMUTnumericalFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2: Get Mean and CI - for timeseriesplot %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
replaceAnalysemeanRangeFile = true;
meanRange_window=50;
meanRange_stepsize=50;

%% Analysis
analyseMeanRangeFile = code.timeseries_analysis.meanRange(NETnr,SETnr,meanRange_window,meanRange_stepsize,CHintMUTnumericalFile,replaceAnalysemeanRangeFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 3: Analyze Tipping Point %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
replaceAnalyseTppFile = true;
analyseTPP_window_dAb=200;
analyseTPP_window_PCA=2000;
analyseTPP_stepsize=200;

%% Analysis
analyseTppFile = code.timeseries_analysis.analyseTPP(NETnr,SETnr,analyseTPP_window_dAb,analyseTPP_window_PCA,analyseTPP_stepsize,CHintMUTnumericalFile,replaceAnalyseTppFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 4: Analyze Time Series %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
replaceAnalyseMakePCAskewFile = true;
makePCAskew_window=2000;
makePCAskew_stepsize=200;

%% Analysis
analyseMakePCAskewFile = code.timeseries_analysis.makePCA(NETnr,SETnr,makePCAskew_window,makePCAskew_stepsize,CHintMUTnumericalFile,analyseTppFile,replaceAnalyseMakePCAskewFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 5: Determine Critical Range - trend in Expl. Var. %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings Analysis 4: Critical range
replaceAnalyseCriticalRangeFile = true;
ExplVartrend_window=10;
ExplVartrend_stepsize=1;

%% Analysis 4: determine and analyse critical range
analyseCriticalRangeFile = code.timeseries_analysis.analyseCriticalRange(NETnr,SETnr,ExplVartrend_window,ExplVartrend_stepsize,analyseTppFile,analyseMakePCAskewFile,replaceAnalyseCriticalRangeFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 6: plot figures %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot time series
code.figs.plot_timeSeries_A(analyseMeanRangeFile,analyseTppFile,analyseMakePCAskewFile,analyseCriticalRangeFile);

%% plot dAb and PC1 at/before Tipping point
code.figs.plot_PCA_dAb_A(analyseTppFile);

%% plot PCA
code.figs.plot_PCA_A(analyseTppFile,analyseMakePCAskewFile,analyseCriticalRangeFile);

%% plot Explained Variance and Similarity
code.figs.plot_ExplVar_Sim_A(analyseTppFile,analyseMakePCAskewFile,analyseCriticalRangeFile);

