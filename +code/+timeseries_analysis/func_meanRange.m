function [noiseseries_MEAN,noiseseries_LOW,noiseseries_HIGH,Xaxis_meanrange]=func_meanRange(noiseseries,RANGE,window,stepsize)

%% make start and end position of moving window
DATAlength=length(noiseseries(:,1));
NRdata=length(noiseseries(1,:));

STARTlist=[1:stepsize:(DATAlength-window+1)]'; %% lengte 90
ENDlist=[window:stepsize:DATAlength]'; %% lengte 91...

NRsteps=length(STARTlist);

%% postition sorted series, confidence interval
LOW=(1-RANGE)./2;
HIGH=1-LOW;
POS_LOW=round(LOW*window);
POS_HIGH=round(HIGH*window);

%% determine mean and range
noiseseries_MEAN=nan(NRsteps,NRdata);
noiseseries_HIGH=nan(NRsteps,NRdata);
noiseseries_LOW=nan(NRsteps,NRdata);
for StepNR=1:NRsteps
    
    %% get the data belonging to this step
    StepSERIES=noiseseries([STARTlist(StepNR,1):ENDlist(StepNR,1)],:);
    
    %% sort series
    SORT_StepSERIES=sort(StepSERIES);
    
    %% add to list
    noiseseries_MEAN(StepNR,:)=mean(StepSERIES);
    noiseseries_LOW(StepNR,:)=SORT_StepSERIES(POS_LOW,:);
    noiseseries_HIGH(StepNR,:)=SORT_StepSERIES(POS_HIGH,:);

end

%% Xaxis_meanrange - used for plotting
Xaxis_meanrange=([0.5*window:stepsize:DATAlength-0.5*window+1]')./DATAlength;