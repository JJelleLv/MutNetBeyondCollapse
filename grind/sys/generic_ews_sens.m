%function generic_ews_sens
% make heatmaps of the sensitivity of the tau values of generic_ews of the
% window size and bandwithd
%
% Usage
%
% generic_ews_sens(timeseries, option pairs)
%
% all options pairs: (with defaults)
% res= generic_ews_sens(timeseries, 'winsize', (5:5:90), 'detrending', 'no', ...
%     'bandwidth', (1:1:20), 'logtransform', false, 'interpolate', false)
% To draw plots:
% generic_ews_sens('-p',res,nrs,unit) = unit can be '%' or something else
% generic_ews_sens('-p',res,1,'yr') = [plot the first tau value = AR1,
% translate % to the unit of the time series
%
% Arguments (only those different from sgeneric_ews)
%
% winsize	   - an array with the window sizes to be tested (in % of the
%                size of the time series)
%
% bandwidth	   - an array with the bandwidths to be tested (in % of the
%                size of the time series)
%
%
% generic_ews_sens returns a struct with the fields:
%
%     bandwidths:   matrix with all bandwidth tested
%     winsizes:     matrix with the window sizes
%     taus:         3D matix with the taus (taus(:,:,1)=AR1,
%                   taus(:,:,2)=acf, taus(:,:,3)=standard dev/CV, taus(:,:,4)=skewness)
%     timeseries:   the original time series
%     title:        title of the time series (default is empty)
%     datalength:   length of the data set (can be used to translate % to real
%                    time values
%
% Author(s)
%
% Vasilis Dakos vasilis.dakos@gmail.com
% MATLAB version by Egbert van Nes
%
% References
%
% Dakos, V., et al (2008). "Slowing down as an early warning signal for abrupt climate change." Proceedings of the National Academy of Sciences 105(38): 14308-14312
%
% Dakos, V., et al (2012)."Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data." PLoS ONE 7(7): e41010. doi:10.1371/journal.pone.0041010
%
%
function [res] = generic_ews_sens(timeseries,varargin)
if nargin>0 && ischar(timeseries) && strncmp(timeseries,'-p',2)
    bandwidths=varargin{1}.bandwidths;
    winsizes=varargin{1}.winsizes;
    taus=varargin{1}.taus;
    if nargin>2
        figs=varargin{2}';
    else
        figs=1:4;
    end
    if nargin>3
        units=varargin{3};
    else
        units='%';
    end
    if ~strcmp(units,'%')
        bandwidths=bandwidths*varargin{1}.datalength/100;
        winsizes=winsizes*varargin{1}.datalength/100;
    end
    titles={'autoregression','acf','standard deviation','skewness'};
    for i=figs
        figure
        contourf(bandwidths,winsizes,taus(:,:,i));
        i_plotdefaults;
        xlabel(sprintf('filtering bandwidth (%s)',units))
        ylabel(sprintf('sliding window (%s)',units));
        h=colorbar;
        i_plotdefaults;
        ylabel(h,titles{i});
        if ~isempty(varargin{1}.title)
            title(varargin{1}.title);
        else
            title(titles{i});
        end
    end
    return;
end
EWSopt = struct('winsize', (5:5:90), ...
    'detrending','gaussian',... %  = c("no", "gaussian", "linear", "first-diff"),
    'bandwidth', (1:1:20), ...
    'silent',false,...
    'logtransform', false, ...
    'cv',false,...
    'arlag', 1, ...
    'interpolate', false, ...
    'AR_n', false, ...
    'powerspectrum' , false,...
    'title','');

f=fieldnames(EWSopt);
for i = 1:2:length(varargin)
    if ~any(strcmp(f,varargin{i}))
        s=sprintf('"%s", ',f{:});
        error('Unknown option for generic_ews_sens: %s\nValid options are: %s',varargin{i},s(1:end-1));
    end
    EWSopt.(varargin{i}) = varargin{i + 1};
end
EWSopt.silent=true;
if EWSopt.powerspectrum||EWSopt.AR_n
    error('Option not yet supported');
end
taus=zeros(numel(EWSopt.winsize),numel(EWSopt.bandwidth),4);
[bandwidths,winsizes]=meshgrid(EWSopt.bandwidth,EWSopt.winsize);
if exist('i_waitbar','file')==0
    initgrind;
end
i_waitbar(0,size(winsizes,1),'sensitivity of generic_ews')
for i=1:size(winsizes,1)
    for j=1:size(winsizes,2)
        EWSopt.winsize=winsizes(i,j);
        EWSopt.bandwidth=bandwidths(i,j);
        [~,~,taus1,opts]=generic_ews(timeseries,EWSopt);
        taus(i,j,:)=taus1;
    end
    i_waitbar(1);
end
i_waitbar([]);
res.bandwidths=bandwidths;
res.winsizes=winsizes;
res.taus=taus;
res.timeseries=timeseries;
res.title=EWSopt.title;
res.datalength=opts.datalength;
