%function generic_ews
% generic_ews is used to estimate statistical moments within rolling
% windows along a timeserie (based on the R early warning signals toolbox,
% but simplified)
%
% Usage
%
% generic_ews(timeseries, option pairs)
%
% all options pairs: (with defaults)
% generic_ews(timeseries, 'winsize', 50, 'detrending', 'no', ...
%     'bandwidth', [], 'logtransform', false, 'interpolate', false)
%
% Arguments
%
% timeseries   - a numeric vector of the observed univariate timeseries values or a numeric
%                matrix where the first column represents the time index and the second the
%                observed timeseries values. Use vectors/matrices with headings.
%
% winsize	   - is the size of the rolling window expressed as percentage of the timeseries length
%                (must be numeric between 0 and 100). Default is 50%.
%
% bandwidth	   - is the bandwidth used for the Gaussian kernel when gaussian filtering is applied.
%                It is expressed as percentage of the timeseries length (must be numeric between 0 and 100)
%                Alternatively it can be given by the optimal bandwidth suggested by Bowman and Azzalini (1997) Default).
%
% detrending   - the timeseries can be detrended/filtered prior to analysis.
%                There are four options: 'no'= no detrending, 'gaussian' = gaussian filtering,
%                'linear' = linear detrending, or 'first-diff' = first-differencing. Default is 'no' detrending.
%
% logtransform - if TRUE data are logtransformed prior to analysis as log(X+1). Default is FALSE.
%
% interpolate  - If TRUE linear interpolation is applied to produce a timeseries of equal length as the original. Default is FALSE (assumes there are no gaps in the timeseries).
% arlag        - The lag for autoregression. Default = 1
%
%not yet supported in the MATLAB version:
% AR_n	       - If TRUE the best fitted AR(n) model is fitted to the data. Default is FALSE.
%
% powerspectrum	-If TRUE the power spectrum within each rolling window is plotted. Default is FALSE.
%
%
% generic_ews returns a table (or matrix for older matlab varsions the statisics toolbox is present) that contains:
%
% tim	- the time index.
% ar1	- the autoregressive coefficient ar(1) of a first order AR model fitted on the data within the
%         rolling window.
% sd	- the standard deviation of the data estimated within each rolling window.
% sk	- the skewness of the data estimated within each rolling window.
% kurt	- the kurtosis of the data estimated within each rolling window.
% cv	- the coefficient of variation of the data estimated within each rolling window.
% returnrate - the return rate of the data estimated as 1-ar(1) cofficient within each rolling window.
% densratio	- not supported in MATLAB:  the density ratio of the power spectrum of the data estimated
%           as the ratio of low frequencies over high frequencies within each rolling window.
% acf1	- the autocorrelation at first lag of the data estimated within each rolling window.
%
% In addition, generic_ews returns a plots. The  plot contains the original data,
% the detrending/filtering applied and the residuals (if selected), and all the moment statistics.
% For each statistic trends are estimated by the nonparametric Kendall tau correlation (if statistics toolbox
% is present).
% Not supported: The second plot, if asked, quantifies resilience indicators fitting AR(n) selected by the
% Akaike Information Criterion. The third plot, if asked, is the power spectrum estimated
% by spec.ar for all frequencies within each rolling window.
%
% Author(s)
%
% Vasilis Dakos vasilis.dakos@gmail.com
% MATLAB version by Egbert van Nes
%
% References
%
% Ives, A. R. (1995). "Measuring resilience in stochastic systems." Ecological Monographs 65: 217-233
%
% Dakos, V., et al (2008). "Slowing down as an early warning signal for abrupt climate change." Proceedings of the National Academy of Sciences 105(38): 14308-14312
%
% Dakos, V., et al (2012)."Methods for Detecting Early Warnings of Critical Transitions in Time Series Illustrated Using Simulated Ecological Data." PLoS ONE 7(7): e41010. doi:10.1371/journal.pone.0041010
%
%
function [inds, trend, taus, EWSopt] = generic_ews(timeseries,varargin)
if isa(timeseries,'table')
    timeseries=timeseries{:,:};
else
    timeseries = double(timeseries);
end

EWSopt = struct('winsize', 50, ...
    'detrending','gaussian',... %  = c("no", "gaussian", "linear", "first-diff"),
    'bandwidth', 10, ...
    'silent',false,...
    'logtransform', false, ...
    'arlag', 1, ...
    'cv',false,...
    'interpolate', false, ...
    'AR_n', false, ...
    'powerspectrum' , false,...
    'title','');
if EWSopt.powerspectrum||EWSopt.AR_n
    error('Option not yet supported');
end
taus=zeros(4,1);
trend = [];
f=fieldnames(EWSopt);
if ~isempty(varargin)&&isstruct(varargin{1})
    %all options can be entered as struct
    if exist('mergestructs','file')==0
        initgrind;
    end
    %replace common fields in EWSopt
    EWSopt=mergestructs(EWSopt,varargin{1});
else
    for i = 1:2:length(varargin)
        if ~any(strcmp(f,varargin{i}))
            s=sprintf('"%s", ',f{:});
            error('Unknown option for generic_ews: %s\nValid options are: %s',varargin{i},s(1:end-1));
        end
        EWSopt.(varargin{i}) = varargin{i + 1};
    end
end
if size(timeseries, 2) < 2
    timeseries = [transpose(linspace(1,size(timeseries,1),size(timeseries,1))), timeseries];
elseif any(diff(timeseries(:,1))<0)
    error('generic_ews:timeseries','Error in generic_ews: time (first column) must be monotoneously increasing')
end
timeorig=timeseries;
if EWSopt.interpolate
    ts = transpose(linspace(timeseries(1, 1), timeseries(1, end), size(timeseries, 1)));
    Ys = interp1(timeseries(:,1),timeseries(:,2), ts);
    timeseries = [ts, Ys];
end
if EWSopt.logtransform
    timeseries(:, 2) = log(timeseries(:, 2) + 1);
end
dataorig = timeseries;
EWSopt.datalength=max(timeseries(:, 1)) - min(timeseries(:, 1));
if ~strcmpi(EWSopt.detrending, 'no')
    if strcmpi(EWSopt.detrending, 'first-diff')
        timeseries(:, 2) = [nan; diff(timeseries(:, 2))];
    elseif strcmpi(EWSopt.detrending, 'linear')
        p=polyfit(timeseries(:,1), timeseries(:,2),1);
        trend=polyval(p,timeseries(:,1));
        timeseries(:,2)=timeseries(:,2)-trend;
    else
        if ~isempty(EWSopt.bandwidth)
            EWSopt.absbandwidth=EWSopt.bandwidth/100*(max(timeseries(:, 1)) - min(timeseries(:, 1)));
        else
            EWSopt.absbandwidth=[];
        end
        trend = ksmooth(timeseries(:,1), timeseries(:,2), EWSopt.absbandwidth, EWSopt.detrending);
        timeseries(:, 2) = timeseries(:, 2) - trend;
    end
end
if EWSopt.winsize > 100
    EWSopt.winsize = 100;
end
if EWSopt.winsize < 0.1
    EWSopt.winsize = 0.1;
end
ndx=~isnan(trend);
timeseries=timeseries(ndx,:);
dataorig=dataorig(ndx,:);
trend=trend(ndx);
EWSopt.avwin = EWSopt.winsize / 100 * (max(timeseries(:, 1)) - min(timeseries(:, 1)));

ars = moving_window(timeseries(:,1), timeseries(:,2), EWSopt.avwin, @autoregression,EWSopt.arlag);
acfs =  moving_window(timeseries(:,1), timeseries(:,2), EWSopt.avwin, @acf);
if EWSopt.cv
    vars = moving_window(timeseries(:,1), timeseries(:,2), EWSopt.avwin, 'cv',dataorig(:,2:end));
else
    vars = moving_window(timeseries(:,1), timeseries(:,2), EWSopt.avwin, @std);
end
skews = moving_window(timeseries(:,1), timeseries(:,2), EWSopt.avwin, @skewness);

hasstats= exist('corr','file')~=0;
xlims=[min(dataorig(:,1)),max(dataorig(:,1))];
if hasstats
    taus(1)=corr(timeseries(:,1), ars, 'type','Kendall','rows','complete');
    taus(2)=corr(timeseries(:,1), acfs, 'type','Kendall','rows','complete');
    taus(3)=corr(timeseries(:,1), vars, 'type','Kendall','rows','complete');
    taus(4)=corr(timeseries(:,1), skews, 'type','Kendall','rows','complete');
end
if ~EWSopt.silent
    figure;
    subplot(3, 2, 1);
    plot(timeseries(:,1), dataorig(:,2), 'k-','Tag','original data');
    xlim(xlims)
    adjustticklabels(gca,'Y');
    if ~isempty(trend)
        hold on;
        plot(timeseries(:,1), trend, 'r-','Tag','trend');
        hold off;
    end
    if ~isempty(EWSopt.title)
        text(0.1,0.9,EWSopt.title, 'units','normalized','fontsize', 11,'fontname','Arial')
        set(gcf,'name',EWSopt.title)
    end
    xlabel('Time');
    ylabel('variable and trend');
    
    subplot(3, 2, 2);
    if ~isempty(trend)
        h = stem(timeseries(:,1), timeseries(:,2),'k.','Tag','residuals');
        adjustticklabels(gca,'Y');
        text(0.1,0.9,'residuals', 'units','normalized','fontsize', 11,'fontname','Arial')
        set(h, 'markersize', 1);
    else
        text(0.1,0.9,'No residuals - no detrending', 'units','normalized','fontsize', 11,'fontname','Arial')
    end
    xlabel('Time');
    xlim(xlims);
    
    xlabel('Time');
    subplot(3, 2, 3);
    plot(timeseries(:,1), ars, 'k-','tag','AR(1)');
    adjustticklabels(gca,'Y');
    text(0.1,0.9,sprintf('ar(%d)',EWSopt.arlag), 'units','normalized','fontsize', 11,'fontname','Arial')
    if hasstats
        text(0.1,0.1,sprintf('{\\tau} = %5.3g',taus(1)), 'units','normalized','fontsize', 11,'fontname','Arial');
    end
    xlabel('Time');
    xlim(xlims);
    
    subplot(3, 2, 4);
    plot(timeseries(:,1), acfs, 'k-','tag','acfs');
    adjustticklabels(gca,'Y');
    text(0.1,0.9,'acf(1)', 'units','normalized','fontsize', 11,'fontname','Arial')
    if hasstats
        text(0.1,0.1,sprintf('{\\tau} = %5.3g',taus(2)), 'units','normalized','fontsize', 11,'fontname','Arial');
    end
    xlabel('Time');
    xlim(xlims);
    
    subplot(3, 2, 5);
    plot(timeseries(:,1), vars, 'k-','tag','std');
    adjustticklabels(gca,'Y');
    if ~EWSopt.cv
        text(0.1,0.9,'standard deviation', 'units','normalized','fontsize', 11,'fontname','Arial')
    else
        text(0.1,0.9,'coefficient of variation', 'units','normalized','fontsize', 11,'fontname','Arial')
    end
    if hasstats
        text(0.1,0.1,sprintf('{\\tau} = %5.3g',taus(3)), 'units','normalized','fontsize', 11,'fontname','Arial');
    end
    
    
    %   i_plotdefaults
    xlabel('Time');
    xlim(xlims);
    %  ylabel('standard deviation');
    subplot(3, 2, 6);
    plot(timeseries(:,1), skews, 'k-','tag','skewness');
    adjustticklabels(gca,'Y');
    %   i_plotdefaults
    xlabel('Time');
    %   ylabel('skewness');
    text(0.1,0.9,'skewness', 'units','normalized','fontsize', 11,'fontname','Arial')
    if hasstats
        text(0.1,0.1,sprintf('{\\tau} = %5.3g',taus(4)), 'units','normalized','fontsize', 11,'fontname','Arial');
    end
    ch = findobj(gcf,'type','axes');
    set(ch,'TickDir','out');
    adjustpos([0, 0.07]);
    removeaxis('x');
    xlim(xlims);
end
if nargout > 0
    if ~isoctave&&verLessThan('matlab','8.4.0')
        inds = [timeseries(:,1), ars, vars, skews];
    else
        inds = table(timeseries(:,1),timeseries(:,2),trend,ars,acfs,vars,skews);
        if EWSopt.cv
            inds.Properties.VariableNames={'Time','Variable','trend','ar1','acf','cv','skewness'};
        else
            inds.Properties.VariableNames={'Time','Variable','trend','ar1','acf','std','skewness'};
        end
        if isempty(EWSopt.bandwidth)
            % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
            n=size(timeorig,1);
            hx = median(abs(timeorig(:,1) - median(timeorig(:,1)))) / 0.6745 * (4 / 3 / n)^0.2;
            hy = median(abs(timeorig(:,2) - median(timeorig(:,2)))) / 0.6745 * (4 / 3 / n)^0.2;
            EWSopt.bandwidth = sqrt(hy * hx);
        end
        inds.Properties.Description=sprintf('data generated with generic_ews\nwinsize = %g, detrending = %s, bandwidth = %g, logtransform = %d, arlag = %d, interpolate = %d',...
            EWSopt.winsize,EWSopt.detrending,EWSopt.bandwidth,EWSopt.logtransform,EWSopt.arlag,EWSopt.interpolate);
    end
end

%end
function b = linregres(x, y) %#ok<DEFNU>
%b=(sum(x.*y)-1/length(x)*sum(x)*sum(y))/(sum(x.^2)-1/length(x)*sum(x).^2);
b = sum(x .* y) ./ (sum(x.^2)); %without intercept
%intercept
%a=mean(y)-b*mean(x);

function res = skewness(x, flag)
if nargin < 2
    flag = 0;
end
x = x - mean(x);
s2 = mean(x.^2); % this is the biased variance estimator
m3 = mean(x.^3);
res = m3 ./ s2.^(1.5);
if flag == 0
    n = length(x);
    if n > 3
        n(n<3) = NaN; % bias correction is not defined for n < 3.
        res = res .* sqrt((n - 1) ./ n) .* n ./ (n - 2);
    end
end
function res = acf(X, alag)
if nargin < 2
    alag = 1;
end
%autocorrelation function like in R
n = length(X);
s = var(X);
mu = mean(X);
Xt = X(1:end - alag);
Xtk = X(1 + alag:end);
res = 1 / (n - 1) / s*sum((Xt - mu) .* (Xtk - mu));

function res = moving_window(ts, A, lag, funct,opt)
% if nargin <4
%    funct = lag;
%    lag = A;
%    A = ts;
%    ts = (1:size(A, 1))';
% end
res = A .* NaN;
i1 = 1;
i = 1;
if isa(funct,'function_handle')
    sfun=func2str(funct);
else
    sfun=funct;
end
if ts(end) - ts(1) > lag
    tlag = ts + lag;
    while i1 < length(ts)
        while (i < length(ts)) && (ts(i) < tlag(i1))
            i = i + 1;
        end
        if ts(i) >= tlag(i1)
            windowA = A(i1:i - 1, :);
            %write at end/start of the window?
            switch sfun
                case 'cv'
                    res(i - 1, :)=std(opt(i1:i - 1,:))./mean(opt(i1:i - 1,:));
                case 'autoregression'
                    if nargin==5
                        res(i - 1, :) = autoregression(windowA,opt);
                    else
                        res(i - 1, :) = autoregression(windowA,1);
                    end
                otherwise
                    res(i - 1, :) = feval(funct, windowA); %end
            end
            %res(i1,:)= feval(funct, windowA); %start
            %res(floor((i1+i)/2),:)= feval(funct, windowA); %centre
        end
        i1 = i1 + 1;
        i = i - 1;
    end
end


function res = ksmooth(x, y, bandwidth, method)
if nargin < 3
    bandwidth = [];
end
if nargin < 4
    method = 'normal';
end
n = length(x);
if isempty(bandwidth)
    % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
    hx = median(abs(x - median(x))) / 0.6745 * (4 / 3 / n)^0.2;
    hy = median(abs(y - median(y))) / 0.6745 * (4 / 3 / n)^0.2;
    h = sqrt(hy * hx);
    fprintf('bandwidth smoothing set to %g\n', h);
    if h < sqrt(eps) * n
        error('GRIND:outf:LittleVar','There is not enough variation in the data. Regression is meaningless.')
    end
else
    h = bandwidth;
end
switch lower(method(1))
    case 'g' %"Normal" = Gaussian kernel function
        h = h * 0.36055512754640; %variance of 0.13 to get quartiles at + /  -  0.25 * h (comparable with ksmooth (R))
        res = ones(size(y, 1), 1);
        for k = 1:n
            xx = abs(x(k) - x) / h;
            z = exp(-xx(xx < 4).^2 / 2); %x < 4 is more efficient for long time series (negligible effect)
            res(k) = sum(z .* y(xx < 4)) / sum(z);
        end
        %remove tails
        x1=(x(end)-x)>h/2&(x-x(1))>h/2;
        res(~x1)=NaN;
    case 'b' %"box" = moving average
        d = h / 2; % 0.25 quartiles
        res =  cell(1, n);
        for k = 1:n
            xx = abs(x(k) - x) / h;
            z = xx < d;
            res(k) = sum(z .* y) / sum(z);
        end
end


function removeaxis(orientation)
hs = getsubplotgrid(gcf);
if nargin == 0
    orientation = 'b';
end
[rows, cols] = size(hs);
if strncmpi(orientation,'b',1) || strncmpi(orientation,'x',1)
    for i = 1:rows - 1
        for j = 1:cols
            if ishandle(hs(i, j))&&(hs(i, j)~=0)
                set(get(hs(i,j),'xlabel'),'string','')
                set(hs(i, j), 'XTickLabel', []);
            end
        end
    end
    for i = 2:rows
        for j = 1:cols
            if ishandle(hs(i, j))&&(hs(i, j)~=0)
                set(get(hs(i,j),'title'),'string','')
            end
        end
    end
end
if strncmpi(orientation,'b',1) || strncmpi(orientation,'y',1)
    for i = 1:rows
        for j = 2:cols
            if ishandle(hs(i, j))&&(hs(i, j)~=0)
                set(get(hs(i,j),'ylabel'),'string','')
                set(hs(i, j), 'YTickLabel', []);
            end
        end
    end
end
set(hs(:),'fontsize',12)
aa=get(get(gcf,'children'),'xlabel');
for i = 1:length(aa)
    set(aa{i}, 'fontsize', 12);
end
aa=get(get(gcf,'children'),'ylabel');
for i = 1:length(aa)
    set(aa{i}, 'fontsize', 12);
end
function adjustpos(spacing)
set(gcf,'PaperPositionMode','auto')
set(gcf,'units','normalized');
figpos = get(gcf, 'position');
vert = 0.5469; %normal size of figure
hor = 0.41;
hs = getsubplotgrid(gcf);
if nargin == 0
    spacing = 0.02;
end
if length(spacing) == 1
    spacing = spacing + zeros(2, 1);
end
[rows, cols] = size(hs);
if cols<4&&rows<=4
    set(gcf, 'position', [figpos(1:2) hor*cols / 4 vert * rows / 4]);
else
    rat = rows / cols;
    if rat < 1
        set(gcf, 'position', [figpos(1:2) hor vert * rat]);
    else
        set(gcf, 'position', [figpos(1:2) hor / rat vert]);
    end
end
for i = 1:rows
    for j = 1:cols
        h=findobj(gcf,'tag',sprintf('subplot(%d,%d)',i,j));
        if any(ishandle(h))&& any(h~=0)
            set(h, 'position', subplotpos(rows, cols,j, i, spacing));
        end
    end
end
function hs = getsubplotgrid(h)
ch = get(h, 'children');
tags = get(ch, 'tag');
ch = ch(~strcmp(tags, 'legend'));
types = get(ch, 'type');
poss = get(ch, 'position');
ipos = zeros(length(ch), 4);
for i = length(ch):-1:1
    if iscell(types) && strcmp(types{i}, 'axes')
        ipos(i, :) = poss{i};
    elseif ischar(types) && strcmp(types, 'axes')
        ipos(i, :) = poss{i};
    else
        ipos(i, :) = [];
    end
end
colpos = sort(unique(sort(ipos(:, 1))), 'ascend');
rowpos = sort(unique(sort(ipos(:, 2))), 'descend');
hs = zeros(length(rowpos), length(colpos));
for i = 1:length(ch)
    if strcmp(types{i}, 'axes')
        arow=find(rowpos == ipos(i, 2));
        acol=find(colpos == ipos(i, 1));
        hs(arow, acol) = ch(i);
        set(ch(i),'tag',sprintf('subplot(%d,%d)',arow,acol));
    end
end

function adjustticklabels(hax,orient)
axis(hax,'tight')
newtick=get(hax,[orient 'tick']);
tickdiff=(newtick(2)-newtick(1));
newlim=[newtick(1)-tickdiff newtick(end)+tickdiff];
axis(hax,'manual')
set(hax,[orient 'lim'],newlim);
set(hax,[orient 'tick'],[newtick(1)-tickdiff newtick]);




