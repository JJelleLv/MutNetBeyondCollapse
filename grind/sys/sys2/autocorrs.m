%AUTOCORRS   Autocorrelation plot of the state variables
%   Create a plot of the autocorrelation of the state variables 
%   with different time lags. The maximum values of the mean of 
%   the correlations of the state variables is displayed.
%
%   Usage:
%   AUTOCORRS - creates a plot with a maximum time lag of g_grind.ndays/2 days
%   AUTOCORRS MAXLAG - creates a plot with a maximum time lag 
%   of MAXLAG days.
%   AUTOCORRS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'maxlag' [number>0] - maximum time lag
%
%
%   See also takens, simtime
%
%   Reference page in Help browser:
%      <a href="matlab:commands('autocorrs')">commands autocorrs</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function [result,A]=autocorrs(varargin)
%(maxlag)
global g_Y g_t t g_grind;
fieldnams={'maxlag', 'n>0', 'maximum time lag',500}';
args=i_parseargs(fieldnams,'maxlag','',varargin);
i_parcheck;
oldstep=g_grind.tstep;
if isnan(g_grind.tstep)||(g_grind.tstep~=g_grind.ndays)
    g_grind.tstep=g_grind.ndays;
end
N0 = i_initvar;
if i_settingschanged(N0)
    i_ru(t, g_grind.ndays, N0, 1);
end
g_grind.tstep=oldstep;
if ~isfield(args,'maxlag')
    args.maxlag = round(size(g_t, 1) / 2)-1;
end
sizY = size(g_Y, 2);
corr = zeros(args.maxlag + 1, sizY)+nan;
for timestep = 0:args.maxlag
    lag = interp1(g_t, g_Y, g_t - timestep);
    k=find(~isnan(lag(:, 1)),1);
    if ~isempty(k)
        try
            cc = corrcoef([lag(k:size(lag, 1), :), g_Y(k:size(g_Y, 1), :)]);
        catch err
            %      err=lasterror;
            disp(cc);
            rethrow(err);
        end
        if size(cc,1)==2*sizY
            for i = 1:sizY
                corr(timestep + 1, i) = cc(i, sizY + i);
            end
        end
    end
end
hfig = i_makefig('autocorr');
set(hfig,'Name','Auto correlation plot');
if g_grind.statevars.dim>1
    meancorr=mean(corr,2);
    plot(transpose(0:args.maxlag), [corr, meancorr]);
else
    meancorr=corr;
    plot(transpose(0:args.maxlag), corr);
end
i_plotdefaults(hfig);
xlabel('Time lag(t)');
ylabel('Correlation coefficient');
ud = get(gca, 'userdata');
ud.meta=struct('func','autocorr','xname','lag','yname','corr','zname','');
set(gca, 'userdata', ud);

s=i_statevars_names;
for i = 1:g_grind.statevars.dim
    s{i} = i_disptext(char(s{i}));
end
if i>1
    s{i+1}='mean';
    legend(s);
end
A = i_maxima([transpose(0:args.maxlag), meancorr],2);
i=0;
for j=1:size(A,1)
    if A(j,size(A,2))>0.5
        i=i+1;
        A(i,:)=(A(j,:));
    end
end
A=A(1:i,:);
if nargout>0
    if length(corr)>1
        result=corr(2:length(corr));
    end
else
    if isempty(A)
        disp('No periods')
    else
        disp('Periods');
        fprintf('%6s %20s\n','Lag','Mean autocorrelation');
        for j=1:size(A,1)
            fprintf('%6g %20.4g\n',A(j,1),A(j,2));
        end
    end
end

function maxima = i_maxima(A, iX)
maxima = ones(size(A, 1), size(A, 2));
imax = 0;
siz = size(A, 1);
for i = 2:siz - 1
    if (A(i - 1, iX) < A(i, iX)) && (A(i, iX) > A(i + 1, iX))
        imax = imax + 1;
        maxima(imax, :) = A(i, :);
    end
end
maxima = maxima(1:imax, :);




