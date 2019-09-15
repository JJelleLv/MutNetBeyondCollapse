%LYAPUNOV   The maximal Lyapunov exponent
%   Calculate the maximal Lyapunov exponent (lambda) by running a model twice 
%   with slightly different initial conditions. This parameter expresses the
%   sensitivity to initial conditions. A positive Lyapunov
%   exponent is a signature of chaos.
%
%   Usage:
%   LYAPUNOV - opens a dialog for the number of days and disturbance of the 
%   initial conditions of all variables.
%   LYAPUNOV NDAYS - calculate NDAYS days with a default disturbance of 1E-8.
%   LYAPUNOV NDAYS DISTURB = Disturb the initial conditions with DISTURB 
%   and run again.
%   LYAPUNOV('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'disturb' [number>0] - size of the perturbation of the second run.
%     'ndays' [number>0] - number of days for the simulations
%     'silent' [logical] - suppresses figures and output
%   LYAPUNOV('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-s' - silent mode, suppersses figures and output.
%
%   See also lyapspect, lorenzmap, takens, poincaremap, poincaresect, z1test
%
%   Reference page in Help browser:
%      <a href="matlab:commands('lyapunov')">commands lyapunov</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [res,nhoriz] = lyapunov(varargin)
global t g_Y g_t g_grind;
defaultdisturb = 1E-8;
if nargin == 0
    lyapunov('??dlgnoopt');
    return;
    %      i_parcheck;
    %      answer = {num2str(g_grind.ndays), num2str(defaultdisturb)};
    %      prompt={'Number of days for both runs','Difference in initial conditions'};
    %      answer = inputdlg(prompt, 'Lyapunov exponent', 1, answer);
    %      args.ndays = str2double(answer{1});
    %      args.disturb = str2double(answer{2});
    %      args.silent = false;
    %      args.opts={};
else
    fieldnams={'ndays', 'n>0', 'number of days for the simulations',g_grind.ndays;...
        'disturb', 'n>0', 'size of the perturbation of the second run.',defaultdisturb;...
        'silent', 'l', 'suppresses figures and output',false}';
    [args,def]=i_parseargs(fieldnams,'ndays,disturb,silent','-s',varargin);
    i_parcheck;
end
args=mergestructs(def,args);
if any(strcmp(args.opts,'-s'))
    args.silent=true;
end
N0 = i_initvar;
oldtstep = g_grind.tstep;
if isnan(oldtstep)
    g_grind.tstep = args.ndays;
end
try
    i_ru(t, args.ndays, N0, 0);
    Y = g_Y;
    % T = g_t;
    N0 = N0 + args.disturb;
    i_ru(t, args.ndays, N0, 1);
    g_grind.tstep = oldtstep;
catch err
    %   err=lasterror;
    g_grind.tstep = oldtstep;
    rethrow(err);
end

[a, b, dist, nhorizon] = calclyap(g_t, g_Y, Y);
if ~args.silent % Plot results if required
    hfig = i_makefig('lyap1');
    plot(g_t, dist, g_t(1:nhorizon), exp(a * g_t(1:nhorizon) + b));
    i_plotdefaults(hfig);
    set(hfig, 'Name', 'Lyapunov plot');
    if length(g_t) - nhorizon < 10
        ch = '> ';
    else
        ch='= ';
    end
    htitle=title(['\lambda_{max} = ' num2str(a) ';  time horizon \tau ' ch num2str(nhorizon) ]);
    set(htitle,'fontweight','normal');
    xlabel('Time (t)');
    ylabel('difference between 2 runs (log scale)');
    set(gca, 'YScale', 'Log');
    ud = get(gca, 'userdata');
    ud.meta=struct('func','lyapunov','xname','t','xlim',[min(g_t) max(g_t)],'yname','difference','zname','');
    set(gca, 'userdata', ud);
    
    hfig = i_makefig('lyap2');
    %H = figure(i_figno('lyap2'));
    plot(g_t, Y(:, 1), g_t, g_Y(:, 1));
    i_plotdefaults(hfig);
    htitle=title('Time plot with 2 slightly different initial settings');
    ud = get(gca, 'userdata');
    ud.meta=struct('func','lyapunov','xname','t','xlim',[min(g_t) max(g_t)],'yname',i_statevars_names(1),'zname','');
    set(gca, 'userdata', ud);
    set(htitle,'fontweight','normal');
    set(hfig, 'Name', 'Lyapunov time plot');
    xlabel('Time (t)');
    ylabel(i_disptext(i_statevars_names(1)));
    legend(sprintf('%s 1^{st} run', i_disptext(i_statevars_names(1))), ...
        sprintf('%s 2^{nd} run', i_disptext(i_statevars_names(1))))
elseif nargout == 0
    fprintf('n horizon: %g\nLyapunov exponent:  %g\n',nhorizon,a);
end

if nargout > 0
    res = a;
    nhoriz = nhorizon;
end

%%%% function calclyap
function [a, b, dist, nhorizon] = calclyap(gg_t, Y1, Y2)
dist = sqrt(sum((Y1 - Y2).^2, 2));
%dists=abs(Y1-Y2);
%dists(dists<1E-16)=1E-16;
maxdist = max(dist);
if maxdist < 1E-3 * max(max(Y1))
    maxdist = 1E-3 * max(max(Y1));
end
nhor = -1;
for tt = 1:length(gg_t)
    if abs(dist(tt)) < 1E-16
        warning('GRIND:lyapunov:nodiff','Small difference between runs: log of (almost) zero');
        dist(tt) = 1E-16;
    end
    if (nhor==-1)&&(dist(tt) > maxdist / 4)
        nhor = tt;
    end
end
if nhor ~= -1
    nhorizon = nhor;
else
    nhorizon = length(gg_t);
end
n = nhorizon;
%tim = transpose([1:ndays]);
logd = log(dist(1:nhorizon));
tt1 = gg_t(1:nhorizon);
sumx = sum(tt1);
sumy = sum(logd);
sumx2 = sum(tt1 .* tt1);
sumxy = sum(tt1 .*  logd);
a = (sumxy - sumx * sumy / n) / (sumx2 - sumx * sumx / n);
b = sumy / n - a * sumx / n;

