%LORENZMAP   Create a Lorenz map
%   Plot subsequent maximums of a variable.
%
%   Usage:
%   LORENZMAP VAR1 plot a Lorenz map of the variable VAR1.
%   LORENZMAP plot the variable of the x axis of the phase plane
%
%   LORENZMAP analyses the results of the last run. (if there is
%   no last run or parameters have changed, it calls RU). Use
%   <a href="matlab:help ru">ru</a> if you want to update the last run.
%   LORENZMAP('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'var' [state variable] - state variable to use in the analysis
%
%   See also ru, takens, lyapunov, poincaremap, poincaresect, makemap
%
%   Reference page in Help browser:
%      <a href="matlab:commands('lorenzmap')">commands lorenzmap</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function lorenzmap(varargin)
%avar
global g_Y g_t g_grind t;
fieldnams={'var', 'v', 'state variable to use in the analysis',g_grind.xaxis.var}';
args=i_parseargs(fieldnams,'var','',varargin);
i_parcheck;
if ~isfield(args,'var')
   args.var = g_grind.xaxis.var;
end
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(t, g_grind.ndays, N0, 1);
end
iX = i_varno(args.var);
if isempty(iX)
   [maxima,ndx]=i_maxima(outfun(g_grind.xaxis.var),1);
else
   [maxima,ndx] = i_maxima(g_Y, iX);
end
hfig=i_makefig('lorenzmap');
set(hfig,'name','Lorenz map')
maxlag = i_lagmap(maxima);
hp=plot(maxima,maxlag,'.','Color',g_grind.pen.color);
ud=get(hp,'userdata');
ud.tdata=g_t(ndx(2:end));
ud.plane = '';
ud.xdata = maxima;
ud.ydata = maxlag;
ud.zdata = [];
ud.avar=args.var;
set(hp,'userdata',ud);
i_plotdefaults(hfig);
xlabel(i_disptext([ args.var '_{max, t}']));
ylabel(i_disptext([ args.var '_{max, t+1}']));
hold on;
rang = [min(maxima), max(maxima)];
plot(rang, rang, 'k');
poincaresect('-addreplaydata',get(hfig,'CurrentAxes'));
hold off;






