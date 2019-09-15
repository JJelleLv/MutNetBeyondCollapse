%TAKENS   Create a Takens plot 
%   Plot of N(t) versus N(t+timestep). TAKENS analyses the results of the 
%   last run. (if there is no last run or if parameters have changed, it calls RU). 
%   Use <a href="matlab:help ru">ru</a> if you want to update the last run.
%
%   Usage:
%   TAKENS uses a default time lag  of 1 time step.
%   TAKENS TIMELAG uses a time lag of TIMELAG time units.
%   TAKENS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'timelag' [number>0] - time lag for the Takens plot
%
%
%   See also lorenzmap, lyapunov, poincaresect, poincaremap, makemap
%
%   Reference page in Help browser:
%      <a href="matlab:commands('takens')">commands takens</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function takens(varargin)
global g_Y g_t t g_grind;
fieldnams={'timelag', 'n>0', 'time lag for the Takens plot',1}';
args=i_parseargs(fieldnams,'timelag','',varargin);
i_parcheck;
if ~isfield(args,'timelag')
   args.timelag = 1;
end
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(t, g_grind.ndays, N0, 1);
end
lag=interp1(g_t,g_Y,g_t+args.timelag);
hfig=i_makefig('takens');
oldhold = ishold;
hold on;
ud = get(gca, 'userdata');
ud.meta=struct('func','takens','xname','statevars','yname','statevars(t+tau)','zname','');  
set(gca, 'userdata', ud);
leg = cell(g_grind.statevars.dim, 1);
Hs = zeros(g_grind.statevars.dim, 1);
for iX = 1:size(g_Y, 2)
   hp = plot(g_Y(:, iX), lag(:, iX), '.');
   set(hp, 'Color', g_grind.pen.color2);
   Hs(iX)=hp;
   leg{iX}=i_statevars_names(iX);
   nextpen;
end
i_plotdefaults(hfig);
legend(Hs,leg);
xlabel('X_t');
ylabel(['X_{t+' num2str(args.timelag) ')}']);
if ~oldhold
   hold off;
end





