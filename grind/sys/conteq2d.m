%CONTEQ2D   Continue an equilibrium in two dimensions
%   Depreciated function, will be removed. Use <a href="matlab:help contbif">contbif engine grind</a> (or engine cocogrind) instead.
%  
%   This command can do a simple 2D bifurcation analysis by repeating a 1D bifurcation
%   analysis several times (NGRID). Use this engine only if coco does not work.
%   There are two ways of doing this:
%   1) engine cocogrind (default) - for the 1D bifuction coco is used.
%   2) engine grind - for the 1D bifuction the GRIND engine is used.
%
%   Usage:
%   CONTEQ2D -  - user is prompted for information.
%   See <a href="matlab:help contbif">CONTBIF</a> for all other command line options
%
%   See also contbif, conteq, open_matcont_gui
%
%   Reference page in Help browser:
%      <a href="matlab:commands('conteq2d')">commands conteq2d</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function conteq2d(varargin)
hasengine=false;
for i=1:length(varargin)
    if ischar(varargin{i})&&strcmp(varargin{i},'engine')
        hasengine=true;
        break;
    end
end
if hasengine
    conteq(varargin{:});
else
    conteq(varargin{:},'engine','cocogrind');
end
% i_parcheck;
% global g_grind g_cont;
% if isfield(g_grind,'pars')&&isempty(g_grind.pars)
%    error('GRIND:conteq2d:NoParameters','no parameters to analyse');
% end
% if isempty(g_cont)
%     g_cont=grind_coco;
% end
% if nargin == 0
%    options={};
%    if isfield(g_grind, 'conteq2d') && ~isempty(g_grind.conteq2d)
%       settings = g_grind.conteq2d;
%       pars=settings.allpars(settings.freepars);
%       if isempty(settings.parranges)
%           settings.parranges={[0, 10] [0, 10]};
%       elseif length(settings.parranges)==1
%           settings.parranges{2}=[0 10];
%       end
%       answer{1} = pars{1};
%       answer{2} = num2str(settings.parranges{1}(1));
%       answer{3} = num2str(settings.parranges{1}(2));
%       answer{4} = num2str(settings.relh);
%       answer{5} = pars{2};
%       answer{6} = num2str(settings.parranges{2}(1));
%       answer{7} = num2str(settings.parranges{2}(2));
%       answer{8} = num2str(settings.ngrid);
%    else
%       answer={'','0','10','1','','0','10','10'};
%    end
%    prompt={'parameter xaxis','Minimum','Maximum','Relative step size (relh) (conteq)','Parameter yaxis','Minimum','Maximum','Number of steps (ngrid)'};
%    answer = inputdlg(prompt, 'conteq 2D', 1, answer);
%    if isempty(answer)
%        disp('Cancel pressed')
%        return;
%    else
%        settings=g_cont.settings;
%        if isempty(settings.parranges)
%            settings.parranges={[nan, nan] [nan, nan]};
%        end
%        settings.freepars=[0 0];
%        settings.freepars(1)=find(strcmp(settings.allpars,answer{1}));
%        settings.parranges{1}(1) = str2double(answer{2});
%        settings.parranges{1}(2) = str2double(answer{3});
%        settings.relh = str2double(answer{4});
%        settings.freepars(2)=find(strcmp(settings.allpars,answer{5}));
%        settings.parranges{2}(1) = str2double(answer{6});
%        settings.parranges{2}(2) = str2double(answer{7});
%        settings.ngrid = str2double(answer{8});
%    end
% else
%    [settings, options] = i_parsearguments(g_cont,g_cont.defaultsettings,varargin,true);
%   
%    g_cont.settings=settings;
%    if ~isfield(settings,'ngrid')
%        settings.ngrid=10;
%    end
% end
% i_run_conteq2d(settings,options)
% function i_run_conteq2d(settings,options)
% global g_cont;
% if any(strcmp(options,'-1'))
%     plotresults;
%     return;
% end
% if ~isfield(settings,'ngrid')
%     settings.ngrid=10;
% end
% oldg_cont=grind_coco;
% oldg_cont.assign(g_cont); %as g_cont is a handle you need to deep copy it
% if isempty(g_cont)
%     g_cont=grind_coco;
% end
% hwait = i_waitbar(0, settings.ngrid+1, 'conteq2d', sprintf('Running conteq 1/%d',settings.ngrid+1), 0.5, 1);
% %The cancelbutton sets the userdata of the waitbar to 1
% %make sure the program reacts to that
% g_grind.conteq2d = settings;
% pars = settings.allpars(settings.freepars);
% conteqcomm=struct('silent',true,'N0',[],'par1',pars{1}, 'parranges',settings.parranges{1},'NPR',100);
% if ishandle(hwait)&&get(hwait,'userdata')==1
%     g_cont.assign(oldg_cont);
%     return;
% end
% i_waitbar(1);
% conteq(conteqcomm);
% vcont_run = g_cont.run;
% n = length(vcont_run.parvalues);
% vparvals=vcont_run.parvalues;
% oldp1 = evalin('base', pars{1});
% oldp2 = evalin('base', pars{2});
% conteqcomm=struct('silent',true,'N0',[],'par1',pars{2}, 'parranges',settings.parranges{2},'NPR',100);
% parvals2=linspace(settings.parranges{1}(1), settings.parranges{1}(2),settings.ngrid);
% try
%    g_cont.conteq2d.points=[];
%    g_cont.conteq2d.settings=settings;
%    for i = 1:settings.ngrid
%        p=parvals2(i);
%        ndxs=find(diff([0;vparvals>p])~=0);
%        i_waitbar(1,sprintf('Running conteq %d/%d',i+1,settings.ngrid+1));
%        for k=1:length(ndxs)
%            p = vcont_run.parvalues(ndxs(k));
%            N1 = transpose(vcont_run.Y(1,:,ndxs(k)));
%            assignin('base', pars{1}, p);
%            conteqcomm.N0=N1;
%            conteq(conteqcomm);
%            if ishandle(hwait)&&get(hwait,'userdata')==1
%                error('grind:conteq2d','Cancel pressed');
%            end
%            g_cont.conteq2d.points=[g_cont.conteq2d.points, g_cont.findbif({'T','F','H'})];
%        end
%    end
%    i_waitbar([]);
%    plotresults;
%    assignin('base', pars{1}, oldp1);
%    assignin('base', pars{2}, oldp2);
%    c2d=g_cont.conteq2d;
%    g_cont.assign(oldg_cont);
%    g_cont.conteq2d=c2d;
% catch err
%    i_waitbar([]);
%    assignin('base', pars{1}, oldp1);
%    assignin('base', pars{2}, oldp2);
%    c2d=g_cont.conteq2d;
%    g_cont.assign(oldg_cont);
%    g_cont.conteq2d=c2d;
%    rethrow(err);
% end
% function plotresults
% global g_cont;
%    [hfig,isnew]=i_makefig('conteq2d');
%    oldhold=ishold;
%    if isnew
%        hold on;
%    end
%    ss={g_cont.conteq2d.points.id};
%    pars=g_cont.conteq2d.settings.allpars(g_cont.conteq2d.settings.freepars);
%    p1=g_cont.conteq2d.settings.freepars(1);
%    p2=g_cont.conteq2d.settings.freepars(2);
%    ndx=strncmp(ss,'H',1);
%    parvalues=[g_cont.conteq2d.points.p0];
%    leg = {};
%    if any(ndx)
%       plot(parvalues(p1, ndx), parvalues(p2, ndx), 'or');
%       hold on;
%       leg = {'Hopf'};
%    end
%    ndx=strncmp(ss,'F',1);
%    if any(ndx)
%       plot(parvalues(p1, ndx), parvalues(p2, ndx), 'ob');
%       hold on;
%       leg = [leg;{'Fold'}];
%    end
%    ndx=strncmp(ss,'T',1);
%    if any(ndx)
%       plot(parvalues(p1, ndx), parvalues(p2, ndx), 'og');
%       hold on;
%       leg = [leg; {'Transcritical'}];
%    end
%    leg=i_addlegend(leg);
%    if length(leg)>1
%       legend(leg);
%    elseif length(leg)==1
%       legend(leg{1});
%    end
%    i_plotdefaults(hfig);
%    xlabel(pars{1});
%    ylabel(pars{2});
%    ud = get(gca, 'userdata');
%    ud.meta=struct('func','conteq2d','xname',pars{1},...
%        'yname',pars{2},'zname','');  
%    set(gca, 'userdata', ud);
%    if ~oldhold
%      hold off;
%    end
