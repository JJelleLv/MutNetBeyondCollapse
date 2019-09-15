% internal function of parameter analyser
function i_paranal(plotprev, conteq_eq)
if nargin < 2
   conteq_eq = 0;
end

%parname, start, nend, nsteps, nstabilizing, ndays, lines, plotprev, outputtype, disturb)
%   i_paranal(answer.par{1}, answer.start(1), answer.nend(1), answer.steps(1), answer.stabil, ...
%      answer.writing,answer.lines,plotprev, answer.outputtype, answer.disturbance);

global g_grind g_Y g_t g_paranal;
oldY = g_Y;
oldt = g_t;
oldstep = g_grind.tstep;
try
   scenario.pars = g_grind.paranal.dlg.par(1);
   scenario.parvalues = transpose(linspace(g_grind.paranal.dlg.start(1), g_grind.paranal.dlg.nend(1), g_grind.paranal.dlg.steps(1) + 1));
   % nend = g_grind.paranal.dlg.nend(1);
   % nsteps = g_grind.paranal.dlg.steps(1);
   % ppars=linspace(start,nend,nsteps+1); %more precise than summing per step
   scenario.nstabilizing = g_grind.paranal.dlg.stabil;
   scenario.initialstabil = g_grind.paranal.dlg.stabil;
   scenario.nwriting = g_grind.paranal.dlg.writing;
   scenario.avoidextiction = true;
   if isfield(g_grind.paranal.dlg, 'perturbcallback')
      if ischar(g_grind.paranal.dlg.perturbcallback)
         g_grind.paranal.dlg.perturbcallback = str2func(g_grind.paranal.dlg.perturbcallback);
      end

      scenario.perturbcallback =  g_grind.paranal.dlg.perturbcallback;
   end

   scenario.reset  = false;
   lines = g_grind.paranal.dlg.lines;
   scenario.outputtype = g_grind.paranal.dlg.outputtype;
   silent = 0;
   if isfield(g_grind.paranal, 'silent') &&( g_grind.paranal.silent==1)
      g_grind.paranal.silent = 0;
      silent = 1;
   end

   if ~silent
      scenario.stepcallback =  @stepcallback;
   end

   %disturb=g_grind.paranal.dlg.disturbance;
   if length(g_grind.paranal.dlg.par) > 1&&~isempty(g_grind.paranal.dlg.par{2})
      scenario.pars{2} = g_grind.paranal.dlg.par{2};
      steps1 = length(scenario.parvalues);
      pval2 = transpose(linspace(g_grind.paranal.dlg.start(2), g_grind.paranal.dlg.nend(2), g_grind.paranal.dlg.steps(2) + 1));
      steps2 = length(pval2);
      pval2 = transpose(repmat(pval2, 1, length(scenario.parvalues)));
      scenario.parvalues = repmat(scenario.parvalues, 1, steps2);
      scenario.parvalues = [scenario.parvalues(:), pval2(:)];
      % for 2d paranal reset at each new value of parameter 2
      scenario.reset = false(size(scenario.parvalues, 1), 1);
      scenario.reset(1:steps1:size(scenario.parvalues, 1)) = true;
      plotsfield = 'plots2d';
   else
      plotsfield = 'plots';
      if lines > 1
         lines = 0;
      end

   end

   
   if ~isfield(g_grind, 'paranal') || ~isfield(g_grind.paranal, 'defaultplots') || g_grind.paranal.defaultplots
      paranal('-defaults');
   end

   parlim = [min(scenario.parvalues); max(scenario.parvalues)];
   for i = 1:length(g_grind.paranal.(plotsfield))
      if strcmp(g_grind.paranal.(plotsfield){i}.xaxis{1}, '<param1>')
         g_grind.paranal.(plotsfield){i}.xlim = parlim(:, 1);
      elseif strcmp(g_grind.paranal.(plotsfield){i}.xaxis{1}, '<param2>')
         g_grind.paranal.(plotsfield){i}.xlim = parlim(:, 2);
      end

      if strcmp(g_grind.paranal.(plotsfield){i}.yaxis{1}, '<param1>')
         g_grind.paranal.(plotsfield){i}.ylim = parlim(:, 1);
      elseif strcmp(g_grind.paranal.(plotsfield){i}.yaxis{1}, '<param2>')
         g_grind.paranal.(plotsfield){i}.ylim = parlim(:, 2);
      end

      if strcmp(g_grind.paranal.(plotsfield){i}.zaxis{1}, '<param1>')
         g_grind.paranal.(plotsfield){i}.zlim = parlim(:, 1);
      elseif strcmp(g_grind.paranal.(plotsfield){i}.zaxis{1}, '<param2>')
         g_grind.paranal.(plotsfield){i}.zlim = parlim(:, 2);
      end

   end

   
   
   outputlist = i_paranaldialog('outputlist');
   if ~ischar(scenario.outputtype)
      scenario.outputtype = lower(outputlist{scenario.outputtype});
   end

   % disterror = 1E-4;
   % if ~isempty(disturb)
   %    if iscell(disturb)
   %       disturb = disturb{1};
   %    end

   %    if ~isempty(disturb) & isempty(strfind(disturb, ';'))
   %       disturb = [disturb ';'];
   %    end

   % end

   g_grind.solver.opt.OutputFcn = [];
   % iY = i_varno(g_grind.paranal.(plotsfield){1}.yaxis{1});
   % iZ = i_varno(g_grind.paranal.(plotsfield){1}.zaxis{1});
   if silent
      hfig = [];
   else
      %   parfig = i_figno('paranal');
      [hfig, isnew] = i_makefig('paranal');
      if isnew
         hold on;
      end

      set(hfig, 'Name', 'Parameter analysis');
      ud = get(gca, 'userdata');
      ud.meta=struct('func','paranal','xname',g_grind.paranal.(plotsfield){1}.xaxis{1},...
          'yname',g_grind.paranal.(plotsfield){1}.yaxis{1},'zname',g_grind.paranal.(plotsfield){1}.zaxis{1});  
      set(gca, 'userdata', ud);

      scenario.hfig = hfig;
      if  ~isempty(g_grind.paranal.(plotsfield){1}.zaxis{1}) && isnew
         set(gca, 'View', [8.5 10])
         hold on;
      elseif isnew
         gca;
      end

      i_plotdefaults(hfig);
      if ~ishold
         delete(get(gca, 'children'))
         set(gca, 'userdata', [])
         legend off
      end

      i_odespeed(0, 0, 'init');
      % plotedit off;
      ch=findobj(gca,'type','line');
      g_grind.paranal.pen = i_nextpen([]);
      for i = 1:length(ch)
         g_grind.paranal.pen = i_nextpen(g_grind.paranal.pen);
      end

      if lines == 1
         ppen = g_grind.paranal.pen.pen2;
      else
         ppen = '.';
      end
   end

   if plotprev ~= 1
      if isempty(g_grind.paranal.(plotsfield))
         scenario.xaxis1 = '';
         scenario.yaxis1 = '';
         scenario.zaxis1 = '';
      else
         scenario.xaxis1 = strrep(g_grind.paranal.(plotsfield){1}.xaxis{1}, '<param1>', scenario.pars{1});
         if ~isempty(g_grind.paranal.(plotsfield){1}.yaxis{1})
            scenario.yaxis1 = strrep(g_grind.paranal.(plotsfield){1}.yaxis{1}, '<param1>', scenario.pars{1});
         else
            scenario.yaxis1 = '';
         end

         if ~isempty(g_grind.paranal.(plotsfield){1}.zaxis{1})
            scenario.zaxis1 = strrep(g_grind.paranal.(plotsfield){1}.zaxis{1}, '<param1>', scenario.pars{1});
         else
            scenario.zaxis1 = '';
         end

         if length(scenario.pars) > 1
            scenario.xaxis1 = strrep(scenario.xaxis1, '<param2>', scenario.pars{2});
            if ~isempty(g_grind.paranal.(plotsfield){1}.yaxis{1})
               scenario.yaxis1 = strrep(scenario.yaxis1, '<param2>', scenario.pars{2});
            end

            if ~isempty(g_grind.paranal.(plotsfield){1}.zaxis{1})
               scenario.zaxis1 = strrep(scenario.zaxis1, '<param2>', scenario.pars{2});
            end

         end

      end

      if length(scenario.pars) > 1
         i_waitbar(0, size(scenario.parvalues,1), 'grind', 'Calculating',0.5)
      end

      g_paranal.run =  i_stepscenario(scenario);
      i_waitbar([]);
      if ishandle(hfig)
         h1=findobj(hfig,'tag','tmpdraw');
         if ishandle(h1)
            delete(h1);
         end

      end

      
   end
   
   %gather all data
   %[ps, funx,funy, iY] = i_paranalfun(1,iY, iZ, outputtype, parname, g_paranal.run.t, g_paranal.run.Y, g_paranal.run.p, disterror);
   %h = i_makefig('paranal');
   if ~silent
      i_odespeed(1, 1, 'done');
      for i = 1:length(g_grind.paranal.(plotsfield))
         nx = length(g_grind.paranal.(plotsfield){i}.xaxis);
         ny = length(g_grind.paranal.(plotsfield){i}.yaxis);
         nz = length(g_grind.paranal.(plotsfield){i}.zaxis);
         [hfig,isnew] = i_makefig('paranal', i - 1);
         i_plotdefaults();
         set(gca, 'YLimMode', 'auto');
         set(gca, 'ZLimMode', 'auto');
         set(hfig, 'userdata', []);
         if isnew
            hold('on');
         end

         for k = 1:max([nx, ny, nz])
            ps = outfun(g_grind.paranal.(plotsfield){i}.xaxis{min(nx, k)}, '-p');
            eq1 = g_grind.paranal.(plotsfield){i}.yaxis{min(ny, k)};
            [funx,ndx] = outfun(eq1,'-p', scenario.outputtype);
            if ~isempty(ndx) && (length(ndx)~=length(ps))
               ps = ps(ndx);
            end

            if ~isempty(g_grind.paranal.(plotsfield){i}.zaxis{min(nz, k)})
               funy = outfun(g_grind.paranal.(plotsfield){i}.zaxis{min(nz,k)},'-p', scenario.outputtype);
               if length(funy) ~= length(funx)
                  funy = outfun(g_grind.paranal.(plotsfield){i}.zaxis{min(nz, k)}, '-p');
                  funy = funy(ndx, :);
               end

            else
               funy = [];
            end

            %      iY = i_varno(g_grind.paranal.(plotsfield){i}.yaxis{1});
            %      iZ = i_varno(g_grind.paranal.(plotsfield){i}.zaxis{1});
            %      [ps, funx,funy, iY] = i_paranalfun(i, iY, iZ, outputtype, parname, g_paranal.run.t, g_paranal.run.Y, g_paranal.run.p, iif(hasperm,g_paranal.run.perm,[]), disterror);
            g_grind.paranal.pen = i_nextpen(g_grind.paranal.pen);
            if lines  >= 2
               [X, Y, Z] = makegrids(ps, funx, funy);
            end

            if lines  == 2||lines==3
               if lines == 2
                  hplot = contourf(X, Y, Z);
                  set(gca, 'view', [0 90])
               else
                  hplot = surf(X, Y, Z);
               end
               htitle = title(sprintf('Mean of %s',mydisptext(g_grind.paranal.(plotsfield){i}.zaxis{k}, scenario.pars)));
               set(htitle,'fontweight','normal');
               colorbar;
            else
               if length(g_paranal.run.pars) > 1
                  [ps, funx, funy] = addbreaks(g_paranal, ps, funx, funy);
               end

               hplot = makeplot(i, ppen, ps, funx, funy);
               %          set(hplot,'EraseMode', 'normal')
               addreplaydata(gca, i);
            end

            nser = max(size(funx, 2), size(funy, 2));
            increasingpar = sign(g_grind.paranal.dlg.nend(1) - g_grind.paranal.dlg.start(1));
            if conteq_eq > 0
               set(hplot,'tag','conteq');
               increasingpar = conteq_eq;
               if any(conteq_eq ==[2 4])
                  sers=findobj(gca,'type','line');
                  if ~isempty(sers)
                     set(sers(1:nser),'LineStyle',':')
                  end

               elseif conteq_eq==5
                  sers=findobj(gca,'type','line');
                  if ~isempty(sers)
                     set(sers(1:nser),'LineStyle','None')
                     set(sers(1:nser),'Marker','o')
                     set(sers(1:nser),'MarkerFacecolor','r','MarkerEdgecolor','r')
                  end
               elseif conteq_eq==6
                  sers=findobj(gca,'type','line');
                  if ~isempty(sers)
                     set(sers(1:nser),'LineStyle','None')
                     set(sers(1:nser),'Marker','o')
                     set(sers(1:nser),'MarkerFacecolor','b','MarkerEdgecolor','b')
                 end

               elseif conteq_eq==7
                  sers=findobj(gca,'type','line');
                  if ~isempty(sers)
                     set(sers(1:nser),'LineStyle','None')
                     set(sers(1:nser),'Marker','s')
                     set(sers(1:nser),'MarkerFacecolor','g','MarkerEdgecolor','g')
                  end
               end

            end

            
            if nser == 1
               if isempty(g_grind.paranal.(plotsfield){i}.zaxis{1})
                  i_grindlegend(2,hplot,mydisptext(sprintf('%s', eq1)),{'g_increasingpar',increasingpar});
               else
                  i_grindlegend(2,hplot,' ',{'g_increasingpar',increasingpar});
               end

            else
               texts=cell(nser,1);
               for j = 1:nser
                  if isempty(g_grind.paranal.(plotsfield){i}.zaxis{1})&&j<=ny
                     texts{i}=mydisptext(sprintf('%s(%d)', eq1,j));
                  else
                     texts{i}=' ';
                  end

                end

                i_grindlegend(2,hplot,texts,{'g_increasingpar',increasingpar});
%                for j = 1:nser
%                   if isempty(g_grind.paranal.(plotsfield){i}.zaxis{1})&&j<=ny
%                      i_grindlegend(2,hplot,mydisptext(sprintf('%s(%d)', eq1,j)),{'g_increasingpar',increasingpar});
%                   else
%                      i_grindlegend(2,hplot,' ',{'g_increasingpar',increasingpar});
%                   end

%                end

            end

         end

         set(hfig, 'Name', sprintf('Parameter analysis - %d',i));
         if ~isempty(g_grind.paranal.(plotsfield){i}.zaxis{1})&&(lines~=2)
            set(gca, 'View', [8.5 10])
         end
         if ~isempty(g_grind.paranal.(plotsfield){i}.xlim)&&~any(isnan(g_grind.paranal.(plotsfield){i}.xlim))&&...
             g_grind.paranal.(plotsfield){i}.xlim(1)~=g_grind.paranal.(plotsfield){i}.xlim(2)
            xlim(real(g_grind.paranal.(plotsfield){i}.xlim));
         end

         xlabel(mydisptext(g_grind.paranal.(plotsfield){i}.xaxis{1}, scenario.pars));
         if ~strcmp(scenario.outputtype, 'unchanged')
            ylabel(mydisptext(sprintf('%s_{%s}', g_grind.paranal.(plotsfield){i}.yaxis{1}, scenario.outputtype),scenario.pars));
            zlabel(mydisptext(sprintf('%s_{%s}',g_grind.paranal.(plotsfield){i}.zaxis{1}, scenario.outputtype),scenario.pars));
         else
            ylabel(mydisptext(g_grind.paranal.(plotsfield){i}.yaxis{1}, scenario.pars));
            zlabel(mydisptext(g_grind.paranal.(plotsfield){i}.zaxis{1}, scenario.pars));
         end

         set(findobj(hfig,'type','line'),'linewidth',g_grind.paranal.pen.linewidth);
         i_plotdefaults(hfig);
         %    plotedit on;
         i_grindlegend(11);
         %   if ~oldhold
         %      hold off;
         %   end

      end

   end

   g_Y = oldY;
   g_t = oldt;
   g_grind.tstep = oldstep;
catch err
   g_Y = oldY;
   g_t = oldt;
   g_grind.tstep = oldstep;
   rethrow(err);
end

function dostop = stepcallback(i, scenario)
global g_grind
if length(scenario.pars) == 2
   plotsfield = 'plots2d';
else
   plotsfield = 'plots';
end
% fieldnams={'fun', 'U1', 'The equation(s) to show';...
%    'datastruc', 'r', 'You can use a data structure like g_paranal.run to extract the data';...
%    'catfun', 'r#s#n', 'Get categorial statistics of the state variable';...
%    'times', 'n', 'The times for outputs';...
%    'runfield', 's', 'The name of the field in the datastruc that contains the run information (e.g. ''run'' or ''prevrun'')'}';
% args=i_parseargs(fieldnams,'fun,catfun','-#ndxofparanal,-n,-p,-b,-c,-m,-t',varargin,false,{@i_is_outfun_equation});

i_waitbar(1);
dostop = false;
hfig = get(0, 'currentfigure');
if ishandle(hfig)
   xdata = outfun(scenario.xaxis1);
   [ydata,ndx] = outfun(struct('fun',scenario.yaxis1,'opts',{{'-n'}},'catfun',scenario.outputtype));
   if ~isempty(ndx) && (length(ndx)~=length(xdata))
      xdata = xdata(ndx);
   end

   if ~isempty(scenario.zaxis1)
      zdata = outfun(struct('fun',scenario.zaxis1,'opts',{{'-n'}},'catfun',scenario.outputtype));
      if length(ydata) ~= length(zdata)
         zdata = outfun(struct('fun',scenario.zaxis1));
         zdata = zdata(ndx, :);
      end

   else
      zdata = [];
   end
   if any(~isreal(ydata))||any(~isreal(ydata))||any(~isreal(xdata))
       warning('grind:paranal:complex','Imaginary parts of complex values ignored');
       xdata=real(xdata);
       ydata=real(ydata);
       zdata=real(zdata);
   end
   if g_grind.drawnow && (hfig == scenario.hfig)
      ud = get(scenario.hfig, 'userdata');
      if ~isempty(ud) && isfield(ud, 'stop') && ud.stop
         dostop = true;
      end

      if i == 1
         if isempty(zdata)
            [h] = makeplot(1,'.', xdata(:,1),ydata(:,1),zdata);
         else
            [h] = makeplot(1,'.', xdata(:,1),ydata(:,1),zdata(:,1));
         end

         set(h,'tag','tmpdraw');
         hax = gca;
         if ~isempty(g_grind.paranal.(plotsfield){1}.xlim)
            lim=g_grind.paranal.(plotsfield){1}.xlim;
            if lim(2)==lim(1)
                lim(2)=lim(1)+1;
            end
            set(hax, 'Xlim', real(lim));
         end

         if ~isempty(g_grind.paranal.(plotsfield){1}.ylim)
            lim=g_grind.paranal.(plotsfield){1}.ylim;
            if lim(2)==lim(1)
                lim(2)=lim(1)+1;
            end

            set(hax, 'Ylim', real(lim));
         end

         if ~isempty(g_grind.paranal.(plotsfield){1}.zlim)
            set(hax, 'Zlim', real(g_grind.paranal.(plotsfield){1}.zlim));
         end

         set(h, 'Color', g_grind.paranal.pen.drawcolor);
         xlabel(mydisptext(g_grind.paranal.(plotsfield){1}.xaxis{1}, scenario.pars));
         ylabel(mydisptext(g_grind.paranal.(plotsfield){1}.yaxis{1},  scenario.pars));
         if ~isempty(g_grind.paranal.(plotsfield){1}.zaxis{1})
            zlabel(mydisptext(g_grind.paranal.(plotsfield){1}.zaxis{1},  scenario.pars));
         end

      else
         h=findobj(scenario.hfig,'tag','tmpdraw');
         if length(h) > 1
            h = h(1);
         end

         if ishandle(h)
            prevx = get(h, 'xdata');
            prevy = get(h, 'ydata');
            if isempty(zdata)
               set(h, 'xdata', [prevx(:);xdata(:,1)], 'Ydata', [prevy(:);ydata(:,1)],  'Color', g_grind.paranal.pen.drawcolor);
            else
               prevz = get(h, 'zdata');
               set(h, 'xdata', [prevx(:);xdata(:,1)], 'Ydata', [prevy(:);ydata(:,1)], 'zdata', [prevz(:);zdata(:,1)], 'Color', g_grind.paranal.pen.drawcolor);
            end

         end

      end

   end

end

drawnow;
function s = mydisptext(s, parname)
if strcmp(s, '<param1>')
   s = parname{1};
elseif strcmp(s, '<param2>')
   s = parname{2};
end
s = i_disptext(s);

function h = makeplot(~,ppen, ps, funx, funy,varargin)
global g_grind
maxn = max([size(ps, 2), size(funx, 2), size(funy, 2)]);
oldversion=~isoctave&&verLessThan('matlab','8.4.0');
oldhold = ishold;
% if ~oldhold
%      delete(get(gca,'children'))
%      set(gca,'userdata',[])
%      legend off
% end

hold on;
apen = g_grind.paranal.pen;
if ~isempty(funy)
   for i = 1:maxn
      h = plot3(ps(:,min(size(ps,2),i)), funx(:,min(size(funx,2),i)), funy(:,min(size(funy,2),i)), ppen, ...
         'MarkerSize', 3,  'Color', apen.color2,varargin{:});
      if ps(1) < ps(end)
         set(h,'tag','paranal+');
      else
         set(h,'tag','paranal-');
      end

      if oldversion
         set(gca, 'drawmode','fast');
      else
         set(gca,'SortMethod','depth');
      end

      apen = i_nextpen(apen);
      if ppen ~= '.'
         ppen = apen.pen2;
      end

      box on;
   end

elseif ~isempty(funx)
   for i = 1:maxn
      h = plot(ps(:,min(size(ps,2),i)), funx(:,min(size(funx,2),i)), ppen, ...
         'MarkerSize', 3, 'Color', apen.color2,varargin{:});
      if ps(1) < ps(end)
         set(h,'tag','paranal+');
      else
         set(h,'tag','paranal-');
      end

      apen = i_nextpen(apen);
      if ppen ~= '.'
         ppen = apen.pen2;
      end

   end

else
   h = [];
   warning('GRIND:paranal','Cannot plot empty series');
   return;
end

if ~oldhold
   hold off;
end


function addreplaydata(h, No)
global g_paranal;
ud = get(h, 'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
ud.replay.onturn = @i_replayparanalturn;
ud.replay.pars = g_paranal.run.parvalues(:, 1);
h1=findobj(h,'tag','hwhere');
if ~isempty(h1)
   delete(h1);
end

ud.replay.hwhere = -1;
tstep = (ud.replay.pars(end) - ud.replay.pars(1)) / length(ud.replay.pars);
ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',size(g_paranal.run.parvalues,1));
ud.replay.opt = No;
set(h,'userdata', ud);


function replaystart(hax)
global g_grind g_paranal;
if ishandle(hax)&&~isempty(g_paranal)&&~isempty(g_paranal.run)
   i_figure(get(hax, 'parent'));
   if length(g_paranal.run.pars) > 1
      plotsfield = 'plots2d';
   else
      plotsfield = 'plots';
   end

   set(hax,'ylimmode','manual');
   set(hax,'zlimmode','manual');
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')&&(isfield(g_paranal,'run')&&~isempty(g_paranal.run.parvalues))
      ud.replay.pars = g_paranal.run.parvalues(:, 1);
      tstep = (ud.replay.pars(end) - ud.replay.pars(1)) / length(ud.replay.pars);
      ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',size(g_paranal.run.parvalues,1));
      No = ud.replay.opt;
      if No > length(g_grind.paranal.(plotsfield))
         No = length(g_grind.paranal.(plotsfield));
      end

      plt = g_grind.paranal.(plotsfield){No};
      outputlist = i_paranaldialog('outputlist');
      ud.replay.outputtype = lower(outputlist{g_grind.paranal.dlg.outputtype});
      ud.replay.n = max([length(plt.yaxis), length(plt.zaxis), length(plt.xaxis)]);
      ud.replay.ydata = cell(ud.replay.n , 1);
      ud.replay.xdata = cell(ud.replay.n , 1);
      for i = 1:ud.replay.n
         if i <= length(plt.xaxis)
            ud.replay.xdata{i} = real(outfun(struct('fun',plt.xaxis{i},'opts',{{'-p'}})));
         end

         if i <= length(plt.yaxis)
            [ud.replay.ydata{i},ndx] = outfun(struct('fun',plt.yaxis{i},'opts',{{'-p'}},'catfun',ud.replay.outputtype));
            if ~isempty(ndx)
               y = ud.replay.ydata{i};
               %(can have more than one columns
               siz=size(g_paranal.run.t);
               siz(2)=size(y,2);
               ydata=zeros(siz(1)*siz(3),siz(2)) + NaN;
               %only fill the ndx points
               ydata(ndx,:)=real(y);
               %permute back to g_paranal style
               ud.replay.ydata{i} =permute(reshape(transpose(ydata),siz(2),siz(1),siz(3)),[2,1,3]);
            end

         end

         if i <= length(plt.zaxis)
            [ud.replay.zdata{i},ndx] = outfun(struct('fun',plt.zaxis{i},'opts',{{'-p'}},'catfun',ud.replay.outputtype));
            if ~isempty(ndx)
               z = ud.replay.zdata{i};
               %(can have more than one columns
               siz=size(g_paranal.run.t);
               siz(2)=size(z,2);
               zdata=zeros(siz(1)*siz(3),siz(2)) + NaN;
               %only fill the ndx points
               zdata(ndx,:)=real(z);
               %permute back to g_paranal style
               ud.replay.zdata{i} =permute(reshape(transpose(zdata),siz(2),siz(1),siz(3)),[2,1,3]);
            end

            
         end

      end

      if ~ishandle(ud.replay.hwhere)
         oldnext = get(hax, 'NextPlot');
         set(hax,'NextPlot','add');
         if ~isempty(ud.replay.zdata{1})
            ud.replay.hwhere = plot3(hax,ud.replay.xdata{1}(end),ud.replay.ydata{1}(end),ud.replay.zdata{1}(end), 'ro');
         else
            ud.replay.hwhere = plot(hax,ud.replay.xdata{1}(end),ud.replay.ydata{1}(end), 'ro');
         end

         set(ud.replay.hwhere, 'MarkerFaceColor', [1 0 0]);
         set(ud.replay.hwhere,'tag','hwhere');
         i_grindlegend(-1, ud.replay.hwhere);
         set(hax, 'NextPlot', oldnext);
      end

      
   end

   set(hax, 'userdata', ud);
   %   set(hax,'drawmode','fast');
   if isfield(g_paranal, 'run')&&~isempty(g_paranal.run.parvalues)
%      xlim(hax, sort(real([g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)])));
 %     set(hax,'xlimmode','manual');
   end

end


function replayend(hax, closedlg)
try
   if closedlg
      replaycallback(hax, '', 1);
   end

catch
end

if ishandle(hax)
   i_figure(get(hax, 'parent'));
   % set(ax,'xlimmode','auto');
   set(hax,'ylimmode','auto');
   set(hax,'zlimmode','auto');
   ud = get(hax, 'userdata');
   ud.replay.xdata = {};
   ud.replay.ydata = {};
   ud.replay.zdata  = {};
   if ishandle(ud.replay.hwhere)
      delete(ud.replay.hwhere);
      ud.replay.hwhere = -1;
   end

   set(hax, 'userdata', ud);
end

function p = replaycallback(hax, avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&isempty(avar)||(isfield(g_paranal, 'run')&&~isempty(g_paranal.run)&&~isempty(g_paranal.run.parvalues)&&strcmp(avar, g_paranal.run.pars{1}))
   ud = get(hax, 'userdata');
   ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
   [~, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
   t = g_paranal.run.t(ndx1);
   p = g_paranal.run.parvalues(stepndx, 1);
   %  t=g_paranal.run.t(1)+relt*(g_paranal.run.t(end)-g_paranal.run.t(1));
   if isfield(ud, 'replay')
      if ~isfield(ud.replay, 'xdata')
         replaystart(hax);
         ud = get(hax, 'userdata');
      end

      ser = findobj(hax, 'tag','paranal+');
      if g_paranal.run.parvalues(1, 1) > g_paranal.run.parvalues(end, 1)
         ser1= findobj(hax, 'tag', 'paranal-');
         if ~isempty(ser1)
            ser = ser1;
         end

      end

      ix = [];
      iy = [];
      iz = [];
      for i = 1:min(length(ser), ud.replay.n)
         if i<=length(ud.replay.xdata)&&~isempty(ud.replay.xdata{i})
            ix = i;
         end

         if i<=length(ud.replay.ydata)&&~isempty(ud.replay.ydata{i})
            iy = i;
         end

         if i<=length(ud.replay.zdata)&&~isempty(ud.replay.zdata{i})
            iz = i;
         end

         if ~isempty(g_paranal.run.t)
            ndx=g_paranal.run.t <= t;
            if ~isempty(iz)&&~isempty(ud.replay.zdata)&&~isempty(ud.replay.zdata{iz})
               set(ser(end-i+1), 'XData', ud.replay.xdata{ix}(ndx), 'YData', ud.replay.ydata{iy}(ndx), ...
                  'ZData', ud.replay.zdata{iz}(ndx));
               if ishandle(ud.replay.hwhere)
                  set(ud.replay.hwhere,'Xdata',ud.replay.xdata{ix}(ndx1),'Ydata',ud.replay.ydata{ix}(ndx1),...
                     'Zdata', ud.replay.zdata{ix}(ndx1));
               end

               
            else
               set(ser(end-i+1), 'XData', ud.replay.xdata{ix}(ndx), 'YData', ud.replay.ydata{iy}(ndx));
               if ishandle(ud.replay.hwhere)
                  set(ud.replay.hwhere,'Xdata',ud.replay.xdata{ix}(ndx1),'Ydata',ud.replay.ydata{ix}(ndx1));
               end

            end

         end

      end

   end

end


function [X, Y, Z] = makegrids(p, p1, Y1)
x = unique(p);
y = unique(p1);
[X, Y] = meshgrid(x, y);
ndx = ~(isnan(p) | isnan(p1) | isnan(Y1));
warning('off','MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
F = TriScatteredInterp(p(ndx), p1(ndx), Y1(ndx));
warning('on','MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
Z = F(X, Y);

function [ps, funx, funy] = addbreaks(g_paranal, ps, funx, funy)
if size(g_paranal.run.parvalues, 2) == 2
   diff1 = g_paranal.run.parvalues(2, :) - g_paranal.run.parvalues(1, :);
   if (diff1(1)==0)&&(diff1(2)~=0)
      ddif=find(diff(ps) ~= 0);
   elseif (diff1(2)==0)&&(diff1(1)~=0)
      ddif=find(diff(funx) ~= 0);
   else
      return;
   end

   for i = length(ddif):-1:1
      ps = [ps(1:ddif(i)); NaN; ps(ddif(i) + 1:end)];
      funx = [funx(1:ddif(i)); NaN; funx(ddif(i) + 1:end)];
      funy = [funy(1:ddif(i)); NaN; funy(ddif(i) + 1:end)];
   end

end


