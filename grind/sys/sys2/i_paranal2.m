% internal function of parameter analyser
function i_paranal2(plotprev, conteq_eq)
if nargin < 2
   conteq_eq = 0;
end

%parname, start, nend, nsteps, nstabilizing, ndays, lines, plotprev, outputtype, disturb)
%   i_paranal(answer.par{1}, answer.start(1), answer.nend(1), answer.steps(1), answer.stabil, ...
%      answer.writing,answer.lines,plotprev, answer.outputtype, answer.disturbance);

global t g_grind g_Y g_t g_paranal;
parname = g_grind.paranal.dlg.par{1};
start = g_grind.paranal.dlg.start(1);
nend = g_grind.paranal.dlg.nend(1);
nsteps = g_grind.paranal.dlg.steps(1);
ppars=linspace(start,nend,nsteps+1); %more precise than summing per step
nstabilizing = g_grind.paranal.dlg.stabil;
ndays = g_grind.paranal.dlg.writing;
lines = g_grind.paranal.dlg.lines;
outputtype = g_grind.paranal.dlg.outputtype;
silent = 0;
if isfield(g_grind.paranal, 'silent') &&( g_grind.paranal.silent==1)
   g_grind.paranal.silent = 0;
   silent = 1;
end

%disturb=g_grind.paranal.dlg.disturbance;

if ~isfield(g_grind, 'paranal') || ~isfield(g_grind.paranal, 'defaultplots') || g_grind.paranal.defaultplots
   paranal('-defaults');
end

parlim = [start, nend];
if start > nend
   parlim = [nend, start];
end

for i = 1:length(g_grind.paranal.plots)
   if strcmp(g_grind.paranal.plots{i}.xaxis{1}, '<param1>')
      g_grind.paranal.plots{i}.xlim = parlim;
   end

   if strcmp(g_grind.paranal.plots{i}.yaxis{1}, '<param1>')
      g_grind.paranal.plots{i}.ylim = parlim;
   end

   if strcmp(g_grind.paranal.plots{i}.zaxis{1}, '<param1>')
      g_grind.paranal.plots{i}.zlim = parlim;
   end

end

if lines > 1
   lines = 0;
end


outputlist = i_paranaldialog('outputlist');
outputtype = lower(outputlist{outputtype});

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
% iY = i_varno(g_grind.paranal.plots{1}.yaxis{1});
% iZ = i_varno(g_grind.paranal.plots{1}.zaxis{1});
if silent
   h = [];
else
   parfig = i_figno('paranal');
   [h, isnew] = i_makefig('paranal');  
   set(h, 'Name', 'Parameter analysis');
   if  ~isempty(g_grind.paranal.plots{1}.zaxis{1}) && isnew
      set(gca, 'View', [8.5 10])
      hold on;
   end

   if ~ishold
        delete(get(gca,'children'))
        set(gca,'userdata',[])
        legend off
   end

   i_odespeed(0, 0, 'init');
  % plotedit off;
   ch=findobj(gca,'type','line');
   g_grind.paranal.pen = i_nextpen([]);
   for i = 1:length(ch)
      g_grind.paranal.pen = i_nextpen(g_grind.paranal.pen);
   end

   if lines
      ppen = g_grind.paranal.pen.pen2;
   else
      ppen = '.';
   end
end

oldpar = evalin('base', char(parname));
oldtstep = g_grind.tstep;
% funx = [];
% funy = [];
% ps = [];
%try
hasperm = ~isempty(g_grind.permanent);
if plotprev ~= 1
   g_paranal.Y = [];
   g_paranal.t = [];
   g_paranal.p = [];
   %   pperm=[];
   if hasperm
      g_paranal.perm = [];
   end

   g_paranal.parname = parname;
   %  wb=waitbar(0,'Calculating...');
   %try
   %   oldi = 1;
   %  newi = 0;
   t1 = t;
   N0 = i_initvar;
   OldN0 = N0;
   nstat = size(N0, 1);
   if isempty(g_grind.paranal.plots)
      xaxis1 = '';
      yaxis1 = '';
      zaxis1 = '';
   else
      xaxis1 = strrep(g_grind.paranal.plots{1}.xaxis{1}, '<param1>', parname);
      if ~isempty(g_grind.paranal.plots{1}.yaxis{1})
         yaxis1 = strrep(g_grind.paranal.plots{1}.yaxis{1}, '<param1>', parname);
      else
         yaxis1 = '';
      end

      if ~isempty(g_grind.paranal.plots{1}.zaxis{1})
         zaxis1 = strrep(g_grind.paranal.plots{1}.zaxis{1}, '<param1>', parname);
      else
         zaxis1 = '';
      end

   end

   for i = 1:nsteps + 1
      %  waitbar(i/nsteps,wb);
      g_grind.solver.opt.OutputFcn = [];
      ppar=ppars(i);
      multassignin('base', char(parname), ppar);
      if ~g_grind.solver.isdiffer
         g_grind.tstep = 2;
      else
         g_grind.tstep = NaN;
      end

      if nstabilizing > 0
         if i == 1 %extra initial stabilizing
            i_ru(t1, 2 * nstabilizing, N0, 0);
         else
            i_ru(t1, nstabilizing, N0, 0);
         end
         N0 = transpose(g_Y(size(g_Y, 1), :));
         drawnow;
      end
      g_grind.tstep = oldtstep;
      if ndays == 0
         g_t = i - 1;
         g_Y = transpose(N0);
      else
         if ~isempty(g_t)
            t1 = g_t(size(g_t, 1));
         else
            t1 = 0;
         end
         i_ru(t1, ndays, N0, 0);
         drawnow;
         t1 = g_t(size(g_t, 1));
      end

      p = repmat(ppar,size(g_t, 1), 1);
      N0 = transpose(g_Y(size(g_Y, 1), :));
      %little invasion to avoid hanging in a trivial equilibrium (tricky
      %if the model can have negative numbers
      if ~g_grind.solver.isinteger
      for j = 1:nstat
         epsil = 0.001;
         if isnan(N0(j)) || ((N0(j) < epsil) && (N0(j) > 0))
            N0(j) = epsil;
         elseif (N0(j) >  - epsil)  && (N0(j) < 0)
            N0(j) =    -epsil;
         elseif N0(j) == 0
            N0(j) = sign(OldN0(j)) * epsil;
         else %very small positive disturbance
            N0(j) = N0(j) + epsil * 0.001 * sign(N0(j));
         end
      end
      end

%      ppar = ppar + (nend - start) / nsteps;
      g_paranal.Y = [g_paranal.Y; g_Y];
      g_paranal.t = [g_paranal.t; g_t];
      g_paranal.p = [g_paranal.p; p];
      if hasperm
         pperm = defpermanent('-g', []);
         g_paranal.perm = [g_paranal.perm; pperm];
      end

      %[xdata, ydata, zdata] = i_paranalfun(1,iY, iZ, outputtype, parname, g_t, g_Y, p, pperm, disterror);
      xdata = outfun(xaxis1);
      [ydata,ndx] = outfun(yaxis1, '-n',outputtype);
      if ~isempty(ndx) && (length(ndx)~=length(xdata))
         xdata = xdata(ndx);
      end

      if ~isempty(zaxis1)
         zdata = outfun(zaxis1, '-n',outputtype);
         if length(ydata) ~= length(zdata)
            zdata = outfun(zaxis1);
            zdata = zdata(ndx, :);
         end

      else
         zdata = [];
      end

      if ~silent && g_grind.drawnow && (gcf == parfig)
         ud = get(parfig, 'userdata');
         if ~isempty(ud) && isfield(ud, 'stop') && ud.stop
            break;
         end

         if i == 1
            if isempty(zdata)
               [h] = makeplot(1,ppen, xdata(:,1),ydata(:,1),zdata);
            else
               [h] = makeplot(1,ppen, xdata(:,1),ydata(:,1),zdata(:,1));
            end
            hax = gca;
            if ~isempty(g_grind.paranal.plots{1}.xlim)
               set(hax, 'Xlim', g_grind.paranal.plots{1}.xlim);
            end

            if ~isempty(g_grind.paranal.plots{1}.ylim)
               set(hax, 'Ylim', g_grind.paranal.plots{1}.ylim);
            end

            if ~isempty(g_grind.paranal.plots{1}.zlim)
               set(hax, 'Zlim', g_grind.paranal.plots{1}.zlim);
            end

            set(h, 'Color', g_grind.paranal.pen.drawcolor);
            xlabel(mydisptext(g_grind.paranal.plots{1}.xaxis{1}, parname));
            ylabel(mydisptext(g_grind.paranal.plots{1}.yaxis{1}, parname));
            if ~isempty(g_grind.paranal.plots{1}.zaxis{1})
               zlabel(mydisptext(g_grind.paranal.plots{1}.zaxis{1}, parname));
            end

         else
            if isempty(zdata)
               set(h, 'xdata', xdata(:,1), 'Ydata', ydata(:,1),  'Color', g_grind.paranal.pen.drawcolor);
            else
               set(h, 'xdata', xdata(:,1), 'Ydata', ydata(:,1), 'zdata', zdata(:,1), 'Color', g_grind.paranal.pen.drawcolor);
            end

            drawnow;
         end

      end

   end

   multassignin('base', char(parname), oldpar);
   if ishandle(h)
      delete(h);
   end

end

%gather all data
%[ps, funx,funy, iY] = i_paranalfun(1,iY, iZ, outputtype, parname, g_paranal.t, g_paranal.Y, g_paranal.p, disterror);
%h = i_makefig('paranal');
if ~silent
   i_odespeed(1, 1, 'done');
   try
      for i = 1:length(g_grind.paranal.plots)
         nx = length(g_grind.paranal.plots{i}.xaxis);
         ny = length(g_grind.paranal.plots{i}.yaxis);
         nz = length(g_grind.paranal.plots{i}.zaxis);
         [h,isnew] = i_makefig('paranal', i - 1);
         set(gca, 'YLimMode', 'auto');
         set(gca, 'ZLimMode', 'auto');
         set(h, 'userdata', []);
         if isnew
             hold('on');
         end

         for k = 1:max([nx, ny, nz])
            ps = outfun(g_grind.paranal.plots{i}.xaxis{min(nx, k)}, '-p');
            eq1=g_grind.paranal.plots{i}.yaxis{min(ny,k)};
            [funx,ndx] = outfun(eq1,'-p', outputtype);
            if ~isempty(ndx) && (length(ndx)~=length(ps))
               ps = ps(ndx);
            end

            if ~isempty(g_grind.paranal.plots{i}.zaxis{min(nz, k)})
               funy = outfun(g_grind.paranal.plots{i}.zaxis{min(nz,k)},'-p', outputtype);
               if length(funy) ~= length(funx)
                  funy = outfun(g_grind.paranal.plots{i}.zaxis{min(nz, k)}, '-p');
                  funy = funy(ndx, :);
               end

            else
               funy = [];
            end

            %      iY = i_varno(g_grind.paranal.plots{i}.yaxis{1});
            %      iZ = i_varno(g_grind.paranal.plots{i}.zaxis{1});
            %      [ps, funx,funy, iY] = i_paranalfun(i, iY, iZ, outputtype, parname, g_paranal.t, g_paranal.Y, g_paranal.p, iif(hasperm,g_paranal.perm,[]), disterror);
            g_grind.paranal.pen = nextpen(g_grind.paranal.pen);
            hplot = makeplot(i, ppen, ps, funx, funy);
         %   set(hplot,'EraseMode', 'normal')
            addreplaydata(gca,i);
            nser = max(size(funx, 2), size(funy, 2));
            increasingpar = sign(g_grind.paranal.dlg.nend(1) - g_grind.paranal.dlg.start(1));
            if conteq_eq > 0
               set(hplot,'tag','conteq');
               increasingpar = conteq_eq;
               if conteq_eq==2
                  sers=findobj(gca,'type','line');
                  if ~isempty(sers)
                     set(sers(1:nser),'LineStyle',':')
                  end

               end

            end

            if nser == 1
               if isempty(g_grind.paranal.plots{i}.zaxis{1})
                  i_grindlegend(2,hplot,mydisptext(sprintf('%s', eq1)),{'g_increasingpar',increasingpar});
               else
                  i_grindlegend(2,hplot,' ',{'g_increasingpar',increasingpar});
               end

            else
               for j = 1:nser
                  if isempty(g_grind.paranal.plots{i}.zaxis{1})&&j<=ny
                     i_grindlegend(2,hplot,mydisptext(sprintf('%s(%d)', eq1,j)),{'g_increasingpar',increasingpar});
                  else
                     i_grindlegend(2,hplot,' ',{'g_increasingpar',increasingpar});
                  end

               end

            end

         end

         set(h, 'Name', sprintf('Parameter analysis - %d',i));
         if ~isempty(g_grind.paranal.plots{i}.zaxis{1})
            set(gca, 'View', [8.5 10])
         end
         if ~isempty(g_grind.paranal.plots{i}.xlim)
            xlim(g_grind.paranal.plots{i}.xlim);
         end

         xlabel(mydisptext(g_grind.paranal.plots{i}.xaxis{1}, parname));
         if ~strcmp(outputtype, 'unchanged')
            ylabel(mydisptext(sprintf('%s_{%s}', g_grind.paranal.plots{i}.yaxis{1}, outputtype),parname));
            zlabel(mydisptext(sprintf('%s_{%s}',g_grind.paranal.plots{i}.zaxis{1}, outputtype),parname));
         else
            ylabel(mydisptext(g_grind.paranal.plots{i}.yaxis{1}, parname));
            zlabel(mydisptext(g_grind.paranal.plots{i}.zaxis{1}, parname));
         end

         set(findobj(gcf,'type','line'),'linewidth',g_grind.paranal.pen.linewidth);
         i_plotdefaults(gcf);
     %    plotedit on;
         i_grindlegend(11);
         %   if ~oldhold
         %      hold off;
         %   end

      end

      g_Y = [];
      g_t = [];
   catch err
      %    err=lasterror;
      %in case of error restore the old parameter (if functions
      %are used, the parameter has become a matrix)
      g_grind.tstep = oldtstep;
      multassignin('base', char(parname), oldpar);
      rethrow(err);
   end

end

g_grind.tstep = oldtstep;
multassignin('base', char(parname), oldpar);

function s = mydisptext(s, parname)
if strcmp(s, '<param1>')
   s = parname;
end

s = i_disptext(s);

function h = makeplot(iplot,ppen, ps, funx, funy)
global g_grind
maxn = max([size(ps, 2), size(funx, 2), size(funy, 2)]);
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
         'MarkerSize', 3,  'EraseMode', 'none',  'Color', apen.color2);
      if ps(1)<ps(end)
          set(h,'tag','paranal+');
      else
          set(h,'tag','paranal-');
      end

      if ~isoctave&&verLessThan('matlab','8.4.0')
          set(gca, 'drawmode','fast');
      else
          set(gca,'SortMethod','depth');
      end

      apen = nextpen(apen);
      if ppen ~= '.'
         ppen = apen.pen2;
      end

      box on;
   end

elseif ~isempty(funx)
   for i = 1:maxn
      h = plot(ps(:,min(size(ps,2),i)), funx(:,min(size(funx,2),i)), ppen, ...
         'MarkerSize', 3, 'EraseMode', 'none', 'Color', apen.color2);
      if ps(1)<ps(end)
          set(h,'tag','paranal+');
      else
          set(h,'tag','paranal-');
      end

      apen = nextpen(apen);
      if ppen ~= '.'
         ppen = apen.pen2;
      end

   end

else
   warning('GRIND:paranal','Cannot plot empty series');
   return;
end

if maxn > 1
   if length(g_grind.paranal.plots{iplot}.yaxis) == 1
      s1 = sprintf('%s(%%d)\n', g_grind.paranal.plots{iplot}.yaxis{:});
      s = str2cell(sprintf(s1, (1:maxn)));
      s = i_addlegend(s);
      legend(s)
   else
      s = {g_grind.paranal.plots{iplot}.yaxis};
      s = i_addlegend(s);
      legend(s);
   end

end

if ~oldhold
   hold off;
end


function addreplaydata(h,No)
global g_paranal;
ud=get(h,'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
ud.replay.onturn = @i_replayparanalturn;
ud.replay.pars=unique(g_paranal.p);
h1=findobj(h,'tag','hwhere');
if ~isempty(h1)
    delete(h1);
end

ud.replay.hwhere = -1;
tstep=(ud.replay.pars(end)-ud.replay.pars(1))/length(ud.replay.pars);
ud.replay.settings=struct('tvar',g_paranal.parname,'tlim',[g_paranal.p(1) g_paranal.p(end)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',length(g_paranal.p));
ud.replay.opt = No;
set(h,'userdata', ud);




function replaystart(hax)
global g_grind g_paranal;
if ishandle(hax)
   i_figure(get(hax,'parent'));
   set(hax,'ylimmode','manual');
   set(hax,'zlimmode','manual'); 
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')&&(isfield(g_paranal,'p')&&~isempty(g_paranal.p))
      ud.replay.pars=unique(g_paranal.p);
      tstep=(ud.replay.pars(end)-ud.replay.pars(1))/length(ud.replay.pars);
      ud.replay.settings=struct('tvar',g_paranal.parname,'tlim',[g_paranal.p(1) g_paranal.p(end)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',length(g_paranal.p));
      No = ud.replay.opt;
      if No>length(g_grind.paranal.plots)
          No=length(g_grind.paranal.plots);
      end

      plt=g_grind.paranal.plots{No};
      ud.replay.n = max([length(plt.yaxis), length(plt.zaxis), length(plt.xaxis)]);
      ud.replay.ydata = cell(ud.replay.n ,1);
      ud.replay.xdata = cell(ud.replay.n ,1);
      for i = 1:ud.replay.n
         if i<=length(plt.xaxis)
             ud.replay.xdata{i}=outfun(plt.xaxis{i},'-p');
         end

         if i<=length(plt.yaxis)
             ud.replay.ydata{i}=outfun(plt.yaxis{i},'-p');
         end

         if i<=length(plt.zaxis)
             ud.replay.zdata{i}=outfun(plt.zaxis{i},'-p');
         end
         if any(~isreal(ud.replay.xdata{i}))||any(~isreal(ud.replay.ydata{i}))||any(~isreal(ud.replay.zdata{i}))
             ud.replay.xdata{i}=real(ud.replay.xdata{i});
             ud.replay.ydata{i}=real(ud.replay.ydata{i});
             ud.replay.zdata{i}=real(ud.replay.zdata{i});
             warning('grind:paranal2:complex','Ignored imaginary parts of complex values')
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
         set(hax,'NextPlot',oldnext);
      end


   end

   set(hax, 'userdata', ud);
%   set(hax,'drawmode','fast');
   if isfield(g_paranal,'p')&&~isempty(g_paranal.p)
      xlim(hax,sort([g_paranal.p(1) g_paranal.p(end)]));
      set(hax,'xlimmode','manual');
   end

end


function replayend(hax,closedlg)
try
    if closedlg
       replaycallback(hax,'',1);
    end

catch
end

if ishandle(hax)
   i_figure(get(hax,'parent'));
  % set(ax,'xlimmode','auto');
   set(hax,'ylimmode','auto');
   set(hax,'zlimmode','auto'); 
   ud = get(hax, 'userdata');
   ud.replay.xdata = {};
   ud.replay.ydata = {};
   ud.replay.zdata ={};
   if ishandle(ud.replay.hwhere)
       delete(ud.replay.hwhere);
      ud.replay.hwhere = -1;
   end

   set(hax, 'userdata', ud);
end

function p=replaycallback(hax, avar, relt)
global g_paranal;
p=[];
if ishandle(hax)&&isempty(avar)||(isfield(g_paranal,'p')&&~isempty(g_paranal.p)&&strcmp(avar,g_paranal.parname))
   ud = get(hax, 'userdata');
    ndx1=floor(relt*(length(g_paranal.t)-1))+1;
    t=g_paranal.t(ndx1);
    p=g_paranal.p(ndx1);
 %  t=g_paranal.t(1)+relt*(g_paranal.t(end)-g_paranal.t(1));
   if isfield(ud, 'replay')
      if ~isfield(ud.replay, 'xdata')
         replaystart(hax);
         ud = get(hax, 'userdata');
      end

      ser = findobj(hax, 'tag','paranal+');
      if g_paranal.p(1)>g_paranal.p(end)
          ser1= findobj(hax, 'tag', 'paranal-');
          if ~isempty(ser1)
              ser=ser1;
          end

      end

      ix=[];
      iy=[];
      iz=[];
      for i = 1:min(length(ser),ud.replay.n)
         if i<=length(ud.replay.xdata)&&~isempty(ud.replay.xdata{i})
             ix=i;
         end

         if i<=length(ud.replay.ydata)&&~isempty(ud.replay.ydata{i})
             iy=i;
         end

         if i<=length(ud.replay.zdata)&&~isempty(ud.replay.zdata{i})
             iz=i;
         end

         if ~isempty(g_paranal.t)
            ndx=g_paranal.t<=t;
            if ~isempty(iz)&&~isempty(ud.replay.zdata)&&~isempty(ud.replay.zdata{iz})
              set(ser(i), 'XData', ud.replay.xdata{ix}(ndx), 'YData', ud.replay.ydata{iy}(ndx), ...
               'ZData', ud.replay.zdata{iz}(ndx));
              if ishandle(ud.replay.hwhere)
                 set(ud.replay.hwhere,'Xdata',ud.replay.xdata{ix}(ndx1),'Ydata',ud.replay.ydata{ix}(ndx1),...
                     'Zdata',ud.replay.zdata{ix}(ndx1));
              end


            else
               set(ser(i), 'XData', ud.replay.xdata{ix}(ndx), 'YData', ud.replay.ydata{iy}(ndx));
              if ishandle(ud.replay.hwhere)
                 set(ud.replay.hwhere,'Xdata',ud.replay.xdata{ix}(ndx1),'Ydata',ud.replay.ydata{ix}(ndx1));
              end

            end

         end

      end

   end

end


