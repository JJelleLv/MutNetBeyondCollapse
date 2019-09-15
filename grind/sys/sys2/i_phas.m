% i_phas initial function to create phase plane
function [] = i_phas(curplot, newrun)
global g_grind g_t g_Y;
i = i_figno('phase1');
done = 0;
if curplot == i
   done = 1;
   [hfig, new] = i_makefig('phase1');
   hax=get(hfig,'currentaxes');
   if new
      set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
      set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
      avar=g_grind.xaxis.var;
   else
      ud=get(hax,'userdata');
      if isfield(ud,'var')
         avar=ud.var;
      else
         avar=g_grind.xaxis.var;
      end

   end

   iX=i_getno(avar);
      if ~iX.isvar
         avar=g_grind.yaxis.var;
         iX=i_getno(avar);
         if ~iX.isvar
           avar=i_statevars_names(1);
        end

      end

   oldhold = ishold;
   hold on;
   if g_grind.solver.isdiffer
      if new
         close(hfig);
         itermap
         hfig=gcf;
         hold on;
      end

      set(hfig, 'Name', 'Iteration map');
      %difference equation: cobwebbing
      if size(g_Y, 1) == 1
         h=plot(get(hfig,'CurrentAxes'),g_Y, g_Y, g_grind.pen.pen, 'Color', g_grind.pen.color, 'linewidth',g_grind.pen.linewidth...
             , 'markersize',g_grind.pen.markersize, 'markerfacecolor',g_grind.pen.markerfacecolor); %plot equilibrium only
         i_grindlegend(-1, h); % Exclude line from legend
      else
         %cobwebbing plotting method for difference equations
         h = cobwebbing(hfig,outfun(avar), g_grind.drawnow);
         addreplaydata(hfig, h, 4,{avar})
         i_grindlegend(-1, h); % Exclude line from legend
      end

   else
      %1D differential equation: draw on x axis
      if new
         close(hfig);
         plotdiff;
         hfig=gcf;
         hold on;
         set(hfig, 'Name','Plot of 1D differential equation');
      end

      YY = outfun(avar);
      if g_grind.drawnow
         i_drawslowly(hfig,YY, zeros(size(YY)), 1);
      end

      h=plot(get(hfig,'CurrentAxes'),YY, zeros(size(YY)), g_grind.pen.pen, 'Color', g_grind.pen.color2,'linewidth',g_grind.pen.linewidth,...
          'markersize',g_grind.pen.markersize,'markerfacecolor',g_grind.pen.markerfacecolor);
      addreplaydata(hfig, h, 1, {avar})
      i_grindlegend(-1, h); % Exclude line from legend
   end

   xlabel(get(hfig,'CurrentAxes'),i_disptext(avar));
   if ~oldhold
      hold off;
   end

elseif newrun
   if ishandle(i)
      set(i, 'visible', 'off');
   end

end

i = i_figno('phase1a'); %plotreldiff
if curplot == i
    done = 1;
    ud=get(gca,'userdata');
    avar=ud.meta.xname;%g_grind.xaxis.var;
    YY = outfun(avar);
    if g_grind.drawnow
        i_drawslowly(gcf,YY(:,1), zeros(1, size(YY, 1)), 1);
    end
    
    h=plot(YY(:,1), zeros(1, size(YY, 1)), g_grind.pen.pen, 'Color', g_grind.pen.color2,'linewidth',g_grind.pen.linewidth, ...
        'markersize',g_grind.pen.markersize, 'markerfacecolor',g_grind.pen.markerfacecolor);
    addreplaydata(gcf, h, 1)
    i_grindlegend(-1, h); % Exclude line from legend
end

i = i_figno('phase2');
if curplot == i
   done = 1;
   [hfig, new] = i_makefig('phase2');
   if new
      set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
      set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
   end

   oldhold = ishold;
   hold on;
   if isempty(g_grind.yaxis.var)
      close(hfig);
      i_phas(i_figno('phase1'), newrun);
      return;
   else
      h=plot(outfun(g_grind.xaxis.var), outfun(g_grind.yaxis.var),g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth,...
          'markersize',g_grind.pen.markersize, 'markerfacecolor',g_grind.pen.markerfacecolor);
      addreplaydata(hfig, h, 2)
      i_grindlegend(-1, h); % Exclude line from legend
   end

   i_plotdefaults(hfig);
   set(hfig, 'Name', 'Phase plane');
   xlabel(i_disptext(g_grind.xaxis.var));
   ylabel(i_disptext(g_grind.yaxis.var));
   if ~oldhold
      hold off;
   end

elseif newrun
   if ishandle(i)
      set(i, 'visible', 'off');
   end

end

if ~isempty(g_grind.zaxis.var)
   i = i_figno('phase3');
   if curplot == i
      done = 1;
      [hfig, fignew] = i_makefig('phase3');
      if fignew
         set(hfig, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
         set(hfig, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
      end

      set(hfig, 'Name', '3D phase space');
      oldhold = ishold;
      hold on;
      %plotedit off;
      h =  plot3(outfun(g_grind.xaxis.var), outfun(g_grind.yaxis.var), outfun(g_grind.zaxis.var) ...
         , g_grind.pen.pen, 'Color', g_grind.pen.color,'markersize',g_grind.pen.markersize, 'markerfacecolor',g_grind.pen.markerfacecolor); %'Erase', 'none', 
      if ~isoctave&&verLessThan('matlab','8.4.0')
          set(gca, 'drawmode','fast');
      else
          set(gca,'SortMethod','depth');
      end

      i_plotdefaults(hfig);
      addreplaydata(hfig, h, 3)
      box on;
      xlabel(i_disptext(char(g_grind.xaxis.var)));
      ylabel(i_disptext(char(g_grind.yaxis.var)));
      zlabel(i_disptext(char(g_grind.zaxis.var)));
      if ~oldhold
         hold off;
      end

      if fignew
         set(gca, 'View', [322.5, 30]);
      end

   elseif newrun
      if ishandle(i)
         set(i, 'visible', 'off');
      end

   end

end


i = i_figno('poinmap');
if curplot == i
   i_figure(i,g_grind.figopts{:});
   ud = get(i, 'userdata');
   if ~isempty(ud)
      poincaremap(ud.avar1, ud.avar, ud.avalue, ud.increasing);
   end

end

i = i_figno('poinsec');
if curplot == i
   i_figure(i,g_grind.figopts{:});
   ud = get(i, 'userdata');
   if ~isempty(ud)
      poincaresect(ud.avar, ud.avalue, ud.increasing);
   end

end

i = i_figno('lorenzmap');
if curplot == i
   i_figure(i,g_grind.figopts{:});
   ud = get(i, 'userdata');
   if ~isempty(ud)
      lorenzmap(ud);
   end

end

i = i_figno('torus');
if curplot == i
   done = 1;
   i_figure(i,g_grind.figopts{:});
   ud = get(i, 'userdata');
   if ~isempty(ud)
      torus(ud.period, ud.increment);
   else
      torus
   end

end

i = i_figno('time');
if curplot == i
   done = 1;
   time('-n');
end

i = i_figno('dirfield');
if curplot == i
   done = 1;
   hfig = i_figure(i,g_grind.figopts{:});
   oldhold=ishold;
   hold on;
   avar=g_grind.xaxis.var;
   iY = i_getno(g_grind.xaxis.var);
   if ~iY.isvar
      avar=g_grind.yaxis.var;
   end

   h=plot(g_t, outfun(avar), g_grind.pen.pen, 'Color', g_grind.pen.color, 'linewidth',g_grind.pen.linewidth, ...
       'markersize',g_grind.pen.markersize, 'markerfacecolor',g_grind.pen.markerfacecolor);
   addreplaydata(hfig, h, 5)
   if ~oldhold
      hold off;
   end

end

% i = i_figno('potent3');
% if curplot == i
%    done = 1;
%    hfig = figure(i);
%    hold on;
%    h=plot3(outfun(g_grind.xaxis.var), outfun(g_grind.yaxis.var), zeros(size(g_Y, 1)), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
%    addreplaydata(hfig, h, 2);
%    if verLessThan('matlab','8.4.0')
%        set(gca, 'drawmode','fast');
%    else
%        set(gca,'SortMethod','depth');
%    end

%    hold off;
% end

% i = i_figno('potent2');
% if curplot == i
%    done = 1;
%    hfig = figure(i);
%    hold on;
%    h=plot(outfun(g_grind.xaxis.var), outfun(g_grind.yaxis.var), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
%    addreplaydata(hfig, h, 2);
%    hold off;
% end

% i = i_figno('potent1');
% if curplot == i
%    done = 1;
%    hfig = figure(i);
%    hold on;
%    yy = outfun(g_grind.xaxis.var);
%    h=plot(yy, 0, g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
%    plot(yy(end), 0, 'o', 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth, 'markersize',g_grind.pen.markersize);
%    addreplaydata(hfig, h, 1);
%    hold off;
% end

i = i_figno('growths');
if ishandle(i + 1)
   done = 1;
   for j = 1:g_grind.statevars.dim
      if curplot == i + j
         hfig = i_figure(i + j,g_grind.figopts{:});
         hold on;
         if ~isoctave&&verLessThan('matlab','8.4.0')
             set(gca,'DrawMode', 'fast');
         else
             set(gca,'SortMethod','depth');
         end
         h=plot(outfun(g_grind.xaxis.var),outfun(g_grind.yaxis.var), g_grind.pen.pen, 'Color', g_grind.pen.color,'linewidth',g_grind.pen.linewidth,...
         'markerfacecolor',g_grind.pen.markerfacecolor,'markersize',g_grind.pen.markersize);
         addreplaydata(hfig, h, 2)
         hold off;
      end

   end

end

% i = i_figno('varcontour');
% if ishandle(i + 1)
%    done = 1;
%    varcontour;
% end


if ~done
   if isempty(g_grind.zaxis.var)
      i_phas(i_figno('phase2'), newrun)
   elseif isempty(g_grind.yaxis.var)
      i_phas(i_figno('phase1'), newrun)
   else
      i_phas(i_figno('phase3'), newrun)
   end

end

if g_grind.pen.cycle
   nextpen;
end


function addreplaydata(H, HLine, atype, vars)
global g_t t g_grind;
if nargin<4
  vars= {g_grind.xaxis.var,g_grind.yaxis.var,g_grind.zaxis.var};
end

hax=findobj(H,'type','axes');
tags = get(hax, 'tag');
hax = hax(~strcmp(tags, 'legend'));
ud = get(hax, 'userdata');
ud.replay.vars=vars;
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
ch=findobj(hax,'tag','HLine');
if ~isempty(ch)
    set(ch,'tag','');
end

ud.replay.HLine = HLine;
set(HLine,'tag','HLine');
h1=findobj(hax,'tag','hwhere');
if ~isempty(h1)
   delete(h1);
end

ud.replay.hwhere = -1;
ud.replay.type = atype;
%  1 = 1dim 2=2dim 3=3dim 4=cobweb 5=2dim time
if isempty(g_t)
   ud.replay.settings=struct('tvar','t','tlim',[t g_grind.ndays],'numt',iif(~isnan(g_grind.tstep),g_grind.tstep,g_grind.ndays));
else
   ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
end

set(hax,'userdata', ud);

function h = cobwebbing(hfig,g_Y, drawnow)
global g_grind;
cobweb = i_getcobweb(g_Y);
if drawnow
   i_drawslowly(hfig,cobweb(:, 1), cobweb(:, 2));
end

h=plot(get(hfig,'CurrentAxes'),cobweb(:, 1), cobweb(:, 2), g_grind.pen.pen, 'Color', g_grind.pen.color2,'linewidth',g_grind.pen.linewidth, ...
    'markersize',g_grind.pen.markersize, 'markerfacecolor',g_grind.pen.markerfacecolor);


function replaystart(hax)
global g_grind g_t;
if ishandle(hax)
   i_figure(get(hax, 'parent'));
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
     if ishandle(ud.replay.HLine(1))
          ax1 = get(ud.replay.HLine(1), 'parent');
          if ax1~=hax
              ch=findobj(hax,'tag','HLine');
              ud.replay.HLine=ch;
          end

      end

      %  hax = get(ud.replay.HLine(1), 'parent');
      oldnext = get(hax, 'NextPlot');
      set(hax,'NextPlot','add');
      ud.replay.tdata = outfun('t');
      if isempty(g_t)
         ud.replay.settings=struct('tvar','t','tlim',[0 g_grind.ndays],'numt',1000);
      else
        ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
        if g_grind.solver.isdiffer&&g_grind.statevars.dim==1
            ud.replay.settings.numt = ud.replay.settings.numt * 2 - 1;
         end

         if ud.replay.type == 4
            cobweb = i_getcobweb(outfun(ud.replay.vars{1}));
            tt = [ud.replay.tdata'-0.5;ud.replay.tdata'];
            tt = tt(:);
            ud.replay.tdata = tt(1:size(cobweb, 1));
            ud.replay.xdata = cobweb(:, 1);
            ud.replay.ydata = cobweb(:, 2);
         else
            ud.replay.xdata = outfun(ud.replay.vars{1});
            if ismember(ud.replay.type, [2, 3])
               ud.replay.ydata = outfun(g_grind.yaxis.var);
               if ud.replay.type == 3
                  ud.replay.zdata = outfun(g_grind.zaxis.var);
               end

            end
         end
          if isfield(ud.replay,'xdata')&&any(~isreal(ud.replay.xdata(:)))
              warning('grind:phase','Imaginary parts of complex values ignored');
              ud.replay.xdata=real(ud.replay.xdata);
          end
          if isfield(ud.replay,'ydata')&&any(~isreal(ud.replay.ydata(:)))
              warning('grind:phase','Imaginary parts of complex values ignored');
              ud.replay.ydata=real(ud.replay.ydata);
          end
          if isfield(ud.replay,'zdata')&&any(~isreal(ud.replay.zdata(:)))
              warning('grind:phase','Imaginary parts of complex values ignored');
              ud.replay.zdata=real(ud.replay.zdata);
          end
            if ~ishandle(ud.replay.hwhere)&&~isempty(ud.replay.xdata)
               switch ud.replay.type
                case 1
                  ud.replay.hwhere = plot(ud.replay.xdata(1), 0, 'ro');
                case {2,4}
                  if ~isempty(ud.replay.ydata)
                     ud.replay.hwhere = plot(ud.replay.xdata(1), ud.replay.ydata(1), 'ro');
                  end

                case 3
                  if ~isempty(ud.replay.ydata)&&~isempty(ud.replay.zdata)
                     ud.replay.hwhere = plot3(ud.replay.xdata(1), ud.replay.ydata(1), ud.replay.zdata(1), 'ro');
                  end
               end
               if ishandle(ud.replay.hwhere)
                   set(ud.replay.hwhere,'MarkerFaceColor',[1 0 0], 'tag','hwhere');
               end

            end
      end

      set(hax, 'NextPlot', oldnext);
   end

   set(hax, 'userdata', ud);
end

%set(get(hax,'parent'),'erasemode','normal');
function replayend(hax,closedlg)
if closedlg
    replaycallback(hax, '', 1);
end

if ishandle(hax)
   i_figure(get(hax, 'parent'));
   ud = get(hax, 'userdata');
   if ishandle(ud.replay.hwhere)
      delete(ud.replay.hwhere);
      ud.replay.hwhere = -1;
   end

   ud.replay.tdata = [];
   ud.replay.xdata = [];
   ud.replay.ydata = [];
   ud.replay.zdata = [];
   set(hax, 'userdata', ud);
end

function t = replaycallback(hax, avar,relt)
t = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, 't')
   ud = get(hax, 'userdata');
   t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
   if isfield(ud, 'replay')
      if ~isfield(ud.replay, 'xdata')||isempty(ud.replay.xdata)
         replaystart(hax);
         ud = get(hax, 'userdata');
      end

      tt = ud.replay.tdata;
      tnow = find(tt >= t, 1);
      if ~isempty(tnow)
          t=tt(tnow);
      end

      if ~isempty(tt)&&all(ishandle(ud.replay.HLine))&&(ishandle(ud.replay.hwhere)||ismember(ud.replay.type, [4, 5]))
         ndx=tt <= t;
         switch ud.replay.type
          case 4
            set(ud.replay.HLine, 'XData', ud.replay.xdata(ndx),'YData', ud.replay.ydata(ndx));
            set(ud.replay.hwhere,'Xdata',ud.replay.xdata(tnow),'Ydata',ud.replay.ydata(tnow));
         case 5
            set(ud.replay.HLine, 'XData', ud.replay.tdata(ndx),'YData', ud.replay.xdata(ndx));
          case 3
            set(ud.replay.HLine, 'XData', ud.replay.xdata(ndx),'YData', ud.replay.ydata(ndx),...
               'ZData', ud.replay.zdata(ndx));
            set(ud.replay.hwhere,'Xdata',ud.replay.xdata(tnow),'Ydata',ud.replay.ydata(tnow),...
               'Zdata', ud.replay.zdata(tnow));
          case 2
            set(ud.replay.HLine, 'XData', ud.replay.xdata(ndx),'YData', ud.replay.ydata(ndx));
            set(ud.replay.hwhere,'Xdata',ud.replay.xdata(tnow),'Ydata',ud.replay.ydata(tnow));
          case 1
            set(ud.replay.HLine, 'XData', ud.replay.xdata(ndx),'Ydata',zeros(sum(ndx),1));
            set(ud.replay.hwhere,'Xdata',ud.replay.xdata(tnow),'Ydata',zeros(size(tnow)));
         end

      end

   end

end

