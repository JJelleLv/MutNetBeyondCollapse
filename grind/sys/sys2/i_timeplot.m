% i_timeplot initial function to create a time plot
function [] = i_timeplot(hol)
global g_grind g_t;
for No = 1:size(g_grind.timevars, 2)
   if ~isempty(g_grind.timevars{No})
      hfig = i_makefig('time', No - 1);
      if hol
          hold('on');
      end

      %plotedit('off');
      oldhold = ishold;
      apen = g_grind.pen;
      apen.i = 1; %  do not cycle colors between runs
      apen.cycle = true;
      if ~oldhold
         h=get(get(0,'currentfigure'),'currentaxes');
         if ~isempty(h)
            delete(h);
         end

      else
         h=get(get(0,'currentfigure'),'currentaxes');
         series = get(h, 'children');
         for i = 1:length(series)
            apen = i_nextpen(apen);
         end

      end

      set(gca,'YLimMode','auto');
      hold('on');
      tvar = g_grind.outt{No};
      tt = outfun(struct('fun',tvar,'opts',{{}}));
      for i = 1:size(g_grind.timevars{No}, 2)
         yvar = g_grind.timevars{No}{i};
         pen = apen.pen;
         if strncmpi(yvar,'observed ',9)|| strcontains(yvar,'observed(')&&strcmp(tvar, 't')
            tt1 = outfun('observed(''t'')');
            pen = 'o';
         else
            tt1 = tt;
         end

       %  res = outfun(yvar);
         res=outfun(struct('fun',yvar,'opts',{{}}));
         if sum(isnan(res))>length(res)/10
             pen='.-';
         end

         if isempty(res)
            warning('GRIND:time:nodata','No data for "%s"\n', g_grind.timevars{No}{i})
         elseif size(res, 2) == 1
             if size(res,1)~=size(tt1,1)&& strcontains(yvar,'observed')
                hplot=plot(outfun('observed(''t'')'),res, pen, 'Color', apen.color,'Linewidth',apen.linewidth);
             else
                hplot=plot(tt1,res, pen, 'Color', apen.color,'Linewidth',apen.linewidth);
             end

            i_grindlegend(1, hplot, g_grind.timevars{No}{i});
         else
            for j = 1:size(res, 2)
               hplot(j)=plot(tt1,res(:,j), apen.pen, 'Color', apen.color,'Linewidth',apen.linewidth);
               apen = i_nextpen(apen);
            end

         %   i_grindlegend(1, hplot, g_grind.timevars{No}{i});
         %   end

         end

         apen = i_nextpen(apen);
      end

      if ~oldhold
         hold('off');
      end

      tnam = g_grind.diffto;
      f = strfind(tnam, '(');
      tim = tnam(1:f(1) - 1);
      set(hfig,'Name',[tim 'plot']);
      hax = get(hfig, 'CurrentAxes');
      i_plotdefaults(hfig);
      lim = get(hax, 'Ylim');
      if ~isempty(lim)
         if (lim(2) ~= 0) && ((lim(2) - lim(1)) / lim(2) < 0.001)
            n = abs(0.01 * lim(2));
            if (lim(1) > n)
               lim(1) = lim(1) - n;
            elseif lim(1) > 0
               lim(1) = 0;
            end

            lim(2) = lim(2) + n;
            if lim(1) == lim(2)
               lim(2) = lim(1) + eps;
            end

            set(hax, 'Ylim', lim);
         end
      end

      if ~strcmp(g_grind.outt{No}, 't')
         %         xlabel([tim '(' g_grind.outt{No} ')']);
         xlabel( g_grind.outt{No});
      else
         xlabel(tnam);
      end

      i_grindlegend(11,hax);
      ud = get(hax, 'userdata');
      ud.replay.callback = @replaycallback;
      ud.replay.onstart = @replaystart;
      ud.replay.onend = @replayend;
      if isempty(g_t)
         ud.replay.settings=struct('tvar','t','tlim',[0 g_grind.ndays],'numt',0);
      else
         ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
      end

      if g_grind.solver.isdiffer&&g_grind.statevars.dim==1
        ud.replay.settings.numt=ud.replay.settings.numt*2-1;
      end

      ud.replay.opt = No;
      set(hax,'userdata', ud);
   %   plotedit('on');
      %   end

   end

end


function replaystart(hax)
global g_grind g_t t;
if ishandle(hax)
   N0 = i_initvar;
   if i_settingschanged(N0, g_grind.ndays)
      i_ru(t, g_grind.ndays, N0, 1);
   end

   i_figure(get(hax,'parent'));
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
       if g_grind.solver.isdiffer&&g_grind.statevars.dim==1 %for cobwebbing!
        ud.replay.settings.numt=ud.replay.settings.numt*2-1;
      end

      No = ud.replay.opt;
      ser = get(hax, 'children');
      ud.replay.ydata = cell(size(ser));
      ud.replay.xdata = cell(size(ser));
      if No>length(g_grind.outt)
          No=length(g_grind.outt);
      end

      tvar = g_grind.outt{No};
      ud.replay.tdata=outfun('t');
      k = size(g_grind.timevars{No}, 2);
      for i = 1:size(g_grind.timevars{No}, 2)
         yvar = g_grind.timevars{No}{i};
         if (strncmpi(yvar,'observed ',9)||strncmpi(yvar,'observed(',9))&&strcmp(tvar, 't')
            tt1 = outfun('observed t');
            ist=false;
         else
            tt1 = outfun(tvar);
            ist=strcmp(tvar,'t');
         end
         if any(~isreal(tt1))
             warning('grind:time','Imaginary parts of complex values ignored');
             tt1=real(tt1);
         end
         res = outfun(yvar);
         if any(~isreal(res))
             warning('grind:time','Imaginary parts of complex values ignored');
             res=real(res);
         end
         if (size(res, 2) > 1) && (size(tt1 ,2)==1)
             tt1 = repmat(tt1, 1, size(res, 2));
         end

         if ist
             ud.replay.xdata{k}=[];
         else
            ud.replay.xdata{k} = tt1;
         end

         ud.replay.ydata{k} = res;
         k = k - 1;
      end

   end

   set(hax, 'userdata', ud);
   if ~isoctave&&verLessThan('matlab','8.4.0')
       set(hax, 'drawmode','fast');
   else
       set(hax,'SortMethod','depth');
   end

   if length(g_t)>1&&strcmp(tvar,'t')
      set(hax,'xlim',[min(ud.replay.tdata) max(ud.replay.tdata)]);
      set(hax,'xlimmode','manual');
   end

end


function replayend(hax,closedlg)
if closedlg
    replaycallback(hax,'t',1);
end

if ishandle(hax)   
   i_figure(get(hax,'parent'));
   set(hax,'xlimmode','auto');
   ud = get(hax, 'userdata');
   ud.replay.xdata = {};
   ud.replay.ydata = {};
   set(hax, 'userdata', ud);
end

function t=replaycallback(hax,avar, relt)
t=[];
if ishandle(hax)&&isempty(avar)||strcmp(avar,'t')
   ud = get(hax, 'userdata');
   t=ud.replay.settings.tlim(1)+relt*(ud.replay.settings.tlim(end)-ud.replay.settings.tlim(1));
   if isfield(ud, 'replay')
      if ~isfield(ud.replay, 'xdata')
         replaystart(hax);
         ud = get(hax, 'userdata');
      end

      ser = get(hax, 'children');
      for i = 1:min(length(ser),length(ud.replay.xdata))
         if isempty(ud.replay.xdata{i})
             tt = ud.replay.tdata;
         else
             tt=ud.replay.xdata{i};
         end

         if length(tt)==length(ud.replay.tdata)
             ndx=ud.replay.tdata<=t;
         else
             ndx=tt<=t; %if t is observed data this is necessary
         end

         if any(ndx)&&length(ndx)==length(ud.replay.ydata{i})
            tt1=tt(ndx);
            set(ser(i), 'XData', tt1, 'YData', ud.replay.ydata{i}(ndx));
         end

      end

   end

end




