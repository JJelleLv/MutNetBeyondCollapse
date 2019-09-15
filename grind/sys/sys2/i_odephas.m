function status = i_odephas(t, y, flag)
%ODEPHAS2 2-D phase plane ODE output function.
%   When the string 'odephas2' is passed to an ODE solver as the 'OutputFcn'
%   property, i.e. options = odeset('OutputFcn','odephas2'), the solver
%   calls ODEPHAS2(T,Y) after every timestep.  The ODEPHAS2 routine plots
%   the first two components of the solution it is passed as it is computed,
%   adapting the axis limits of the plot dynamically.  To plot two
%   particular components, specify their indices in the 'OutputSel' property
%   passed to the ODE solver.
%
%   At the start of integration, a solver calls ODEPHAS2(TSPAN,Y0,'init') to
%   initialize the output function.  After each integration step to new time
%   point T with solution vector Y the solver calls STATUS = ODEPHAS2(T,Y).
%   If the solver's 'Refine' property is greater than one (see ODESET), then
%   T is a column vector containing all new output times and Y is an array
%   comprised of corresponding column vectors.  The STATUS return value is 1
%   if the STOP button has been pressed and 0 otherwise.  When the
%   integration is complete, the solver calls ODEPHAS2([],[],'done').
%
%   See also ODEPLOT, ODEPHAS3, ODEPRINT, ODE45, ODE15S, ODESET.

%   Mark W. Reichelt and Lawrence F. Shampine, 3-24-94
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.18 $  $Date: 1997/11/21 23:31:07 $
global g_grind;
status = 0;                             % Assume stop button wasn't pushed.
chunk = 128;                            % Memory is allocated in chunks.
nvar = g_grind.statevars.dim;
hfig = get(0, 'CurrentFigure');
oldversion=isnumeric(hfig);
if nargin < 3 || isempty(flag)           % odephas2(t, y)
   if ~isempty(hfig)
      ud = get(hfig, 'UserData');
      if isempty(ud) || (length(fieldnames(ud)) < 10)
         ud = [];
         i_odephas([], [], 'done');
      else
         % Append y to ud.y, allocating if necessary.
         nt = length(t);
         chunk = max(chunk, nt);
         rows = size(ud.y, 1);
         oldi = ud.i;
         newi = oldi + nt;
         if newi > rows
            ud.y = [ud.y; zeros(chunk, nvar)];
         end
         ud.y(oldi + 1:newi, :) = real(transpose(y));
         ud.ncalled = ud.ncalled + 1;
         ud.i = newi;
         if g_grind.truncate
            minx = min(ud.y(1:newi, ud.iX));
            maxx = max(ud.y(1:newi, ud.iX));
            miny = min(ud.y(1:newi, ud.iY));
            maxy = max(ud.y(1:newi, ud.iY));
            if (minx < g_grind.xaxis.lim(1)) || (maxx > g_grind.xaxis.lim(2)) || ...
                  (miny < g_grind.yaxis.lim(1)) || (maxy > g_grind.yaxis.lim(2))
               ud.stop = true;
            end

         end

         set(hfig, 'UserData', ud);
         if ud.stop == 1                       % Has stop button been pushed?
            status = 1;
            i_odephas([], [], 'done');
         else
            % Rather than redraw all of the data every timestep, we will simply move
            % the line segments for the new data, not erasing.  But if the data has
            % moved out of the axis range, we redraw everything.
            %    xlim = get(gca,'xlim');
            %    ylim = get(gca,'ylim');
            % Replot everything if out of axis range or if just initialized.
            if oldversion
            if (oldi < 2)
               %|| ...
               %          (min(y(1,:)) < xlim(1)) || (xlim(2) < max(y(1,:))) || ...
               %          (min(y(2,:)) < ylim(1)) || (ylim(2) < max(y(2,:)))
               if ~ud.p3
                  set(ud.lines, ...
                     'Xdata',ud.y(1:newi,ud.iX), ...
                     'Ydata', ud.y(1:newi, ud.iY));
                  set(ud.line, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata', ud.y(oldi:newi, ud.iY));
               else
                  set(ud.lines, ...
                     'Xdata',ud.y(1:newi,ud.iX), ...
                     'Ydata', ud.y(1:newi, ud.iY), ...
                     'Zdata', ud.y(1:newi, ud.iZ));
                  set(ud.line, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata', ud.y(oldi:newi, ud.iY), ...
                     'Zdata', ud.y(oldi:newi, ud.iZ));
               end

            elseif ishandle(ud.lines)
               % Plot only the new data.
               %   set(ud.line,'Color',[1,0,0]...     % "erase" old segment
               %       'Xdata',ud.y(oldi:newi,1), ...
               %       'Ydata',ud.y(oldi:newi,2), ...
               if ~ud.p3
                  set(ud.lines, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata',ud.y(oldi:newi,ud.iY), ...
                     'Color', g_grind.pen.drawcolor);
                  if newi > 25
                     set(ud.line, ...
                        'Xdata',ud.y(newi - 25:newi,ud.iX), ...
                        'Ydata',ud.y(newi - 25:newi,ud.iY), ...
                        'Color', [0, 0, 0]);
                  else
                     set(ud.line, ...
                        'Xdata',ud.y(1:newi,ud.iX), ...
                        'Ydata',ud.y(1:newi,ud.iY), ...
                        'Color', [0, 0, 0]);
                  end

               else
                  set(ud.lines, ...
                     'Xdata',ud.y(oldi:newi,ud.iX), ...
                     'Ydata',ud.y(oldi:newi,ud.iY), ...
                     'Zdata', ud.y(oldi:newi, ud.iZ), ...
                     'Color', g_grind.pen.drawcolor);
                  if newi > 25
                     set(ud.line, ...
                        'Xdata',ud.y(newi - 25:newi,ud.iX), ...
                        'Ydata',ud.y(newi - 25:newi,ud.iY), ...
                        'Zdata', ud.y(newi - 25:newi, ud.iZ), ...
                        'Color', [0, 0, 0]);
                  end

               end
            end

            else
                %MATLAB Version 2014b or more
               if ~ud.p3
                  ud.lines.addpoints(ud.y(oldi:newi,ud.iX),ud.y(oldi:newi, ud.iY));          
                  set(ud.line, 'xdata',ud.y(newi,ud.iX),'ydata',ud.y(newi,ud.iY));
               else
                   ud.lines.addpoints(ud.y(oldi:newi,ud.iX),ud.y(oldi:newi, ud.iY),ud.y(oldi:newi, ud.iZ));
                   set(ud.line, 'xdata',ud.y(newi,ud.iX),'ydata',ud.y(newi,ud.iY),'zdata',ud.y(newi,ud.iZ));
               end

             %  drawnow('nocallbacks');
            end
         end
      end
   end

else
   switch(flag)
    case 'init'                           % odephas2(tspan,y0,'init')
      ud = get(gcf, 'userdata');
      ud.y = zeros(chunk, nvar);
      ud.i = 1;
      ud.ncalled = 1;
      ud.stop = false;
      ud.y(1, :) = real(transpose(y));
      ud.iX = i_varno(g_grind.xaxis.var);
      ud.iY = i_varno(g_grind.yaxis.var);
      ud.iZ = i_varno(g_grind.zaxis.var);
      % Rather than redraw all data at every timestep, we will simply move
      % the last line segment along, not erasing it.
      if ~isempty(hfig)
         hax = get(hfig, 'CurrentAxes');
         
         isnew = length(get(hfig, 'Children')) < 2;
         if isnew
            set(hax, 'Xlim', g_grind.xaxis.lim);
            set(hax, 'Ylim', g_grind.yaxis.lim);
         end

         ud.p3=hfig == i_figno('phase3');
         oldhold=ishold;
         hold on;
         if oldversion
             if ud.p3
                 ud.lines = plot3(hax,y(ud.iX),y(ud.iY),y(ud.iZ),'-','Color',g_grind.pen.drawcolor,'EraseMode','none','tag','animate');
                 ud.line = plot3(hax,y(ud.iX),y(ud.iY),y(ud.iZ),'-','Color',[0,0,0],'EraseMode','xor','tag','animate');
                 set(hax, 'drawmode','fast');
                 if isnew
                     set(hax, 'Zlim', g_grind.zaxis.lim);
                     set(hax,'View',[322.5, 30]);
                 end

             else
                 ud.lines = plot(hax,y(ud.iX),y(ud.iY),'-','Color',g_grind.pen.drawcolor,'EraseMode','none','tag','animate');
                 ud.line = plot(hax,y(ud.iX),y(ud.iY),'-','Color',[0,0,0],'EraseMode','xor','tag','animate');
            end

         else
             %after MATLAB 2014b
             if ud.p3
                 %note MATLAB R2014 cannot have a handle in the first
                 %position, use 'parent' instead
                 ud.lines = animatedline(y(ud.iX),y(ud.iY),y(ud.iZ),'parent',hax,'LineStyle','-','Color',g_grind.pen.drawcolor,'tag','animate');
                 ud.line=plot3(hax,y(ud.iX),y(ud.iY),y(ud.iZ),'Color',g_grind.pen.drawcolor,'marker','.','markersize',10,'tag','animate');
                 i_grindlegend(-1,ud.lines)
                 i_grindlegend(-1,ud.line);
                 set(hax,'SortMethod','depth');
                 if isnew
                     set(hax, 'Zlim', g_grind.zaxis.lim);
                     set(hax,'View',[322.5, 30]);
                 end

             else
                 %note MATLAB R2014 cannot have a handle in the first position
                 ud.lines = animatedline(y(ud.iX),y(ud.iY),'parent',hax,'LineStyle','-','Color',g_grind.pen.drawcolor,'tag','animate');
                 ud.line=plot(hax,y(ud.iX),y(ud.iY),'Color',g_grind.pen.drawcolor,'marker','.','markersize',10,'tag','animate');
                 i_grindlegend(-1,ud.lines)
                 i_grindlegend(-1,ud.line);
               %  ud.line =animatedline(y(ud.iX),y(ud.iY),'LineStyle','-','Color',g_grind.pen.drawcolor, 'MaximumNumPoints',3);
          %       ud.line = animatedline(y(ud.iX),y(ud.iY),'LineStyle','-','Color',[0 0 0]);
             end

         end

             
         if ~oldhold
           hold off
         end

         if oldversion
            set(hax,'DrawMode','fast');
         else
            set(hax,'SortMethod','depth');
         end

         % The STOP button.
         h = findobj(hfig,'Tag','stop');
         if isempty(h)
            ud.stop = 0;
            pos = get(0, 'DefaultUicontrolPosition');
            pos(1) = pos(1) - 15;
            pos(2) = pos(2) - 15;
        %    str = 'ud=get(gcf,''UserData''); ud.stop=1; set(gcf,''UserData'',ud);set(findobj(gcf,''Tag'',''stop''),''Visible'',''off'')';
            uicontrol( ...
               'Style','push', ...
               'String','Stop', ...
               'Position',pos, ...
               'Callback',@stopcallback, ...
               'Visible','on',...
               'Tag','stop');
            set(hfig, 'DeleteFcn', @stopcallback);
         else
            set(h,'Visible','on');            % make sure it's visible
            if ishold
               oud = get(hfig, 'UserData');
               ud.stop = oud.stop;             % don't change old ud.stop status
            else
               ud.stop = 0;
            end
         end
         set(hfig, 'UserData', ud);
         drawnow;
      end

    case 'done'                           % odephas2([],[],'done')
      if ~isempty(hfig)
         ud = get(hfig, 'UserData');
         ud.stop=false;
         set(hfig,'userdata',ud);
         %  ud.y = ud.y(1:ud.i,:);
         %  set(f,'UserData',ud);
         %  set(ud.lines, ...
         %      'Xdata',ud.y(:,1), ...
         %      'Ydata',ud.y(:,2));
         if ~isempty(ud) && isfield(ud, 'line')
            lines=findobj(hfig,'tag','animate');
            delete(lines);
            if ~ishold
               set(findobj(hfig,'Tag','stop'),'Visible','off');
               refresh;                          % redraw figure to remove marker frags
            end
         end

      end
      set(findobj(hfig,'Tag','stop'),'Visible','off');
      %
   end
end
if exist('ud','var')&&~isempty(ud)&&isfield(ud,'ncalled')
   if (g_grind.slowdown < 0.01)
      if imod(ud.ncalled, 10) == 0
          drawnow;
     %    pause(g_grind.slowdown)
     % else
     %    drawnow;
      end

   elseif (g_grind.slowdown < 0.02)
      if imod(ud.ncalled, 5) < 2
         pause(g_grind.slowdown / 2)
      else
         drawnow;
      end

   else
      pause(g_grind.slowdown / 10);
   end


end
function stopcallback(hObject,~)
hfig=getparentfig(hObject);
ud=get(hfig,'UserData'); 
ud.stop=1; 
set(hfig,'UserData',ud);
set(findobj(hfig,'Tag','stop'),'Visible','off');


function m = imod(x, y)
%fast simple mod (for integers)
m = x - y .* floor(x ./ y);
