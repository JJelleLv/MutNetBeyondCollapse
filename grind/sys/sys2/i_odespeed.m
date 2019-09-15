function status = i_odespeed(~, y, flag)
global g_grind;
status = 0;                             % Assume stop button wasn't pushed.
if nargin < 3 || isempty(flag)      % odephas2(t, y)
   f = get(0, 'CurrentFigure');
   if ~isempty(f)
      ud = get(f, 'UserData');
      if g_grind.truncate
         minx = min(y(ud.iX, :));
         maxx = max(y(ud.iX, :));
         miny = min(y(ud.iY, :));
         maxy = min(y(ud.iY, :)); 
         if ~isempty(minx) && ~isempty(miny) && ((minx < g_grind.xaxis.lim(1)) || (maxx > g_grind.xaxis.lim(2)) || ...
               (miny < g_grind.yaxis.lim(1)) || (maxy > g_grind.yaxis.lim(2)))
            ud.stop = 1;
         end

      end

      set(f, 'UserData', ud);
      if ~isempty(ud) && isfield(ud, 'stop') && ud.stop                       % Has stop button been pushed?
         status = 1;
         i_odespeed([], [], 'done');
      end
   end

else
   switch(flag)
    case 'init'                      % odephas2(tspan,y0,'init')
      
      f = get(0, 'CurrentFigure');
      if ~isempty(f)
         ud = get(f, 'userdata');
         ud.iX = i_varno(g_grind.xaxis.var);
         ud.iY = i_varno(g_grind.yaxis.var);
         ud.iZ = i_varno(g_grind.zaxis.var);
         % The STOP button.
         
         h = findobj(f,'Tag','stop');
         if isempty(h)
            ud.stop = 0;
            pos = get(0, 'DefaultUicontrolPosition');
            pos(1) = pos(1) - 15;
            pos(2) = pos(2) - 15;
            str = 'ud=get(gcf,''UserData''); ud.stop=1; set(gcf,''UserData'',ud);set(findobj(gcf,''Tag'',''stop''),''Visible'',''off'')';
            uicontrol( ...
               'Style','push', ...
               'String','Stop', ...
               'Position',pos, ...
               'Callback',str, ...
               'Tag','stop');
         else
            set(h,'Visible','on');            % make sure it's visible
            if ishold
               oud = get(gcf, 'UserData');
               if ~isempty(oud) && isfield(oud, 'stop')
                  ud.stop = oud.stop;             % don't change old ud.stop status
               else
                  ud.stop = 0;
               end

            else
               ud.stop = 0;
            end
            drawnow;
         end
         set(f, 'UserData', ud);
      end

    case 'done'                           % odephas2([],[],'done')
      f = get(0, 'CurrentFigure');
      if ~isempty(f)
         ud = get(f, 'UserData');
         if ~isempty(ud) && isfield(ud, 'stop')
            set(findobj(f,'Tag','stop'),'Visible','off');
         end

      end

      %
   end
end

