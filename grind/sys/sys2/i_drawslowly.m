function i_drawslowly(hfig,x, y, hashead)
%COMET  Comet-like trajectory.
%   COMET(Y) displays an animated comet plot of the vector Y.
%   COMET(X,Y) displays an animated comet plot of vector Y vs. X.
%   COMET(X,Y,p) uses a comet of length p*length(Y).  Default is p = 0.10.
%
%   Example:
%       t = -pi:pi/200:pi;
%       comet(t,tan(sin(t))-sin(tan(t)))
%
%   See also COMET3.

%   Charles R. Denham, MathWorks, 1989.
%   Revised 2-9-92, LS and DTP; 8-18-92, 11-30-92 CBM.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.8 $  $Date: 1997/11/21 23:46:03 $
global g_grind;
if nargin == 0
   error('GRIND:drawslowly:ArgError','Not enough input arguments.');
end
if nargin < 3
   y = x;
   x = 1:length(y);
end
if nargin < 4
   hashead = 0;
end

if any(~isreal(x))||any(~isreal(y))
    x=real(x);
    y=real(y);
    warning('grind:ignorecomplex','Imaginary parts of complex values ignored');
end
hax=get(hfig,'CurrentAxes');
%ax = gca;
%if ~ishold,
%   axis([min(x(isfinite(x))) max(x(isfinite(x))) ...
%      min(y(isfinite(y))) max(y(isfinite(y)))])
%end

if ~isoctave&&verLessThan('matlab','8.4.0')
   %if ~g_grind.version.isoctave
   ud = get(hfig, 'userdata');
   h = findobj(hfig,'Tag','stop');
   if isempty(h)
      ud.stop = 0;
      pos = get(0, 'DefaultUicontrolPosition');
      pos(1) = pos(1) - 15;
      pos(2) = pos(2) - 15;
      str = 'ud=get(gcf,''UserData''); ud.stop=1; set(gcf,''UserData'',ud);';
      uicontrol( ...
         'Style','push', ...
         'String','Stop', ...
         'Position',pos, ...
         'Callback',str, ...
         'Tag','stop');
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
   %end

   xlim1  = get(hax, 'Xlim');
   xlim2  = [max(x(:)), min(x(:))];
   xlim([min([xlim1(1), xlim2(1)]), max([xlim1(2), xlim2(2)])]);
   if hashead
      head = line('parent',hax,'color',[0 0 0],'marker','o','erase','xor', ...
         'xdata',x(1),'ydata',y(1));
   end
   body = line('parent',hax,'color',[0 1 1],'linestyle','-','erase','none', ...
      'xdata',[],'ydata',[]);
   tail = line('parent',hax,'color',g_grind.pen.drawcolor,'linestyle','-','erase','none', ...
      'xdata',[],'ydata',[]);
   m = length(x);
   k = 5;
   if k > m - 1
      k = m - 1;
   end

   % Grow the body
   for i = 2:k + 1
      j = i - 1:i;
      if hashead
         set(head,'xdata',x(i),'ydata',y(i));
      end

      set(body,'xdata',x(j),'ydata',y(j))
      drawnow
   end
   set(hfig, 'pointer', 'watch');
   
   % Primary loop
   for i = k + 2:m
      j = i - 1:i;
      if hashead
         set(head,'xdata',x(i),'ydata',y(i))
      end

      set(body,'xdata',x(j),'ydata',y(j))
      set(tail,'xdata',x(j-k),'ydata',y(j-k))
      slowdown(i)
      %   if ~g_grind.version.isoctave
      ud = get(hfig, 'UserData');
      if ud.stop
         break;
      end

      %   end

   end
   
   % Clean up the tail
   for i = m + 1:m + k
      j = i - 1:i;
      set(tail,'xdata',x(j-k),'ydata',y(j-k))
      drawnow
   end
   if hashead
      delete(head);
   end

   delete(body);
   delete(tail);
   %if ~g_grind.version.isoctave
   set(findobj(hfig,'Tag','stop'),'Visible','off');
   set(hfig, 'pointer', 'arrow');
   %end

else
   %MATLAB >2014
   %note MATLAB R2014 cannot have a handle in the first position
   h = animatedline(x(1), y(1),'parent',hax, 'Color', g_grind.pen.drawcolor);
   if hashead
      h2=plot(hax,x(1),y(1),'ro','MarkerFaceColor','r');
   end

   for i = 2:length(x)
      h.addpoints(x(i), y(i));
      if hashead
         set(h2,'xdata',x(i),'ydata',y(i));
      end

      if mod(i, 10) == 0
         drawnow;
      end

      slowdown(i)
   end

   if hashead
      delete(h2);
   end

   delete(h);
end

function slowdown(ncalled)
global g_grind;
if (g_grind.slowdown < 0.01)
   if mod(ncalled, 10) == 0
      pause(g_grind.slowdown)
%    else
%       drawnow;
   end

elseif (g_grind.slowdown < 0.02)
   if mod(ncalled, 1) < 2
      pause(g_grind.slowdown / 2)
   else
      drawnow;
   end

else
   pause(g_grind.slowdown / 10);
end


