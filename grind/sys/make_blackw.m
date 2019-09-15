function make_blackw(minim,blackw)
if nargin<2
    blackw=1;
end
if nargin < 1
   minim = 0.005;
else
   minim = i_checkstr(minim);
end
box on;
lines = get(gca, 'children');
YLim = get(gca, 'ylim');
XLim = get(gca, 'xlim');
set(gca, 'linewidth', 1.5);
%l = -1;
for i = length(lines):-1:1
   if ~isempty(lines(i)) && (strcmp(get(lines(i),'Type'),'line'))
      ys = get(lines(i), 'Ydata');
      xs = get(lines(i), 'Xdata');
      toosmall = 1;
      for j = 1:length(ys)
         if (ys(j) / (YLim(2) - YLim(1)) > minim) && (xs(j) / (XLim(2) - XLim(1)) > minim)
            toosmall = 0;
         end
      end
      if toosmall
         delete(lines(i))
         disp(['Deleted line ' num2str(i)]);
      else
         l = length(lines)-i+1;
         switch mod(l, 4)
          case 0
            linestyle = '-';
          case 1
            linestyle = ':';
          case 2
            linestyle = '-.';
          otherwise
            linestyle = '--';
         end
         set(lines(i), 'linestyle', linestyle);
         if blackw
             if l < 16
                set(lines(i), 'Color', [0, 0, 0]);
             else
                set(lines(i), 'Color', [0.7, 0.7, 0.7]);
             end
         end
         if l < 4
            linewidth = 1;
         elseif l < 8
            linewidth = 1.5;
         elseif l < 12
            linewidth = 2;
         elseif l < 16
            linewidth = 1;
         elseif l < 20
            linewidth = 1.5;
         else
            linewidth = 2;
         end
         disp([num2str(l) ' ' num2str(max(ys)) linestyle ' ' num2str(linewidth)]);
         set(lines(i), 'linewidth', linewidth)
      end
   end
end



