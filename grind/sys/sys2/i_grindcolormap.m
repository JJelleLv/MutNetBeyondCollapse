function c=i_grindcolormap(m, col)
%Default colormap of grind
%
global g_grind;
if nargin < 1
   g = get(0, 'currentfigure');
   if isempty(g)
      m=64; 
   else 
      m = size(get(g,'colormap'),1);
   end
end

if nargin < 2
   col=g_grind.pen.colormap; 
end

s=size(col,1);
col=rgb2hsv(col);
%c=ones(m,3);
x=1:s;
x2=1:(s-1)/(m-1):s;
c=interp1(x,col,x2);
c=hsv2rgb(c);

