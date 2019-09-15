%SETCOLORMAP   set the colormap to a color gradient
%    Use this function to adapt the colorbar in 3D graphs and to set the default 
%    colormap in GRIND. (It can also be used outside GRIND) You set the bottom 
%    color (i.e. the color of the lowest value) and the top color and the range 
%    is made based on a gradient in the hue-saturation-value of both colors.
%  
%    Usage: 
%    SETCOLORMAP - sets the color map of the current figure. First you are asked to select
%    the bottom color, thereafter the top color.
%    SETCOLORMAP N - sets the size of the color map to N (default N=64).
%    SETCOLORMAP [R1 G1 B1] [R2 G2 B2] - sets the colormap to a gradient in 2 RGB colors.
%    Res=SETCOLORMAP(N,[R1 G1 B1],[R2 G2 B2]) - saves the colormap (size N) to Res
%    SETCOLORMAP('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'colorb' [color code] - color of the lowest value (bottom)
%     'colort' [color code] - color of the highest value (top)
%     'n' [integer>0] - number of classes in the colormap (64 by default)
%
%
%See also viewcells, <a href="matlab:help surf">surf</a>, <a href="matlab:help contourf">contourf</a>, <a href="matlab:help colormap">colormap</a>, <a href="matlab:help colorbar">colorbar</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setcolormap')">commands setcolormap</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function result = setcolormap(varargin)
%(m, bottomcolor, topcolor)
%hsv = 1;
fieldnams={'n', 'i>0', 'number of classes in the colormap (64 by default)',64;...
   'colorb', 'U1', 'color of the lowest value (bottom)', [0 1 1];...
   'colort', 'U1', 'color of the highest value (top)',[1 0 1]}';
args=i_parseargs(fieldnams,'if(nargs==2),deffields=''colorb,colort'';else,deffields=''n,colorb,colort'';end','',varargin,false,{@i_iscolor});
g = get(0, 'currentfigure');
global g_grind;
if isempty(g)
   if ~isempty(g_grind)
      cm = i_grindcolormap;
   else
      cm = zeros(64, 3);
   end
else
   cm = get(gcf, 'colormap');
end
if ~isfield(args,'n')
   args.n = size(cm, 1);
end
if (length(args.n)==1)&&(args.n<2)&&(nargin<3)
   args.n=[args.n,args.n,args.n]; %grays assumed;
end
   if ~isfield(args,'colorb')
      args.colorb = uisetcolor(cm(args.n, :), 'Select the bottom color');
   else
      if (length(args.colorb)==1)&&args.colorb<2
         args.colorb=[args.colorb,args.colorb,args.colorb];
      end
   end
   if ~isfield(args,'colort')
      args.colort = uisetcolor(cm(args.n, :), 'Select the top color');
   else
      if (length(args.colort)==1)&&args.colort<2
         args.colort=[args.colort,args.colort,args.colort];
      end
   end
col = [args.colorb; args.colort];
if ~isempty(g_grind) && isfield(g_grind, 'pen')
   g_grind.pen.colormap = col;
end
map = i_grindcolormap(args.n, col);
if (nargout < 1)&& ~isempty(g)
   colormap(map);
else
   if nargout<1
      disp('No model selected, cannot select colormap for grind');
   end
   result = map;
end
