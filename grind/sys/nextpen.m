%NEXTPEN   Set next pen for phase plane
%   Assign the next color for g_grind.pen (the color of the next trajectory in the
%   phase plane).
%
%   Usage:
%   NEXTPEN - assign the next color
%   NEXTPEN APEN - select APEN as pen
%   NEXTPEN('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'color' [color code or empty] - color code
%     'cycle' [logical] - change the colors
%     'linestyle' [- | : | -. | -- or empty] - line style
%     'linewidth' [number>0] - Line width (pixels)
%     'markerfacecolor' [color code or none | auto] - marker face color code
%     'markersize' [number>0] - marker size (pixels)
%     'markerstyle' [. | o | x | + | * | s | d | v | ^ | < | > | p | h or empty] - marker style
%     'pen' [pen] - short code for line and marker style (for instance "g:")
%
%   See also setpen
%
%   Reference page in Help browser:
%      <a href="matlab:commands('nextpen')">commands nextpen</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function pen = nextpen(varargin)
global g_grind;
if nargin == 0
   if isempty(g_grind)
      args.apen=[];
   else
      args.apen = g_grind.pen;
   end
else
   if nargin==1&&isstruct(varargin{1});
       pen=varargin{1};
       fields=  {'color','cycle', 'linestyle', 'linewidth','markersize','markerstyle','markerfacecolor','pen'};
       f=fieldnames(pen);
       ndx=ismember(f,fields);
       varargin{1}=rmfield(pen,f(~ndx));
   end
   args.apen=setpen(varargin{:});
%    fieldnams={'apen', '', 'select a pen (string) or use a structure to select the next pen',''}';
%    args=i_parseargs(fieldnams,'apen','',varargin);
end
pen1=i_nextpen(args.apen);
if nargout>0
    pen=pen1;
end