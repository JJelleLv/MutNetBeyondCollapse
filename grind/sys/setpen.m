%SETPEN   Set the pen for the plots
%   The user is prompted for the following information:
%   - Color of the pen - select from list
%   - Cycle pen in plots, if 'Yes' is selected each trajectory gets another
%     color
%   - Line style - select from list
%   - Style of marker of data points - select from list
%
%   Usage:
%   SETPEN - User is prompted for information.
%   SETPEN PEN - short code for line, marker style and color (for instance "k--o").
%   SETPEN PEN CYCLE - pen code and CYCLE (can be 1 (=true) or 0 (false)).
%   SETPEN('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'color' [color code or empty] - color code
%     'cycle' [logical] - change the colors
%     'linestyle' [- | : | -. | -- or empty] - line style
%     'linewidth' [number>0] - Line width (pixels)
%     'markersize' [number>0] - marker size (pixels)
%     'markerstyle' [. | o | x | + | * | s | d | v | ^ | < | > | p | h or empty] - marker style
%     'markerfacecolor' [color code or none | auto] - marker face color code
%     'pen' [pen] - short code for line and marker style (for instance "g:")
%
%
%      color               marker style             line style
%      y     yellow        .     point              -     solid
%      m     magenta       o     circle             :     dotted
%      c     cyan          x     x-mark             -.    dashdot 
%      r     red           +     plus               --    dashed   
%      g     green         *     star
%      b     blue          s     square
%      w     white         d     diamond
%      k     black         v     triangle (down)
%                          ^     triangle (up)
%                          <     triangle (left)
%                          >     triangle (right)
%                          p     pentagram
%
%                          h     hexagram
%
% 
%   See also RU, PHAS, NEXTPEN 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setpen')">commands setpen</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function thepen=setpen(varargin)
%(acolor, linestyle, markerstyle,cycle )
global g_grind;
if isfield(g_grind,'pen')
    apen=g_grind.pen;
    res=parsepen(apen.pen);
    apen.markerstyle=res.markerstyle;
    apen.color=g_grind.pen.color;
    apen.linestyle=res.linestyle;
else
    apen=struct('color','b','linestyle','-','markerstyle','','markerfacecolor','none','linewidth',1,'markersize',5,'cycle',true);
end
fieldnams={'color', 'U2#E', 'color code',apen.color;...
    'linestyle', 'e[-|:|-.|--]#E', 'line style',apen.linestyle;...
    'markerstyle', 'e[.|o|x|++|*|s|d|v|^|<|>|p|h]#E', 'marker style',apen.markerstyle;...
    'linewidth', 'n>0', 'Line width (pixels)',apen.linewidth;...
    'pen','U1','short code for line and marker style (for instance "g:")',invparsepen(apen);...
    'markerfacecolor', 'U2#e[none|auto]', 'marker face color code',apen.markerfacecolor;...
    'markersize', 'n>0', 'marker size (pixels)',apen.markersize;...
    'cycle', 'l', 'change the colors',apen.cycle}';
args=i_parseargs(fieldnams,'pen,cycle','',varargin,false,{@validpen,@i_iscolor});
i_parcheck;
oldpen=parsepen(g_grind.pen.pen);
colorcodes=struct('code',{'b',   'g',   'r',     'c',   'm',    'y',    'k'},...
    'color',{[0 0 1],[0 1 0],[1 0 0],[0 1 1],[1 0 1],[1 1 0],[0 0 0]});
if isfield(args,'pen')
    res=parsepen(args.pen);
    args.linestyle=res.linestyle;
    args.markerstyle=res.markerstyle;
    if ~isempty(res.color)
        args.color=res.color;
    end
end
if isfield(args,'linewidth')
    g_grind.pen.linewidth=args.linewidth;
end
if isfield(args,'markersize')
    g_grind.pen.markersize=args.markersize;
end
if ~isfield(args,'linestyle')
    args.linestyle=oldpen.linestyle;
end
if ~isfield(args,'markerstyle')
    args.markerstyle=oldpen.markerstyle;
end
if isfield(args,'markerfacecolor')
    g_grind.pen.markerfacecolor=args.markerfacecolor;
end
if ~isfield(args,'color')
    args.color=g_grind.pen.color;
end
if ~isfield(args,'cycle')
    args.cycle=g_grind.pen.cycle;
end
colorlist={'b   blue', 'g   green', 'r   red', 'c   cyan', 'm   magenta', 'y   yellow', 'k   black', '    other'};
linestylelist = {'    no line','-     solid',':     dotted','-.    dashdot','--    dashed'};
markstylelist = {'      no symbol', '.     point', 'o     circle', 'x     x-mark',...
    '+     plus', '*     star', 's     square', 'd     diamond', 'v     triangle (down)',...
    '^     triangle (up)', '<     triangle (left)', '>     triangle (right)', 'p     pentagram',...
    'h     hexagram'};
if nargin  == 0
    if isnumeric(args.color)&&length(args.color)==3
        for i=1:length(colorcodes)
            if all(colorcodes(i).color==args.color)
                args.color=colorcodes(i).code;
                break
            end
        end
    end
    if isempty(args.color)
        icolor=1;
    elseif isnumeric(args.color)
        icolor=8;
    else
        icolor=find(strncmp(args.color,colorlist,1));
    end
    icolor = listdlg('ListString', colorlist, ...
        'SelectionMode', 'single', ...
        'PromptString', 'Select a color', ...
        'InitialValue', icolor);
else
    if ischar(args.color)
        icolor=find(strncmp(args.color,colorlist,1));
    elseif isnumeric(args.color)
        icolor=[];
        g_grind.pen.color=args.color;
    end
end
if ~isempty(icolor)
    if icolor ~= 8
        oldcycle = g_grind.pen.cycle;
        g_grind.pen.color=colorcodes(icolor).color;
        g_grind.pen.i = icolor - 1;
        g_grind.pen.cycle = 1;
        nextpen;
        g_grind.pen.cycle = oldcycle;
    else
        g_grind.pen.i = 8;
        c = uisetcolor;
        if length(c) > 1
            g_grind.pen.color = c;
        end
    end
end
if nargin  == 0
    space=' ';
    linestyle=[args.linestyle, space(ones(1,2-length(args.linestyle)))];
    iline=find(strncmp(linestylelist,linestyle,2));
    if isempty(iline)
        iline=2;
    end
    iline = listdlg('ListString', linestylelist, ...
        'SelectionMode','single',...
        'PromptString','Select a line style',...
        'InitialValue', iline);
    if iline > 0
        args.linestyle = strtrim(linestylelist{iline}(1:2));
    end
elseif ~isfield(args,'linestyle')
    args.linestyle = '-';
end
if nargin == 0
    imark=find(strncmp(markstylelist,args.markerstyle,1));
    if isempty(imark)
        imark=1;
    end
    i = listdlg('ListString', markstylelist, ...
        'SelectionMode', 'single', ...
        'PromptString', 'Select a data point symbol', ...
        'InitialValue', imark);
    if i > 0
        args.markerstyle = markstylelist{i}(1);
    end
end

g_grind.pen.pen =  invparsepen(args);
g_grind.pen.pen2 =g_grind.pen.pen; %this is used to cycle line styles in time plots
g_grind.pen.color2=g_grind.pen.color;
if nargin  == 0
    c=questdlg('Cycle colors in plots?','GRIND');
    if ~strcmp(c, 'Cancel')
        g_grind.pen.cycle = strcmp(c, 'Yes');
    end
elseif isfield(args,'cycle')
    g_grind.pen.cycle = args.cycle;
end
if nargout>0
  thepen=g_grind.pen;
end
function res=invparsepen(args)
if isfield(args,'markerstyle')
    if isempty(args.linestyle)
        res =  strtrim(args.markerstyle);
    else
        res = [strtrim(args.linestyle) strtrim(args.markerstyle)];
    end
else
        res =strtrim(args.linestyle);
end

function [v, errormsg] = validpen(x)
if nargin==0 %helptext
    v='pen';
    return;
end
v=x;
errormsg='';
if ~ischar(x)
    errormsg = 'argument should be char';
else
    res=parsepen(x);
    lenpen=length(res.linestyle)+length(res.markerstyle)+length(res.color);
    if lenpen<length(x)
        errormsg=sprintf('Unknown characters in pen "%s"',x);
    end
end


function res=parsepen(apen)
acolor=regexp(apen,'[bgrcmykw]','match','once');
alinestyle=regexp(apen,'[-][-]','match','once');
if isempty(alinestyle)
    res=struct('color',acolor,'markerstyle',regexp(apen,'((?<![-])[\.])|([ox+*sdv^<>ph])','match','once'),...
        'linestyle',regexp(apen,'([-][-])|([-][.])|[-:]','match','once'));
else
    res=struct('color',acolor,'markerstyle',regexp(apen,'[.ox+*sdv^<>ph]','match','once'),'linestyle',alinestyle);
end



