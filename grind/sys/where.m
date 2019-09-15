%WHERE   Show initial conditions
%   display the current initial condition in the current plot.
%
%   Usage:
%   WHERE - shows the initial conditions as triangle (^)
%   WHERE START - shows the initial conditions as triangle (^)
%   WHERE MARK - shows the initial conditions with a marker MARK
%   . = point; o = circle; x = x-mark; + = plus; * = star; s = square;
%   d = diamond; v = triangle (down); ^ = triangle (up); < = triangle (left);
%   > = triangle (right); p = pentagram; h = hexagram
%   WHERE MARK FACECOLOR - sets a face color (MATLAB color code).
%   WHERE END - shows the final conditions as triangle (v)
%   WHERE END MARK FACECOLOR - shows the final conditions with a marker MARK
%   WHERE('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'blink' [logical] - blink once before showing the point
%     'marker' [pen] - kind of marker
%     'markerfacecolor' [color code or none | auto] - the color to fill of the marker
%     'time' [start | end | sblink] - time of the run to show (start/end)
%  
%  
%   See also null, ru 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('where')">commands where</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function h1=where(varargin)
%(apen, apen1)
global g_Y g_t g_grind t;
fieldnams={'time', 'e[start|end|sblink]', 'time of the run to show (start/end)','end';...
   'marker', 'U1', 'kind of marker','.';...
   'blink', 'l', 'blink once before showing the point',false;...
   'markerfacecolor','U2#e[none|auto]','the color to fill of the marker','none'}';
args=i_parseargs(fieldnams,'if argtype(1,''e[start|end|sblink]''),deffields=''time,marker,markerfacecolor'';else,deffields=''marker,markerfacecolor'';end;','',varargin,false,{@validpen,@i_iscolor});
if isfield(args,'time')&&~any(strcmpi(args.time,{'start','sblink','end'}))
    if isfield(args,'mark')
        h1=args.marker;
        args.marker=args.time;
        args.time=h1;
    else
       args.marker=args.time;
       args.time='end';
    end
end
if ~isfield(args,'time')
    args.time='end';
end
if ~isfield(args,'markerfacecolor')
    args.markerfacecolor='none';
end
if ~isfield(args,'blink')
    args.blink=false;
end
wherefinish = 0;
h=[];
if (nargin == 0) || (strcmpi(args.time, 'start'))
   if ~isfield(args,'marker')
      args.marker= '^';
   end
elseif strcmpi(args.time, 'sblink')
   args.blink=true;
   if ~isfield(args,'marker')
      args.marker= '^';
   end
elseif strcmpi(args.time, 'end')
   wherefinish = 1;
   if ~isfield(args,'marker')
      args.marker= 'v';
   end
end
i_parcheck;
oldY = g_Y;
oldt = g_t;
if wherefinish
   N0 = i_initvar;
   if i_settingschanged(N0)
      i_ru(t, g_grind.ndays, N0, 1);
   end
   oldY=g_Y;
   g_Y=g_Y(end,:); 
   g_t=g_t(end); 
else
   g_Y = transpose(i_initvar); 
   g_t=t;
end
oldpen = g_grind.pen;
try
   g_grind.pen.pen=args.marker;
   g_grind.pen.markersize=5;
   g_grind.pen.markerfacecolor = args.markerfacecolor;
%    if (length(args.marker)==2)   
%        g_grind.pen.pen=args.marker(2);
%        if (args.marker(1)=='k')
%          fcolor=[0 0 0];
%        elseif args.marker(1)=='g'
%          fcolor=[0.8 0.8 0.8];
%        else
%          fcolor=[1 1 1];
%        end
%    else
%       fcolor='none';
%       g_grind.pen.pen = args.marker;
%    end
   if args.blink
      g_grind.pen.color = [1 0 0];
      currfig=get(0,'currentfigure');
      if ~isempty(currfig)
          i_phas(currfig, 0);
          ch=get(gca,'children');
          if ~isempty(ch)
              h=ch(1);
              set(h,'tag','where')
              set(h,'userdata',[]);
              set(h,'markerfacecolor',[1 0 0]);
              pause(0.5);
              set(h,'markerfacecolor',args.markerfacecolor,'color',[0 0 0]);
          end
      end
  else
      g_grind.pen.color = [0 0 0];
      i_phas(gcf, 0);
      ch=get(gca,'children');
      h=ch(1);
   %   set(h,'markerfacecolor',fcolor); 
  end
   g_grind.pen = oldpen;
   g_Y = oldY;
   g_t = oldt;
catch err
   g_grind.pen = oldpen;
   g_Y = oldY;
   rethrow(err);
end
if nargout==1
    h1=h;
end

function [v, errormsg] = validpen(x)
if nargin==0 %helptext
    v='pen';
    return;
end
errormsg='';
v=x;
if ~ischar(x)
    errormsg = 'argument should be char';
else
    res=parsepen(x);
    lenpen=length(res.linestyle)+length(res.markerstyle)+length(res.color);
    if lenpen<length(x)
        errormsg=sprintf('Unknown or double characters in pen "%s"',x);
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