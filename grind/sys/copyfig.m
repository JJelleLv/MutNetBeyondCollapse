%COPYFIG   Make a copy of a figure in MATLAB. 
%   Note that when you copy a figure to a figure with a high number, 
%   it is not erased with the command <a href="matlab:help era">era</a>. Subsequently, the figure can
%   be combined with the command <a href="matlab:help combfig">combfig</a>.
%
%   Usage:
%   COPYFIG TONR - Copy the current figure to a new figure with number TONR.
%   COPYFIG FROMNR TONR - Copy figure FROMNR to TONR.
%   COPYFIG('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'from' [handle] - original figure handle.
%     'to' [number] - number of the new figure.
%
%   See also combfig, era
%
%   Reference page in Help browser:
%      <a href="matlab:commands('copyfig')">commands copyfig</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function copyfig(varargin)
%(from,to)
fieldnams={'from','h','original figure handle.',get(0,'CurrentFigure');...
   'to','n','number of the new figure',1000}';
args=i_parseargs(fieldnams,'if nargs==1,deffields=''to'';else,deffields=''from,to'';end;','',varargin);
if ~exist('i_checkstr','file')
    addpath([grindpath filesep 'sys2']);
end    
if ~isfield(args,'from')
   args.from=get(0,'CurrentFigure');
end
if ~isfield(args,'to')
   prompt={'Figure number to copy (gcf=current):','Copy to:'};
   answer=inputdlg(prompt,'Copy figure',1,{'',''});
   if isempty(answer{1})
       answer{1}='gcf';
   end
   args=i_parseargs(fieldnams,'if nargs==1,deffields=''to'';else,deffields=''from,to'';end;','',answer); 
end
if ishandle(args.to)
   delete(args.to);
end
h=i_figure(args.to);
if ishandle(args.from)
%   ax=findobj(args.from,'type','axes');
   if strcmp(get(args.from,'type'),'axes')
       hax=args.from;
   else
       hax=get(args.from,'children');
   end
   copyobj(hax, h);
   axnew=findobj(h,'type','axes');
   set(axnew,'position', [0.1300 0.1100 0.7750 0.8150], 'Units','normalized');      

else
   error('GRIND:copyfig:UnknownFig','Copyfig: source figure doesn''t exist');
end

