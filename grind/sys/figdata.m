%FIGDATA   Extract xyz data from a figure
%  Use this function to export data from the currently selected 
%  figure to a matrix. 
%
%
%  Usage:
%  M=FIGDATA - copies the first series to matrix M. If there is a series 
%  selected, that series is copied.
%  M=FIGDATA(SERNO) - copies the SERNOth series to M.
%  FIGDATA SERNO - shows the SERNOth series.
%  FIGDATA('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'hax' [handle] - function handle to get the data (default gca)
%     'serno' [integer>0] - number of the series to copy
%   FIGDATA('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - gives an overview of the series of the current figure
%     '-a' - -all, gets one matrix of all series.
%     '-m' [SERNO SERNO]- merges two selected series.
%
%    
%  See also copyfig, varcopy
%
%   Reference page in Help browser:
%      <a href="matlab:commands('figdata')">commands figdata</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function A = figdata(varargin)
%(serno,ax)
fieldnams={'serno', 'i>0', 'number of the series to copy',1;...
   'hax', 'h', 'function handle to get the data (default gca)',[]}';
args=i_parseargs(fieldnams,'if(hasoption(''-a'')),deffields=''hax'';else,deffields=''serno,hax'';end;','-?,-m,-a',varargin);
if ~isfield(args,'hax')
    hfig=get(0,'CurrentFigure');   
    if isempty(hfig)
      error('GRIND:figdata:NoFig','No figure available');
   else
      args.hax=get(hfig,'CurrentAxes');
   end
end

series = get(args.hax, 'children');
ndx=~strcmp(get(series,'type'),'text');
series=series(ndx);
if ~isfield(args,'serno')
   sel=strcmp(get(series,'selected'),'on');
   if any(sel)
      args.serno=find(sel);
      fprintf('Using selected series, series %d\n',args.serno(1));
   else
      disp('Using series 1');
      args.serno=1;
   end
end
if any(strcmp(args.opts,'-m'))
   if length(args.serno)~=2
      error('GRIND:figdata:SelSeries','select two series to merge');
   end
   X1=getv(series(args.serno(1)),'xdata');
   Y1=getv(series(args.serno(1)),'ydata');
   Z1=getv(series(args.serno(1)),'zdata');
   X2=getv(series(args.serno(2)),'xdata');
   Y2=getv(series(args.serno(2)),'ydata');
   Z2=getv(series(args.serno(2)),'zdata');
   if X1(1)==X2(1)
      X1=flipud(X1);
      Y1=flipud(Y1);
      Z1=flipud(Z1);
   end
   set(series(args.serno(1)),'xdata',[X1;X2]); 
   set(series(args.serno(1)),'ydata',[Y1;Y2]); 
   set(series(args.serno(1)),'zdata',[Z1;Z2]);
   fprintf('series %d and %d merged\n',args.serno(1),args.serno(2));
   return
end
if any(strcmp(args.opts,'-a'))
   args.serno=1:length(series);
end
if any(strcmp(args.opts,'-?'))% ~isempty(flag)&&(flag == '?')
   fprintf('The current figure has %d series\n', length(series));
   for i=1:length(series)
      X=getv(series(i),'xdata');
      Y=getv(series(i),'ydata');
      Z=getv(series(i),'zdata');
      s=sprintf('Length series%3d: X: %6d, Y: %6d, Z: %6d',i,length(X),length(Y),length(Z));
      if strcmp(get(series(i),'selected'),'on')
         s=sprintf('%s <== selected',s);
      end
      disp(s);
   end   
   return;
end
ndx=args.serno<=length(series)&args.serno>0;
if any(ndx)
   curser = series(args.serno(ndx));
else
   error('GRIND:figdata:UnknownSeries','Series %s doesn''t exist',mat2str(args.serno));
end
X1 = getv(curser(1), 'xdata');
Y1 = getv(curser(1), 'ydata');
Z1 = getv(curser(1), 'zdata');
YZ1=[Y1,Z1];
for i=2:length(curser)
    YZ2=[getv(curser(i),'ydata'),getv(curser(i),'zdata')];
    X2=getv(curser(i),'xdata');
    [X1,YZ1]=i_concatdata(X1,YZ1,X2,YZ2);
end
A=[X1,YZ1];

function Ser=getv(h,dataax)
Ser=get(h,dataax);
if size(Ser,1)==1
   Ser=transpose(Ser);
end
