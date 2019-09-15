%DIRFIELD   Create direction field of a simple differential equation
%   Create a direction field for studying single first order
%   differential equations. You can click in the figure to run 
%   from that point.
%
%   Usage:
%   DIRFIELD - Create a direction field for 100 days and for the
%   variable which was selected for the X axis (see also <a href="matlab:help ax">ax</a>).
%   DIRFIELD VAR YLIM - Create a direction field for VAR with a
%   range of YLIM.
%   DIRFIELD('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'var' [state variable] - the state variable to analyse
%     'ylim' [number and length(number)==2] - the range of the state variable: two values [min max].
%
%  
%   Examples:
%   DIRFIELD X [0 10]
%
%   See also ru, time, null
%
%   Reference page in Help browser:
%      <a href="matlab:commands('dirfield')">commands dirfield</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function dirfield(varargin)
global g_grind t;
% avar=varargin{1};
% if nargin>1
%     alim=varargin{2};
% end
fieldnams={'var', 'v', 'the state variable to analyse';...
   'ylim', 'n&length(n)==2', 'the range of the state variable: two values [min max].'}';
args=i_parseargs(fieldnams,'if ~isempty(str2num(args{1})),deffields=''ylim'';else,deffields=''var,ylim'';end;','',varargin);
i_parcheck;
npoints = 35;
if g_grind.ndays > 100
   ndays = 100;
else
   ndays = g_grind.ndays;
end
if ~isfield(args,'var')
   args.var = g_grind.xaxis.var;
end
if ~isfield(args,'ylim')
   args.ylim = g_grind.xaxis.lim;
end
iY = i_getno(args.var);
if ~iY.isvar
   args.var = g_grind.yaxis.var;
   args.ylim = g_grind.yaxis.lim;
   iY = i_getno(args.var);
end
if isempty(iY.no)
   i_errordlg('Cannot create direction field if there are no state variables on the ''x''axis.');
   error('GRIND:dirfield:NoStateVariable','Cannot create direction field if there are no state variables on the ''x''axis.');
end
xlim1 = [t, t + ndays];
N = i_initvar;
[X,Y] = meshgrid(linspace(xlim1(1),xlim1(2),npoints),linspace(args.ylim(1) , args.ylim(2),npoints));
N0 = reshape(repmat(transpose(N), npoints * npoints, 1), npoints, npoints, g_grind.statevars.dim);
N0(:, :, iY.no) = Y;
Nres = i_runsinglestep(X, N0, true);
if g_grind.solver.isdiffer
    Nres=Nres(:,:,iY.no)-N0(:,:,iY.no);
end
[hfig, new] = i_makefig('dirfield');
if new
   set(hfig, 'WindowButtonDown', @(hobj,ev)i_callb('mdown2',hobj,args.var,ndays));
   set(hfig, 'WindowButtonMotionFcn', @(hobj,ev)i_callb('mmove',hobj));
   hold off;
end
set(hfig,'Name','Direction field');
oldhold = ishold;
%delta=ndays/npoints;
h = quiver(X, Y, ones(size(X)), Nres(:, :, iY.no));
i_plotdefaults;
hax=get(hfig,'currentaxes');
set(hax, 'ButtonDown', @(hobj,ev)i_callb('mdown2',hobj,args.var,ndays));
ud = get(hax, 'userdata');
ud.meta=struct('func','dirfield','xname','t','xlim',xlim1,'yname',args.var,'ylim',args.ylim,'zname','');  
set(hax, 'userdata', ud);
set(h, 'Color', [0.5 0.5 0.5]);
set(gca, 'XLim', xlim1);
set(gca, 'YLim', args.ylim);
set(gca, 'ZLim',[-1 1]);
xlabel('t');
ylabel(i_disptext(args.var));
if ~oldhold
   hold off;
end

