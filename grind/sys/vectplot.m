%VECTPLOT   special plots for vector state variables
%   Special function for vector state variable. Different ways of 
%   visualizing vector state variables (contour or surface plot 
%   with on the y-axis the reference number of the vector and on the x axis time,
%   or a movie.  A dialog screen (<a href="matlab:help replayall">replayall</a>) is opened to play the movie.
%   It is not possible to use this function for scalar variables.
%
%
%   Usage:
%   VECTPLOT - makes the current plots.
%   VECTPLOT NDAYS - runs the model for NDAYS days.
%   VECTPLOT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'ndays' [number>0] - number of days for the run
%   VECTPLOT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - add: continue with the previous run.
%     '-d' - set back to default output
%     '-l' - list the output
%     '-n' - nocheck: do not check for changed option, just show last results.
%     '-o' opens dialog box to select the plot.
%     '-o' NO [XAXIS YAXIS ZAXIS] ATYPE - the same but now with 3 axes.
%     '-o' N0 [XAXIS YAXIS] ATYPE - make a 2D plot/movie (if there 
%         is no time on the axes a movie is made). The axes can be any function of the
%         state variables or time t. ATYPE is type of plot (see below).
%     '-p' - makes it possible to replay a paranal
%     '-p1' or '-pa' - show the results of the last run of PARANAL.
%     '-p2' or '-paranal2' - show the results of the last two runs of PARANAL
%     '-r' - rerun the model always.
%     '-s' - silent (-s) do not plot the results.
%
%
%   ATYPE can be (see MATLAB help on the functions):
%   surface
%   pcolor
%   bar (or bar3)
%   stem (or stem3)
%   scatter
%   line or plot (or plot3)
%
%
%   See also time, model, viewcells, replayall
%
%   Reference page in Help browser:
%      <a href="matlab:commands('vectplot')">commands vectplot</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function vectplot(varargin)
%function vectplot(flag, no,xyzax, atype, varargin)
global g_grind t g_t g_Y;
%colno = 1;
%oldversion = 0;

if (nargin == 1) && strncmpi(varargin{1}, '-o', 2)
    i_parcheck;
   i_vectplotdlg;
   return;
elseif (nargin > 1) &&  strncmpi(varargin{1}, '-o', 2)
   i_parcheck;
   no = i_checkstr(varargin{2});
   g_grind.vectplot.currno = no;
   f=[strfind(varargin{3}, '[') strfind(varargin{3},']')];
   if length(f) == 2
      xyzax = varargin{3}(f(1) + 1:f(2) - 1);
   else
      xyzax = varargin{3};
   end
   f = strfind(xyzax,',');
   if isempty(f)
      f = strfind(xyzax, ' ');
   end
   if isempty(f)
      g_grind.vectplot.vars{no}.xax = sprintf('cellno(%s)', xyzax);
      g_grind.vectplot.vars{no}.yax = xyzax;
      g_grind.vectplot.vars{no}.zax = '';
   elseif length(f) == 1
      g_grind.vectplot.vars{no}.xax = xyzax(1:f - 1);
      g_grind.vectplot.vars{no}.yax = xyzax(f + 1:end);
      g_grind.vectplot.vars{no}.zax = '';
   else
      g_grind.vectplot.vars{no}.xax = xyzax(1:f(1) - 1);
      g_grind.vectplot.vars{no}.yax = xyzax(f(1) + 1:f(2) - 1);
      g_grind.vectplot.vars{no}.zax = xyzax(f(2) + 1:end);
   end
   if nargin > 3
      atype  = varargin{4};
      if strncmpi(atype, 'su', 2)
         g_grind.vectplot.vars{no}.type = 'surface';
      elseif strncmpi(atype, 'pc', 2)
         g_grind.vectplot.vars{no}.type = 'pcolor';
      elseif strncmpi(atype, 'mc', 2)
         g_grind.vectplot.vars{no}.type = 'mcolor';
      elseif strncmpi(atype, 'ba', 2)
         g_grind.vectplot.vars{no}.type = 'bar';
      elseif strncmpi(atype, 'st', 2)
         g_grind.vectplot.vars{no}.type = 'stem';
      elseif strncmpi(atype, 'sc', 2)
         g_grind.vectplot.vars{no}.type = 'scatter';
      elseif strncmpi(atype, 'l', 1)||strncmpi(atype, 'pl', 2)
         g_grind.vectplot.vars{no}.type = 'plot';
      else
         g_grind.vectplot.vars{no}.type = lower(atype);
      end
      if strcontains(atype, '3')&&isempty(strfind(g_grind.vectplot.vars{no}.type,'3'))
         g_grind.vectplot.vars{no}.type = [g_grind.vectplot.vars{no}.type '3'];
      end
   else
      g_grind.vectplot.vars{no}.type = 'pcolor';
   end
   if length(varargin) > 4
      g_grind.vectplot.vars{no}.axsettings = varargin(5:end);
   end
   return;
end
fieldnams={'ndays', 'n>0', 'number of days for the run',g_grind.ndays}';
res = i_time_options(i_time_options(i_parseargs(fieldnams,'ndays',...
    {'-r','-a','-n','-s','-p2|-paranal2','-p1|-pa','-d','-l','-p','-o'},varargin)));

i_parcheck;
%type = 0;
%log1 = 0;
%xa = '';
%ya = '';
for i = 1:length(res.opts)
   if strncmpi(res.opts(i), '-d', 2) %-defaults
      if ~isempty(g_grind.statevars.vectnames)
         if ~isfield(g_grind, 'vectplot')
            g_grind.vectplot.currno = 1;
            no = 1;
            for k = 1:length(g_grind.statevars.vectnames)
               if g_grind.statevars.dims{k}.dim1 * g_grind.statevars.dims{k}.dim2 > 1
                  g_grind.vectplot.vars{no}.xax = sprintf('cellno(%s)', g_grind.statevars.vectnames{k});
                  g_grind.vectplot.vars{no}.yax = g_grind.statevars.vectnames{k};
                  g_grind.vectplot.vars{no}.zax = '';
                  g_grind.vectplot.vars{no}.type = 'stem';
                  no = no + 1;
               end
            end
         end
      else
         if ~isfield(g_grind, 'vectplot')
            g_grind.vectplot.currno = 1;
            no = 1;
            for k = 1:length(g_grind.statevars.names)
               g_grind.vectplot.vars{no}.xax = 't';
               g_grind.vectplot.vars{no}.yax = g_grind.statevars.names{k};
               g_grind.vectplot.vars{no}.zax = '';
               g_grind.vectplot.vars{no}.type = 'plot';
               no = no + 1;
            end
         end
         warning('GRIND:vectplot:NoVectors','Vectplot is designed to be used for vector state variables');
      end
      if nargin ~= 0
         return;
      end
   end
   if  strncmpi(res.opts(i), '-l', 2)
      displayout;
      return;
   elseif  strncmpi(res.opts(i), '-p', 2)
      addparanalreplaydata;
      return;
   end
end
if ~isfield(g_grind, 'vectplot')
   vectplot('-d');
end
for no = length(g_grind.vectplot.vars):-1:1
   if all([isempty(g_grind.vectplot.vars{no}.xax), isempty(g_grind.vectplot.vars{no}.yax), isempty(g_grind.vectplot.vars{no}.xax)])
      g_grind.vectplot.vars(no) = [];
   end
end
if res.rerun
   %     disp('running');
   i_ru(t, res.ndays, res.N0, 1);
end
if ~res.silent
   makevectplot;
end
if res.adding
   g_grind.solver.addmode = false;
end
if ~isempty(res.OldY)
   g_Y = res.OldY;
   g_t = res.Oldt;
end
function makevectplot
global g_grind;
for no = 1:length(g_grind.vectplot.vars)
   xax = i_getno(g_grind.vectplot.vars{no}.xax);
   xax.name = g_grind.vectplot.vars{no}.xax;
   xax.t = strcmp(xax.name, 't');
   yax = i_getno(g_grind.vectplot.vars{no}.yax);
   yax.name = g_grind.vectplot.vars{no}.yax;
   yax.t = strcmp(yax.name, 't');
   %if yasis is time then we exchange the axis with the x-axis and use a
   %rotation to put the x-asis on the y-axis
   zax = i_getno(g_grind.vectplot.vars{no}.zax);
   zax.name = g_grind.vectplot.vars{no}.zax;
   zax.t = strcmp(zax.name, 't');
   % replay is not possible
   type = g_grind.vectplot.vars{no}.type;
   if strcmp(yax.name, 't')
      exchangexy = 1;
      h = xax;
      xax = yax;
      yax = h;
   else
      exchangexy = 0;
   end
   
   
   X = getvalues(xax);
   Y = getvalues(yax);
   Z = getvalues(zax);
   if any(~isreal(X))||any(~isreal(Y))||any(~isreal(Z))
       X=real(X);
       Y=real(Y);
       Z=real(Z);
       warning('grind:vectplot:complex','Ignored imaginary parts of complex values');
   end
   if isempty(Z)
      if any(strcmp(type,{'plot3','scatter3','stem3'}))
         type = type(1:end - 1);
      end
   else
      if any(strcmp(type,{'plot','scatter','stem'}))
         type = sprintf('%s3', type);
      end
   end
   notime =  ~xax.t && ~yax.t && ~zax.t;
   
   %i = 0;
   %H = i_figno('vectplot') + i;
   %while ishandle(H + i)
   %   i = i + 1;
   %end
   hfig = i_makefig('vectplot', no);
   set(hfig, 'Name', 'Vector plot');
   set(hfig,'DoubleBuffer','on');
   if ~isoctave&&verLessThan('matlab','8.6')
      set(hfig,'renderer','ZBuffer'); %in fact this is solving a window bug
   end
   set(gca,'XLimMode','auto','YLimMode','auto','ZLimMode','auto','CLimMode','auto');
   %   if ~oldversion
   if notime
      X1 = X;
      Y1 = Y;
      Z1 = Z;
      X = minmax(X1);
      Y = minmax(Y1);
      if ~isempty(Z)
         Z = minmax(Z1);
      end
   end
   if (size(Y, 2) > 1) && (size(X ,2)==1)
      X = repmat(X, 1, size(Y, 2));
   end
   if (size(X, 2) > 1) && (size(Y ,2)==1)
      Y = repmat(Y, 1, size(X, 2));
   end
   if isempty(Z)
      if any(strcmp({'surface','mcolor','pcolor','bar3'},type))
         Z = Y;
         if size(Y, 1) == 1
            Z = repmat(Y, 3, 1);
            X  = repmat(X, 3, 1);
         end
         if size(Y, 2) == 1
            Z = repmat(Y, 1, 3);
            X = repmat(X, 3, 1);
         end
         Y = repmat(1:size(X, 2), size(X, 1), 1);
         h = feval(str2func(type),transpose(X),transpose(Y), transpose(Z));
         i_plotdefaults(hfig);
         set(hfig, 'colormap', i_grindcolormap);
         hbar = colorbar;
         ylabel(hbar, yax.name);
         i_plotdefaults(hfig);
         set(h,'tag','vectplot');
         shading flat;
         ylabel(sprintf('Cell number %s', yax.name));
         xlabel(xax.name);
         zlabel(yax.name);
      else
         h = feval(str2func(type), X, Y);
         i_plotdefaults(hfig);
         set(h,'tag','vectplot');
         xlabel(xax.name);
         ylabel(yax.name);
      end
   else
      %       if any(strcmp({'scatter3','plot3'},type))
      %           X=X(:);
      %           Y=Y(:);
      %           Z=Z(:);
      %       end
      h = feval(str2func(type), X, Y, Z);
      set(h,'tag','vectplot');
      xlabel(xax.name);
      ylabel(yax.name);
      zlabel(zax.name);
   end
   if exchangexy
      if isempty(Z)||any(strcmp({'mcolor','pcolor'},type))
         view([  - 90, 90]);
      else
         view([322.5 - 90 30])
      end
      set(gca,'YDir','reverse')
   end
   if ~zax.t
      addreplay(hfig, no)
   end
   if notime
      ser = findobj(gca, 'tag', 'vectplot');
      set(gca,'XLimMode','manual','YLimMode','manual','ZLimMode','manual','CLimMode','manual');
      if ~isempty(Z1)
         set(ser, 'XData', X1(1,:), 'YData',Y1(1,:), 'ZData',Z1(1,:));
      elseif any(strcmp({'surface','pcolor','mcolor','bar3'},type))
         Y1 = Y1(1, :);
         if size(Y1, 1) == 1
            Y1 = repmat(Y1, 3, 1);
         end
         if size(Y1, 2) == 1
            Y1 = repmat(Y1, 1, 3);
         end
         if ~any(strcmp({'pcolor','mcolor'}, type))
            set(ser,'CData', transpose(Y1),'Zdata', transpose(Y1));
         elseif strcmp('mcolor', type)
            Y1 = transpose(Y1);
            Y1(end + 1, end + 1) = 0;
            set(ser, 'CData', Y1, 'Zdata', zeros(size(Y1)));
         else
            set(ser, 'CData', transpose(Y1), 'Zdata', zeros(size(transpose(Y))));
         end
      else
         set(ser, 'XData', X1(1,:), 'YData',Y1(1,:));
      end
      replayall;
   end
   if isfield(g_grind.vectplot.vars{no}, 'axsettings')
      set(gca, g_grind.vectplot.vars{no}.axsettings{:});
   end
end

disp('use <a href="matlab:vectplot -out">vectplot -out</a> to change the output/add plots');
function X = minmax(X1)
X = max(X1, [], 1);
absmin = min(min(X1));
minx = min(X);
if (length(X) > 1)&&absmin < minx
   i=find(X == minx, 1);
   X(i) = absmin;
end
function h = surface(varargin) %#ok<DEFNU>
h = surf( varargin{:});
function h = scatter3(X, Y, Z) %#ok<DEFNU>
h = plot3(X, Y, Z, '.');
function h = scatter(X, Y) %#ok<DEFNU>
h = plot(X, Y, '.');
function X = getvalues(aax)
if isempty(aax.name)
   X = [];
   return
else
   X = outfun(aax.name);
end
function displayout
global g_grind;
for no = 1:length(g_grind.vectplot.vars)
   if isempty(g_grind.vectplot.vars{no}.zax)
      fprintf('vectplot -out %d [%s %s] %s\n', no, g_grind.vectplot.vars{no}.xax, ...
         g_grind.vectplot.vars{no}.yax, g_grind.vectplot.vars{no}.type);
   else
      fprintf('vectplot -out %d [%s %s %s] %s\n', no, g_grind.vectplot.vars{no}.xax, ...
         g_grind.vectplot.vars{no}.yax, g_grind.vectplot.vars{no}.zax, g_grind.vectplot.vars{no}.type);
   end
end

function addreplay(H, No)
global g_t;
hax=findobj(H,'type','axes');
tags = get(hax, 'tag');
hax=hax(~(strcmp(tags,'legend')|strcmp(tags,'Colorbar')));
ud = get(hax, 'userdata');
ud.replay.callback = @replaycallback;
ud.replay.onstart = @replaystart;
ud.replay.onend = @replayend;
if ~isempty(g_t)
   ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
end
ud.replay.opt = No;
set(hax,'userdata', ud);

function replaystart(hax)
global g_grind g_t;
if ishandle(hax)
   %    N0 = i_initvar;
   %    if i_settingschanged(N0, g_grind.ndays)
   %       i_ru(t, g_grind.ndays, N0, 1);
   %    end
   hfig = get(hax, 'parent');
   i_figure(hfig);
   set(hax,'XLimMode','manual','YLimMode','manual','ZLimMode','manual','CLimMode','manual');
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      ud.replay.settings=struct('tvar','t','tlim',[g_t(1) g_t(end)],'numt',length(g_t));
      no = ud.replay.opt;
      if no <= length(g_grind.vectplot.vars)
         xax = i_getno(g_grind.vectplot.vars{no}.xax);
         xax.name = g_grind.vectplot.vars{no}.xax;
         xax.t = strcmp(xax.name, 't');
         yax = i_getno(g_grind.vectplot.vars{no}.yax);
         yax.name = g_grind.vectplot.vars{no}.yax;
         yax.t = strcmp(yax.name, 't');
         zax = i_getno(g_grind.vectplot.vars{no}.zax);
         zax.name = g_grind.vectplot.vars{no}.zax;
         zax.t = strcmp(zax.name, 't');
         ud.replay.xdata = getvalues(xax);
         ud.replay.ydata = getvalues(yax);
         ud.replay.zdata = getvalues(zax);
         type = g_grind.vectplot.vars{no}.type;
         if ~isoctave&&(verLessThan('matlab','8.4.0')&&strcmp(type, 'stem')) %fixing a bug in MATLAB
            ser = get(hax, 'children');
            ser1 = findobj(hax, 'tag', 'vectplot');
            ser2  = get(ser1, 'baseline');
            if ~any(ser2 == ser)
               ud1 = get(ser1, 'userdata');
               v = get(hax, 'view');
               delete(ser);
               if isempty(ud.replay.zdata)
                  h = stem(hax, ud.replay.xdata(end, :), ud.replay.ydata(end, :));
               else
                  h = stem3(hax, ud.replay.xdata(end, :), ud.replay.ydata(end, :), ud.replay.zdata(end, :));
               end
               i_plotdefaults(hfig)
               set(hax, 'view', v);
               set(h,'userdata',ud1,'tag','vectplot');
            end
         end
         if (size(ud.replay.ydata, 2) > 1) && (size(ud.replay.xdata ,2)==1)
            ud.replay.xdata = repmat(ud.replay.xdata, 1, size(ud.replay.ydata, 2));
         end
         if (size(ud.replay.xdata, 2) > 1) && (size(ud.replay.ydata ,2)==1)
            ud.replay.ydata = repmat(ud.replay.ydata, 1, size(ud.replay.xdata, 2));
         end
         if any(~isreal(ud.replay.xdata(:)))||any(~isreal(ud.replay.ydata(:)))
             ud.replay.xdata=real(ud.replay.xdata);
             ud.replay.ydata=real(ud.replay.ydata);
             warning('grind:vecplot:complex','Ignored imaginary parts of complex values')
         end
         if xax.t
            ud.replay.timeax = 1;
         elseif yax.t
            ud.replay.timeax = 1;
            h = ud.replay.xdata;
            ud.replay.xdata = ud.replay.ydata;
            ud.replay.ydata = h;
         else
            ud.replay.timeax = [];
         end
         set(hax, 'userdata', ud);
         if ~isoctave&&verLessThan('matlab','8.4.0')
            set(hax, 'drawmode','fast');
         else
            set(hax,'SortMethod','depth');
         end
      end
   end
end

function replayend(hax, closedlg)
if closedlg
   replaycallback(hax, 't', 1);
end
if ishandle(hax)
   i_figure(get(hax, 'parent'));
   %    ax =  findobj(hax,'type','axes');
   %    tags = get(ax, 'tag');
   %    ax = ax(~strcmp(tags, 'legend'));
   ud = get(hax, 'userdata');
   %  set(ax,'XLimMode','auto','YLimMode','auto','ZLimMode','auto');
   ud.replay.xdata  = [];
   ud.replay.ydata = [];
   ud.replay.zdata = [];
   set(hax, 'userdata', ud);
end
function t = replaycallback(hax,avar, relt)
global g_t;
t = [];
if ishandle(hax)&&isempty(avar )||strcmp(avar, 't')
   ud = get(hax, 'userdata');
   if isfield(ud.replay, 'timeax') %replaystart must have been run
      t = ud.replay.settings.tlim(1) + relt * (ud.replay.settings.tlim(end) - ud.replay.settings.tlim(1));
      [~, ndx] = min(abs(g_t - t));
      t = g_t(ndx);
      if ~isempty(ud.replay.timeax)
         ndx=g_t <= t;
      end
      updatevectplot(hax, ndx, ud);
   end
end
function updatevectplot(hax, ndx, ud)
global g_grind;
no = ud.replay.opt;
type = g_grind.vectplot.vars{no}.type;
if isfield(ud, 'replay')
   ser = findobj(hax, 'tag', 'vectplot');
   if length(ser) > 1
      for i = 1:length(ser)
         if ~isempty(ud.replay.zdata)
            set(ser(i), 'XData', ud.replay.xdata(ndx,i), 'YData',ud.replay.ydata(ndx,i), 'ZData',ud.replay.zdata(ndx,i));
         else
            set(ser(i), 'XData', ud.replay.xdata(ndx,i), 'YData',ud.replay.ydata(ndx,i))
         end
      end
      return;
   end
   if any(strcmp({'surface','pcolor','mcolor','bar3'},type))
      Z = ud.replay.ydata;
      if isempty(ud.replay.timeax)
         Z = Z(ndx, :);
      else
         Z(~ndx, :) = NaN;
      end
      if size(Z, 1) == 1
         Z = repmat(Z, 3, 1);
      end
      if size(Z, 2) == 1
         Z = repmat(Z, 1, 3);
      end
      %       if ~any(strcmp({'mcolor','pcolor'}, type))
      %          set(ser,'CData', Y','Zdata', Y');
      %       else
      Z = transpose(Z);
      C = Z;
      X = get(ser, 'Xdata');
      Y = get(ser, 'Ydata');
      if strcmp('mcolor', type)
         C(end + 1, end + 1) = 0;
         Z = zeros(size(C));
      elseif strcmp('pcolor', type)
         Z = zeros(size(C));
      else
         Z = C;
      end
      set(ser,'Xdata',X,'Ydata',Y, 'Zdata',Z, 'CData', C);
      %     end
      return;
   else
      X = ud.replay.xdata(ndx, :);
      Y = ud.replay.ydata(ndx, :);
   end
   makevect= any(strcmp({'stem','stem3','scatter3'},type));
   if makevect
      X = X(:);
      Y = Y(:);
   end
   if ~isempty(ud.replay.zdata)
      Z = ud.replay.zdata(ndx, :);
      if makevect
         Z = Z(:);
      end
      set(ser, 'XData', X, 'YData',Y, 'ZData',Z);
   else
      set(ser, 'XData', X, 'YData',Y);
   end
end
function addparanalreplaydata
global g_paranal g_grind;
if isempty(g_paranal.run)
   error('grind:time:paranal','Paranal should be run before using this option');
end
vectplot();
for No = 1:size(g_grind.vectplot.vars, 2)
   if ~isempty(g_grind.vectplot.vars{No})
      [h] = i_figno('vectplot') +  No;
      hax=findobj(h,'type','axes');
      tags = get(hax, 'tag');
      hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
      ud = get(hax, 'userdata');
      ud.replay.callback = @paranalreplaycallback;
      ud.replay.onstart = @paranalreplaystart;
      ud.replay.onend = @paranalreplayend;
      ud.replay.onturn = @i_replayparanalturn;
      ud.replay.pars = g_paranal.run.parvalues(1, :);
      tstep = (ud.replay.pars(end) - ud.replay.pars(1)) / length(ud.replay.pars);
      ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',numel(g_paranal.run.t));
      ud.replay.opt = No;
      set(hax,'userdata', ud);
   end
end
replayall('variable', g_paranal.run.pars{1});

function paranalreplaystart(hax)
global g_grind g_paranal;
if ishandle(hax)
   i_figure(get(hax, 'parent'));
   set(hax,'xlimmode','auto');
   set(hax,'ylimmode','auto');
   set(hax,'zlimmode','auto');
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      ud.replay.pars = g_paranal.run.parvalues(:, 1);
      tstep = (ud.replay.pars(end) - ud.replay.pars(1)) / length(ud.replay.pars);
      ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',numel(g_paranal.run.t));
      No = ud.replay.opt;
      ud.replay.ydata = [];
      ud.replay.xdata = [];
      ud.replay.zdata = [];
      if No > length(g_grind.vectplot.vars)
         No = length(g_grind.vectplot.vars);
      end
      xvar = g_grind.vectplot.vars{No}.xax;
      yvar = g_grind.vectplot.vars{No}.yax;
      zvar = g_grind.vectplot.vars{No}.zax;
      
      ud.replay.xdata = outfun(xvar, '-p');
      if ~isempty(yvar)
         ud.replay.ydata = outfun(yvar, '-p');
      end
      if isempty(zvar)
         ud.replay.zdata = outfun(zvar, '-p');
      end
      if strcmp(xvar, 't')
         ud.replay.timeax = 1;
      elseif strcmp(yvar, 't')
         ud.replay.timeax = 1;
         h = ud.replay.xdata;
         ud.replay.xdata = ud.replay.ydata;
         ud.replay.ydata = h;
      else
         ud.replay.timeax = [];
      end
   end
   set(hax, 'userdata', ud);
   %   set(hax,'drawmode','fast');
   %    if ~isempty(g_paranal.run.t)
   %       tt1=g_paranal.run.t(g_paranal.run.p(1) == g_paranal.run.p);
   %       xlim(hax,sort([0 tt1(end) - tt1(1)]));
   %       ylim(hax,sort([0 maxy]));
   %
   %       set(hax,'xlimmode','manual');
   %    end
end

function paranalreplayend(hax, closedlg)
try
   if closedlg
      replaycallback(hax, '', 1);
   end
catch
end
if ishandle(hax)
   i_figure(get(hax, 'parent'));
   %    set(hax,'xlimmode','auto');
   %    set(hax,'ylimmode','auto');
   %  set(hax,'zlimmode','auto');
   ud = get(hax, 'userdata');
   ud.replay.xdata = [];
   ud.replay.ydata = [];
   ud.replay.zdata  = [];
   set(hax, 'userdata', ud);
end

function p = paranalreplaycallback(hax,avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&isempty(avar)&&(~isempty(g_paranal.run)&&strcmp(avar, g_paranal.run.pars{1}))
   ud = get(hax, 'userdata');
   if isfield(ud, 'replay')
      ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
      [~, ~, stepndx] = ind2sub(size(g_paranal.run.t), ndx1);
      p = g_paranal.run.parvalues(stepndx);
      if ~isempty(ud.replay.timeax)
         ndx1 =g_paranal.run.t <= g_paranal.run.t(ndx1);
      end
      updatevectplot(hax, ndx1, ud)
   end
end
