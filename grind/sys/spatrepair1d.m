function spatrepair1d(err, distloc, refloc, disturbfunc, apar, p_start, p_end, npoints)
global g_grind;
if ~g_grind.statevars.vector
   error('GRIND:spatrepair1d:NoVector','designed for vector/matrix variables only')
end
if ~((nargin==1) && (i_checkstr(err)==-1))
   if nargin < 8
      %vertraging voor de diffusion factor
      centr = round([g_grind.statevars.dims{1}.dim1, g_grind.statevars.dims{1}.dim2] ./ 2);
      if isfield(g_grind, 'spatrepair1d')
         answer = g_grind.spatrepair1d;
      else
         answer={'0.01' sprintf('[%g,%g]',centr) '[1,1]' '','','0','10','50'};
      end
      prompt={'max allowed difference','location of disturbance','reference',...
         'disturbance function','parameter','start','end','nsteps'};
      answer = inputdlg(prompt,'spatial repair 1D', 1,answer);
      g_grind.spatrepair1d = answer;
      err = str2double(answer{1});
      distloc = str2num(answer{2});  %#ok<ST2NM>
      refloc = str2num(answer{3});  %#ok<ST2NM>
      disturbfunc = answer{4};
      apar = answer{5};
      p_start = str2double(answer{6});
      p_end = str2double(answer{7});
      npoints = str2double(answer{8});
   end
   params = p_start:(p_end - p_start) / (npoints - 1):p_end;
   isrepairing = zeros(size(params));
   times=zeros(1,length(params));
   for i = 1:length(params)
      assignin('base', apar, params(i));
      evalin('base', disturbfunc);
      simtime 0 500 500
      [Reptime, Repairing] = spatrepair(err,distloc,refloc,10);
      isrepairing(i) = Repairing;
      times(i) = Reptime;
   end
   titl = apar;
   save tmp isrepairing times params titl;
end
load tmp isrepairing times params titl
hfig=figure(1);
plot(params, times);
i_plotdefaults(hfig);
hold on;
xlabel(titl);
ylabel('Repair time (d)');
ylims = get(gca, 'YLim');
set(gca, 'Xlim', [min(params) max(params)]);
maxim = ylims(2);
isreturning = ~isnan(isrepairing) & ~isrepairing;
plotyesno(isreturning, params, maxim, [1 0 0]);
ispattern = isinf(times);
plotyesno(ispattern, params, maxim, [0 1 0]);
istoolong = isnan(times) & ~isnan(isrepairing) & isrepairing;
plotyesno(istoolong, params, maxim, [0.8 0.8 0.8]);
%[params', times',isrepairing']
hold off;
%saveas(gcf,sprintf('fig%s.fig',figname),'fig');
%saveas(gcf,sprintf('fig%s.ai',figname),'ai');

function plotyesno(yndata, pars, maxY, color)
f=find(yndata == 1);
i = 2;
deltapar = (pars(2) - pars(1)) / 2;
if length(f) == 1
   aX = [pars(f(1)) - deltapar pars(f(1)) - deltapar pars(f(1)) + deltapar pars(f(1)) + deltapar];
   aY = [0 maxY maxY 0];
   area(aX, aY, 'facecolor', color);
else
   while i < length(f)
      istart = f(i - 1);
      while f(i) == f(i - 1) + 1 && (i < length(f))
         i = i + 1;
      end
      iend = f(i);
      aX = [pars(istart) - deltapar pars(istart) - deltapar pars(iend) + deltapar pars(iend) + deltapar];
      aY = [0 maxY maxY 0];
      area(aX, aY, 'facecolor', color);
   end
end
