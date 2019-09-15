%QUASIPOT   Create quasi potentials for higher dimensional systems
%    It is only possible to calculate potential landscapes for "gradient systems". All 
%    1 equation systems are gradient. The 2D system:
%    dX/dt=f(x,y)
%    dY/dt=g(x,y)
%    is only "gradient" if
%    d(f(x,y))/dy=d(g(x,y))/dx
%   
%    This is a very restrictive condition that in practice only occurs in symmetric models. So most models are non-gradient.
%    Therefore we use an approximation, the quasi-potential (Vq). We use here two algorithms.
%    The first one (Rodríguez-Sánchez, et al. in prep) is a simple fast method where you decompose the Jacobian in two parts,
%    the gradient part and the curl part. The landscape is constructed with the gradient part, the curl is interpreted as error.
%    This method is exact for gradient systems.
%    The other method is the slower method of Bhattacharya, et al. 2011. To calculate this we solve the following equation
%    in a grid of initial conditions:
%    d(Vq)/dt=-((d(f(x,y))/dt)^2 + (d(g(x,y))/dt)^2).
%    Subsequently all trajectories are aligned and the difference in potential between
%    the equilibria is determined by comparing the potentials of trajectories close to the
%    separatrix.
%   
%    Algorithms:
%    Rodríguez-Sánchez, P., Scheffer, M. & E.H. van Nes (in prep.) Climbing Escher’s ladder, or how to 
%       compute stability landscapes for weakly non-gradient systems. 
%    Bhattacharya, S., Q. Zhang, and M. E. Andersen. 2011. A deterministic map of Waddington's
%      epigenetic landscape for cell fate specification. BMC Systems Biology 5:1-12.
%
%   Usage:
%   QUASIPOT [NPOINTS NPOINTS] - number of grid points in x and y direction (default [300 300] (or [20 20] for bhattacharya))
%   QUASIPOT NPOINTS - run for NDAYS in a grid of NPOINTSxNPOINTS     
%   res=QUASIPOT(..) - saves results in a structure with fields res.X, res.Y and res.pot (grid values of the quasi-potential).
%   QUASIPOT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'algorithm' [bhattacharya | rodriguez] - The kind of algorithm used (default rodriguez).
%     'init_sep' [number] - values of the separatrix (only for bhattacharya)
%     'maxpot' [number] - maximum value of the potential in the plot (default empty)
%     'ndays' [number>0] - max number of time units per run (only for bhattacharya).
%     'ndays_sep' [integer>0] - number of timesteps for separatrix (only for bhattacharya)
%     'npoints' [integer>0 and length(integer)<=2] - number of grid points for each axis (can be a 2 value vector)
%     'silent' [logical] - set to true to suppress figures
%     'tstep' [integer>0] - number of output values per run
%   QUASIPOT('-opt1','-opt2',...) - Valid command line options:
%     '-n' - suppress recalculation of the quasi potentials, even if the parameters have changed.
%     '-r' - recalculate the quasi potentials, even if the parameters have not changed.
%     '-s' - silent mode, figures suppressed.
%     
%
%   See also potential, marbleplot, fokkerplanck
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('quasipot')">commands quasipot</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function res = quasipot(varargin)
global g_grind g_quasipot;
if isempty(g_quasipot)
    evalin('base','global g_quasipot');
end

%number of gridcells
fieldnams={'ndays', 'n>0', 'max number of time units per run (only for bhattacharya).',g_grind.ndays;...
   'npoints', 'i>0&length(i)<=2', 'number of grid points for each axis (can be a 2 value vector)',[300 300];...
   'tstep', 'i>0', 'number of output values per run',g_grind.tstep;...
   'maxpot', 'n', 'maximum value of the potential in the plot (default empty)',[];...
   'silent', 'l', 'set to true to suppress figures', false;...
   'init_sep', 'n', 'values of the separatrix (only for bhattacharya)',[];...
   'ndays_sep', 'i>0', 'number of timesteps for separatrix (only for bhattacharya)',[];...
   'algorithm', 'e[b+hattacharya|r+odriguez]', 'The kind of algorithm used (default rodriguez).','rodriguez'}';
args=i_parseargs(fieldnams,'npoints,ndays,init_sep',{'-r','-s','-n'},varargin);
%default arguments


if ~isfield(args, 'maxpot')
    args.maxpot = nan;
end
if ~isfield(args, 'algorithm')
    args.algorithm='rodriguez';
end

%number of time steps per run
if ~isfield(args, 'ndays')
    args.ndays = g_grind.ndays;
end

if ~isfield(args, 'npoints')
    if strcmp(args.algorithm,'rodriguez')
        args.npoints=[300 300];
    else
        args.npoints = [20 20];
    end
end

if numel(args.npoints) == 1
    args.npoints = [args.npoints args.npoints];
end

if ~isfield(args, 'init_sep')
    args.init_sep=[];
    if ~isfield(args,'ndays_sep')
        args.ndays_sep=args.ndays*2;
    end
else
    error('grind:quasipot','Analysing user-defined separatrices not yet implemented');
end


if ~isfield(args, 'tstep')
    args.tstep = g_grind.tstep;
end

if any(strcmp(args.opts, '-s'))
    args.silent = true;
end

if ~isfield(args, 'silent')
    args.silent = false;
end

if isnan(args.tstep)
    args.tstep = 1000;
end

if strcmp(args.algorithm,'rodriguez')
    res1=i_quasipot_rodriguez(args);
    if nargout>0
        res=res1;
    end
    doplot(args,res1.X,res1.Y,res1.pot,res1.curl);
    return;
end

opts = args.opts;
args = rmfield(args, 'opts');
g_quasipot.run.pars = {g_grind.xaxis.var, g_grind.yaxis.var};
if ~any(strcmp(opts,'-n'))&&settingschanged(args)||any(strcmp(opts, '-r'))
    g_quasipot.args = args;
    [x, y] = meshgrid(linspace(g_grind.xaxis.lim(1) + 0.001, g_grind.xaxis.lim(2), args.npoints(1)), ...
        linspace(g_grind.yaxis.lim(1) + 0.001, g_grind.yaxis.lim(2), args.npoints(2)));
    g_quasipot.run.parvalues = [x(:) y(:)];
    g_quasipot.run.Y = zeros(args.tstep, g_grind.statevars.dim, size(g_quasipot.run.parvalues, 1));
    g_quasipot.run.t = transpose(linspace(0, args.ndays, args.tstep));
    g_quasipot.run.pot = zeros(args.tstep, 1, size(g_quasipot.run.parvalues, 1));
    g_quasipot.run_separ.parvalues= args.init_sep;
    if ~isempty(g_quasipot.run_separ.parvalues)
        g_quasipot.run_separ.Y = zeros(args.tstep, g_grind.statevars.dim, size(g_quasipot.run.parvalues, 1));
        g_quasipot.run_separ.t = transpose(linspace(0, ndays_sep, args.tstep));
        g_quasipot.run_separ.pot = zeros(args.tstep, 1, size(g_quasipot.run.parvalues, 1));
    end
    
    doupdate(args.silent)
end


nruns = size(g_quasipot.run.parvalues, 1);
%Find the attractors
endy = g_quasipot.run.Y(end, :, :);
iX = i_getno(g_quasipot.run.pars{1});
iY = i_getno(g_quasipot.run.pars{2});

endy = reshape(permute(endy, [3, 2, 1]), nruns, g_grind.statevars.dim);
[~, ndx] = sortrows(endy);

%Distance between sorted endpoints
dist = sum(diff(endy(ndx,:)).^2, 2);
%if the distance is >twice the mean, different attractor
ndx1 = [0; find(dist > 2 * mean(dist)); nruns];
attract = zeros(nruns, 1);
uattract = zeros(length(ndx1) - 1, g_grind.statevars.dim);
pattract = zeros(size(uattract, 1), 1);
for i = 1:length(ndx1) - 1
    attract(ndx(ndx1(i) + 1:ndx1(i + 1))) = i;
    uattract(i, :) = endy(ndx(ndx1(i) + 1), :);
end

attract = reshape(attract, g_quasipot.args.npoints(1), g_quasipot.args.npoints(2));
x = reshape(g_quasipot.run.parvalues(:, 1), g_quasipot.args.npoints(1), g_quasipot.args.npoints(2));
y = reshape(g_quasipot.run.parvalues(:, 2), g_quasipot.args.npoints(1), g_quasipot.args.npoints(2));
% figure(10);
% pcolor(x,y,attract);

%Find pairs of points on separatrix
sep1=[zeros(1, size(attract, 2)); diff(attract) ~= 0] > 0;%vertical separatrix point 1
sep2=[diff(attract) ~= 0; zeros(1, size(attract, 2))] > 0;%vertical separatrix point 2
sep3=[zeros(1, size(attract, 2)); diff(attract')~=0]' > 0;%horizontal separatrix point 1
sep4=[diff(attract')~=0;zeros(1,size(attract,2))]' > 0;%horizontal separatrix point 2
sep1 = [find(sep1); find(sep3)];
sep2 = [find(sep2); find(sep4)];
g_quasipot.seppairs = [sep1 sep2];
attrpair = reshape(attract(g_quasipot.seppairs(:)), size(g_quasipot.seppairs));
endpot = g_quasipot.run.pot(end, 1, :);
endpot = reshape(permute(endpot, [3, 2, 1]), nruns, 1);
potdiff = reshape(endpot(g_quasipot.seppairs(:)), size(g_quasipot.seppairs));
potdiff = (potdiff(:, 2) - potdiff(:, 1)) .* sign(attrpair(:, 2) - attrpair(:, 1));
attrpair1 = sort(attrpair, 2);
upairs = unique(attrpair1, 'rows');
upairs(end+1,1)=upairs(1,1);
upairs(end,2)=upairs(end-1,2);
%[~,ndx2]=sort(rand(size(upairs,1),1));
%weight=1;
for i = 1:size(upairs, 1)
    ndx1=all(attrpair1 == repmat(upairs(i, :), size(attrpair1, 1), 1), 2);
    pattract(upairs(i, 2)) = pattract(upairs(i, 1)) + median(potdiff(ndx1));
end

pattract=pattract(2:end);
[pattract,err]=fminsearch(@(pattract)optimizeatt(pattract,g_quasipot,attract),pattract);
pattract=[0;pattract];

for i = 1:length(pattract)
    fprintf('Attractor %d:\n', i);
    for j = 1:g_grind.statevars.dim
        fprintf('%s = %g\n', i_statevars_names(j), uattract(i, j));
    end
    
    fprintf('Quasi-potential = %g\n\n', pattract(i));
end

fprintf('\nError at separatrices: %g\n',err);
pot = g_quasipot.run.pot;
for i = 1:size(pot, 3)
    pot(:, 1, i) = pot(:, 1, i) - pot(end, 1, i);
end

for i = 1:length(pattract)
    pot(:, :, attract(:) == i)=pot(:, :, attract(:) == i) + pattract(i);
end


%figure;
y1 = g_quasipot.run.Y(:, iX.no, :);
y2 = g_quasipot.run.Y(:, iY.no, :);
pot2 = pot(:, 1, :);
x = reshape(g_quasipot.run.parvalues(:, 1), g_quasipot.args.npoints(1), g_quasipot.args.npoints(2));
y = reshape(g_quasipot.run.parvalues(:, 2), g_quasipot.args.npoints(1), g_quasipot.args.npoints(2));
z1 = griddata(y1(:), y2(:), pot2(:), x, y);
if nargout >0
    res.X = x;
    res.Y = y;
    res.pot = z1;
    res.attract=attract;
end
doplot(args,x,y,z1);



function doplot(args,x,y,z1,curl)
global g_grind
if ~args.silent
    z = z1;
    if ~isnan(args.maxpot)
        z(z > args.maxpot) = nan;
    end
    hfig=i_makefig('quasipot');
    i_plotdefaults;
    hold on;
    surfl(x, y, z);
    shading interp
    colormap bone;
    xlabel(g_grind.xaxis.var)
    ylabel(g_grind.yaxis.var)
    alg=args.algorithm;
    alg(1)=upper(alg(1));
    zlabel(sprintf('Potential (%s)',alg))
    set(get(hfig,'CurrentAxes'),'view' ,[-44 68]);
    % -- Plot contour lines --
    hold on;
    contour3(x,y,z, 30);
    hold off;
    axis tight
    if ~isnan(args.maxpot)
        zlim([min(z1(:)) min(max(z1(:)), args.maxpot)]);
    else
        zlim([min(z1(:)) max(z1(:))])
    end
    if nargin>4
        hfig=i_makefig('quasipot',1);
        if ~isoctave&&verLessThan('matlab','8.6')
            set(hfig,'renderer','ZBuffer'); %in fact this is solving a window bug
        end
        pcolor(x,y,curl);
        h=colorbar;
        h1=ylabel(h,'Curl factor');
        set(h1, 'FontSize', g_grind.pen.fontsize);
        i_plotdefaults;
        shading flat
        xlabel(g_grind.xaxis.var)
        ylabel(g_grind.yaxis.var)
        zlabel(sprintf('Curl factor'))
    end
    
end

% for i=1:size(g_quasipot.run.Y,3)
%     potz=pot(:,1,i);
%     if ~isnan(args.maxpot)
%         potz(potz>args.maxpot)=nan;
%     end

%     plot3(y1(:,1,i),y2(:,1,i),potz,'r-');
% end

function changed = settingschanged(args)
global g_quasipot;
if isfield(g_quasipot, 'args')
    [sameargs, fields] = struccmp(args, g_quasipot.args);
    if ~isempty(fields)
        ndx=strcmp(fields,'silent');
        fields=fields(~ndx);
        ndx=strcmp(fields,'maxpot');
        fields=fields(~ndx);
        if isempty(fields)
            sameargs=true;
        end
        
    end
    
else
    sameargs = false;
    %    fields = [];
end

changed=~isfield(g_quasipot,'settings')||~struccmp(par('-v',0),g_quasipot.settings)||~sameargs;

function doupdate(silent)
%update the grid, result in
global g_grind g_quasipot;
stabil;
nm = nmaxima;
if nm.cyclic
    error('grind:quasipot', 'Quasi potential is not meaningful for cyclic models');
end

nruns = size(g_quasipot.run.parvalues, 1);

iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
N0 = [i_initvar; 0];
quasipotfun = i_getodehandle('quasipot');
if ~silent
    i_waitbar(0, nruns+size(g_quasipot.run_separ.parvalues, 1),'quasipot' ,'Calculating',0.5)
end

for i = 1:nruns
    %*** Init conds for given (x,y) ***
    %Initialize coords.
    if ~silent
        i_waitbar(1);
    end
    
    N0(iX.no) = g_quasipot.run.parvalues(i, 1);
    N0(iY.no) = g_quasipot.run.parvalues(i, 2);
    
    %Initialize global arrays (time t = 0 counts as "time step #1")
    [~,yy] = feval(str2func(g_grind.solver.name),quasipotfun, g_quasipot.run.t, N0 , g_grind.solver.opt);
    g_quasipot.run.Y(:, :, i) = yy(:, 1:end - 1);
    g_quasipot.run.pot(:, :, i) = yy(:, end);
end

for i = 1:size(g_quasipot.run_separ.parvalues, 1)
    %*** Init conds for given (x,y) ***
    %Initialize coords.
    if ~silent
        i_waitbar(1);
    end
    
    N0(iX.no) = g_quasipot.run_separ.parvalues(i, 1);
    N0(iY.no) = g_quasipot.run_separ.parvalues(i, 2);
    
    %Initialize global arrays (time t = 0 counts as "time step #1")
    [~,yy] = feval(str2func(g_grind.solver.name),quasipotfun, g_quasipot.run_separ.t, N0 , g_grind.solver.opt);
    g_quasipot.run_sep.Y(:, :, i) = yy(:, 1:end - 1);
    g_quasipot.run_separ.pot(:, :, i) = yy(:, end);
end

g_quasipot.settings = par(struct('opts', {'-v'},'statevars',0));
if ~silent
    i_waitbar([]);
end

function err=optimizeatt(pattract,g_quasipot,attract)
pattract=[0;pattract];
pot = g_quasipot.run.pot;
for i = 1:size(pot, 3)
    pot(:, 1, i) = pot(:, 1, i) - pot(end, 1, i);
end

for i = 1:length(pattract)
    pot(:, :, attract(:) == i)=pot(:, :, attract(:) == i) + pattract(i);
end

differ=[pot(1,1,g_quasipot.seppairs(:,1)),pot(1,1,g_quasipot.seppairs(:,2))];
differ=permute(differ,[3,2,1]);
err=sum((differ(:,1)-differ(:,2)).^2);





