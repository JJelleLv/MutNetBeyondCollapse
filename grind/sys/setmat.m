%SETMAT   Set the values of a matrix (disturbance/gradient)
%   Use this command to set easily the values of a matrix.
%  
%   Usage:
%   SETMAT - opens a window to set the values
%   SETMAT A VAL - sets all values of the matrix A to VAL
%   SETMAT A [L1,L2] SIZE VAL - sets in matrix A at [L1,L2] a square of SIZE cells to VAL
%   SETMAT A -centre SIZE VAL - disturbance VAL,SIZE at the centre
%   SETMAT A -random SIZE VAL - disturbance VAL,SIZE at a random position
%   SETMAT A [L1,L2] [S1,S2] VAL - sets a rectangle of S1xS2 to VAL
%   SETMAT A LOC SIZE [V1,V2] - sets outside the disturbance to V1, inside to V2
%   SETMAT A -xgradient [V1,V2] - sets an x gradient from V1 to V2
%   SETMAT A -ygradient [V1,V2] - sets an y gradient from V1 to V2
%   SETMAT A -uniform [V1,V2] - draws at random values between V1 and V2 (default [0 1])
%   SETMAT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'location' [-centre | -random or all(integer>0) and length(integer)==2] - location of perturbation
%     'num' [integer>=0] - number of (random) perturbations
%     'size' [all(integer>=0) and length(integer)<=2] - size of the perturbation
%     'value' [number] - value for the perturbation
%     'var' [string or number] - name the matrix
%   SETMAT('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - change a centre rectangle
%     '-r' - random location
%     '-u' - uniform draws at random values between V1 and V2 (default [0 1])
%     '-x' - gradient in x direction
%     '-xh' - half of the field in the x direction
%     '-y' - gradient in the y direction
%     '-yh' - half of the field in the y direction
%  
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setmat')">commands setmat</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [A, locs] = setmat(varargin)
%(var, loc, siz, value, num)
global g_grind;
fieldnams={'var', 's#n', 'name the matrix';...
   'location', 'e[-c+entre|-r+andom]#all(i>0)&length(i)==2', 'location of perturbation';...
   'size', 'all(i>=0)&length(i)<=2', 'size of the perturbation';...
   'value', 'n', 'value for the perturbation';...
   'num', 'i>=0', 'number of (random) perturbations'}';
args=i_parseargs(fieldnams,['if((nargs==2)||hasoption(''-x,-y,-u'')),deffields=''var,value'';',...
    'elseif(hasoption(''-c'')||hasoption(''-r'')),deffields=''var,size,value,num'';else,deffields=''var,location,size,value'';end'],...
    '-xh,-yh,-c,-r,-u,-x,-y',varargin);
if ~isfield(args,'var')
    if ~isempty(g_grind)
        var = g_grind.statevars.vectnames{1};
    else
        var = '';
    end
    if ~isempty(g_grind) && isfield(g_grind, 'setmat')
        answer = g_grind.setmat;
    else
        answer = {var, '-centre', '[3 3]', '', ''};
    end
    prompt={'Name of matrix variable','Location of disturbance/option ([x y] -centre -random)',...
        'Size of disturbance (number of cells) or [nx ny]','New value of disturbance or [outside inside]', 'Divide area in n parts (only with -random)'};
    answer = inputdlg(prompt, 'set matrix', 1, answer);
    if isempty(answer)
        disp('Cancelled');
        return;
    end
    args.var = answer{1};
    args.location = answer{2};
    if (args.location(1) ~= '-')
        args.location = i_checkstr(args.location);
    end
    args.size = i_checkstr(answer{3});
    args.value = answer{4};
    args.numparts = i_checkstr(answer{5});
    if ~isempty(g_grind)
        g_grind.setmat = answer;
    end
end
if ~isfield(args,'location')&&numel(args.opts)==1
    args.location=args.opts{1};
end
if isfield(args,'size')
    args.size=round(args.size);
end
if ischar(args.var)
    E = str2num(args.var);  %#ok<ST2NM>
    if isempty(E)
        E = evalin('base', args.var);
    end
    n = size(E);
else
    E = args.var;
    args.var = '';
    n = size(E);
end
if isfield(args,'location')&&ischar(args.location)
    if strncmpi(args.location, '-c', 2)
        args.location = round([n(1), n(2)] ./ 2);
        E = setdisturbance(E, n, args.location, args.size, args.value);
    elseif strncmpi(args.location, '-xh', 3)
        siz=size(E);
        E(:,1:round(siz(1)/2))=i_checkstr(args.value);
    elseif strncmpi(args.location, '-yh', 3)
        args.size=size(E);
        E(1:round(args.size(1)/2),:)=i_checkstr(args.value);
    elseif strncmpi(args.location, '-r', 2)
        if ~isfield(args,'num')
            args.num=1;
        end
        sizes = intdivide(args.size, args.num);
        found=0;
        % set args.num non-overlapping not bordering random locations
        Etmp=E;
        Eopt=[];locsopt=[];
        maxsize=0;
        j=0;
        while ~found
            nuls=zeros(size(E));
            locs=[0 0];
            E=Etmp;
            for k = 1:length(sizes)
                args.location = floor([rand(1) * n(1), rand(1) * n(2)]) + 1;
                locs(k,:)=args.location;
                E = setdisturbance(E, n, args.location, sizes(k), args.value);
                nuls = setdisturbance(nuls, n, args.location, sizes(k), '1');
            end
            s=powerlaw(nuls);
            ssum=sum(sum(nuls));
            if ssum>maxsize
                Eopt=E;
                locsopt=locs;
                maxsize=ssum;
            end
            j=j+1;
            found=((ssum==args.size) && (max(s)<=max(sizes))) || (j>200);
        end
        if j>200
            warning('GRIND:setmat:distoverlap','Could not create non-overlapping "disturbances"');
            E=Eopt;
            locs=locsopt;
        end
    elseif strncmpi(args.location, '-u', 2)||strncmpi(args.location, '-n', 2)
        if isempty(args.size)
            args.size = i_checkstr(args.value);
        end
        if length(args.size) == 1
            if nargin == 3
                args.size = [0 args.size];
            else
                args.size = [0 1];
            end
        end
        if isnan(args.size(1)) %mean or minimum is NaN  - > add noise to current args.value
            m = E;
        else
            m = args.size(1);
        end
        if strncmpi(args.location, '-u', 2)
            E = rand(size(E)) * args.size(2) + m;
        else
            E = randn(size(E))  * args.size(2) + m;
            E(E < 0) = 0;
        end
    elseif strncmpi(args.location, '-x', 2)||strncmpi(args.location, '-y', 2)
        if ~isfield(args,'value')
            args.value = [0 1];
        end
        if length(args.value)==1
            args.value=[0 args.value];
        end
        [X,Y] = meshgrid(args.value(1):(args.value(2) - args.value(1)) / (n(1) - 1):args.value(2), args.value(1):(args.value(2) - args.value(1)) / (n(2) - 1):args.value(2));
        if strncmpi(args.location, '-x', 2)
            E = X;%assignin('base', args.var, X); plotmat(X,args.var);
        else
            E = Y;%assignin('base', args.var, Y); plotmat(Y,args.var);
        end
    elseif strncmpi(args.location, '-s', 2)
        E = ones(size(E)) * evalin('base', args.value);
        %      assignin('base', args.var, E); plotmat(E,args.var);return;
    end
elseif ~isfield(args,'location')&&~isfield(args,'size')
    E(:)=args.value;
else
    E = setdisturbance(E, n, args.location, args.size, args.value);
end

if ~isempty(args.var)&&ischar(args.var)
    try
       assignin('base', args.var, E);
    catch err
        if ~strcmp(err.identifier,'MATLAB:assigninInvalidVariable')
           rethrow(err)
        end
    end
end
if nargout > 0
    A = E;
else
    plotmat(E, args.var);
end

function E = setdisturbance(E, n, loc, siz, value)
if siz > 0
    [xx, yy] = getdistloc(siz);
    xx = xx + loc(1);
    xx(xx <= 0)=xx(xx <= 0) + n(1);
    xx(xx > n(1)) = xx(xx > n(1)) - n(1);
    yy = yy + loc(2);
    yy(yy <= 0)=yy(yy <= 0) + n(2);
    yy(yy > n(2)) = yy(yy > n(2)) - n(2);
    if ischar(value)
        aval = evalin('base', value);
    else
        aval=value;
    end
    if length(aval) == 2
        E = zeros(size(E)) + aval(1);
        aval = aval(2);
    end
    E(sub2ind(size(E), xx, yy)) = aval;
end

function plotmat(E, var)
global g_grind;
hfig=i_makefig('setmat');
if min(size(E)) == 1
    plot(1:length(E), E);
    if ~isempty(g_grind)
        htitle=title(sprintf('initial conditions %s', var));
        i_plotdefaults(hfig);
    else
        htitle=title(sprintf('matrix %s', var));
    end
    set(htitle,'fontweight','normal');
    xlabel('Row nr');
else
    pcolor(E);
    if ~isempty(g_grind)
        colormap('i_grindcolormap');
        i_plotdefaults(hfig);
        htitle=title(sprintf('initial conditions %s', var));
    else
        htitle=title(sprintf('matrix %s', var));
    end
    set(htitle,'fontweight','normal');
    set(gcf, 'renderer', 'zbuffer');
    daspect([1 1 1]);
    if ~isoctave&&verLessThan('matlab','8.4.0')
        set(gca, 'drawmode','fast','ydir','reverse');
    else
        set(gca,'SortMethod','depth','ydir','reverse');
    end
    xlabel('Column');
    ylabel('Row');
    shading flat;
    colorbar;
    i_plotdefaults(hfig);
end

%get x's and y's of a disturbance of n cells (round zero)
%keep square as much as possible
function [Xs, Ys] = getdistloc(n)
if length(n) == 2
    if sum(abs(n - round(n))) > 0.05
        [Xs, Ys] = getdistloc(prod(n));
    else
        [Xs, Ys] = meshgrid(-floor(n(1) / 2):n(1)-floor(n(1) / 2)-1,...
            -floor(n(2) / 2):n(2)-floor(n(2) / 2)-1);
        Xs = transpose(Xs(:));
        Ys = transpose(Ys(:));
    end
else
    n = round(n);
    root = floor((floor(sqrt(n)) + 1) / 2) * 2 - 1;
    siz = (root - 1) / 2;
    [Xs, Ys] = meshgrid(-siz:siz, -siz:siz);
    borderX = -[ones(1, siz * 2 + 1) * (-siz-1) (-siz-1:1:siz + 1)  ones(1, siz * 2 + 1) * (siz + 1) (siz + 1:-1:-siz) (-siz-1)];
    borderY = [siz:-1:-siz-1 ones(1, siz * 2 + 1) * (-siz-1) (-siz-1:1:siz + 1) ones(1, siz * 2 + 2) * (siz + 1)];
    Xs = [transpose(Xs(:)) borderX(1:n - length(Xs(:)))];
    Ys = [transpose(Ys(:)) borderY(1:n - length(Ys(:)))];
end

function A = intdivide(n, parts)
h = 0;
part = floor(n / parts);
addh = n / parts - part;
A = zeros(1,parts);
for i = 1:parts
    h = h + addh;
    if h > 0.99
        h1 = floor(h + 0.001);
        h = h - h1;
    else
        h1 = 0;
    end
    A(i) = part + h1;
end

%function loc=stata(i,num)
%if num<4
%   loc=rand(1,2);
%elseif num<8
%   loc(2)=rand(1);
%   if i<num/2
%      loc(1)=rand(1)*0.5;
%   else
%      loc(2)=rand(1)*0.5+0.5;
%   end
%elseif num<12
%   loc=rand(1,2);
%   if i>num/2
%      loc(1)=loc(1)+0.5;
%   end
%   if (i>num/4)&(i<3*num/4)
%      loc(2)=loc(2)+0.5;
%   end
%end



