% OUTF   Plot descriptive statistics of vector or matrix state variable(s)
%  Various function on vector or matrix state variables.
%
%  Usage:
%  OUTF('FUN','VARNAME') - function FUN of vector VARNAME.
%  OUTF('FUN','VARNAME(20:30,20:30)') - function FUN of the block 20:30,20:30 of the matrix VARNAME.
%  OUTF('FUN',FROMNO,TONO) - function FUN of state variable FROMNO to TONO.
%  OUTF('SHANNON',FROMNO,TONO, THRESHOLD) - set the threshold of absence for Shannon index to THRESHOLD. 
%  The default threshold is 0.01.
%
%  In <a href="matlab:help out">out</a> and <a href="matlab:help ax">ax</a> you can shortcut OUTF as _FUN(VARNAME) or 
%  _FUN(VARNAME(20:30:20:30)) or _FUN(FROMNO,TONO) 
%  Underscore is used to identify that it is not the original MATLAB function (not to be used in 
%  combination with other functions)
%      
%          
%  Functions:
%    cover = number of cells above threshold p (matrix variables)
%    cover< = number of cells below threshold p (matrix variables)
%    cv = plot the coefficient of variation of vector state variable
%    domeigen = dominant eigenvalue
%    gini = Gini coefficient of vector state variable(s)
%    imageigen = imaginary part of the eigenvalues
%    max = maximum of vector state variable(s)
%    maxa = only plot maxima of a state variable (else NaN)
%    mean = average of vector state variables(s)
%    median = median of vector state variable(s)
%    min = minimum of vector state variable(s)
%    mina = only plot minima of a state variable (else NaN)
%    minmaxa = only plot minima+maxima of a state variable (else NaN)
%    mov_autoregression = moving autocorrelation after detrending with kernel smoothing.
%    mov_mean = moving average
%    obserror,X = add white noise (Normally distributed) to the result, X=standard deviation of the noise
%    perc5 = 5% percentile of vector state variable(s)
%    perc95 = 95% percentile of vector state variable(s)
%    percX = X% percentile of vector state variable(s)
%    realeigen = real part of the eigenvalues
%    shannon = Shannon Weaver index of state variables (sometimes called Shannon Wiener index):
%		H` = -sum(p(i) * ln(p(i))) In which: p(i)=X(i)/sum(X)
%    smooth = Kernel Gaussian Smoothing
%    std = standard deviation of vector state variables
%    sum = sum of vector state variables
%    tmaxdiff = time with the maximum increase or decrease of a state variable
%    totmean, totstd, totvar = only plot the absolute mean/std/var of the state variable in time
%    totmin, totmax = only plot the absolute minimum/maximum of the state variable in time
%    var = variance of vector state variables          
%
%
%  See also out, ax
%
%   Reference page in Help browser:
%      <a href="matlab:commands('outf')">commands outf</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function res = outf(varargin)

global g_Y g_t g_permanent g_func g_grind;
%select variables
fun = lower(varargin{1});
if strcmp(fun, 'changeshortcut')
   res = varargin{2};
   if ~isempty(res) && (res(1) == '_')
      funcs={'_cover','_cover<','_cv','_cumsum','_domeigen','_imageigen','_max','_maxa',...
         '_mean','_median','_min','_mina','_minmaxa','_mov','_obserror','_gini',...
         '_perc','_realeigen','_shannon','_smooth','_std','_sum','_tmaxdiff',....
         '_totmax','_totmean','_totmedian','_totmin','_totstd','_totvar',...
         '_var','_variance'};
      if strncmpi(res, '_perc', 5)
         i = 6;
         while strcontains('0123456789.', res(i))
            i = i + 1;
         end
         funcs{1} = res(1:i - 1);
      elseif strncmpi(res, '_mov', 4)
         i = strfind(res, '(');
         if isempty(i)
            funcs{2} = res;
         else
            funcs{2} = res(1:i - 1);
         end
      end
      for i = 1:length(funcs)
         if strncmpi(res, funcs{i}, length(funcs{i}))
            if length(res) == length(funcs{i})
               res=sprintf('outf(''%s'')',res(2:length(res)));
               return;
            elseif strcontains('(0123456789', res(length(funcs{i}) + 1))
               res=['outf(''' res(2:length(funcs{i})) ''',' res(length(funcs{i})+2:end)];
               f = strfind(res, ')');
               if ~isempty(f)
                  f = f(end);
               end
               f = sort([f strfind(res,',')]);
               f2=[strfind(res,''')') strfind(res,''',')]+1;
               f3=[strfind(res,')') strfind(res,',''')];
               f4 = strfind(res, ':');
               if ~isempty(f4) && (length(f4) == 2)
                  f4(2) = f4(2) + 1;
                  while res(f4(2)) ~= ')'
                     f4(2) = f4(2) + 1;
                  end
               end
               for i1 = length(f):-1:1
                  if isempty(f4) || length(f4) ~= 2 || f(i1) < f4(1) || f(i1) > f4(2)
                     if isempty(f3) || ~any(f3 == f(i1))  %isempty(find(f3 == f(i1)))
                        res=[res(1:f(i1)) '''' res(f(i1)+1:end)];
                     end
                     if isempty(f2) || ~any(f2 == f(i1))
                        res=[res(1:f(i1)-1) '''' res(f(i1):end)];
                     end
                  end
               end
            end
         end
      end
   end
   return;
end
if nargin > 1
   avar = varargin{2};
else
   avar = [];
end
if nargin >= 3
   sto = i_checkstr(varargin{3});
else
   sto = [];
end

if nargin >= 4
   par1 = i_checkstr(varargin{4});
else
   par1 = [];
end
if nargin >= 5
   par2 = i_checkstr(varargin{5});
else
   par2 = [];
end
if nargin == 1
   Ys = g_Y;
else
   if ischar(avar)
      sfrom = str2num(avar); %#ok<ST2NM> %  leave this way
   else
      sfrom = avar;
   end
   if isempty(sfrom)
      if isempty(par1)
         par1 = sto;
      end
      if isempty(par2)
         par2 = par1;
         par1 = sto;
      end
      if strcontains(avar, ':')
         f = strfind(avar, '(');
         if ~isempty(f)
            iX = i_getno(avar(1:f(1) - 1));
         end
      else
         iX = i_getno(avar);
      end
      if iX.isvar
         [sfrom, sto] = i_statevarnos(avar);
      elseif iX.isperm
         [sfrom, sto] = permanentvarnos(avar);
      elseif iX.isfun
         [sfrom, sto] =  funvarnos(avar);
      else
         Ys = outfun(avar);
         %disp('Error:outf cannot evaluate functions')
      end
   end
   if isempty(sto)
      L = sfrom;
   else
      L = sfrom:sto;
   end
   if iX.isvar
      Ys = g_Y(:, L);
   elseif iX.isperm
      Ys = g_permanent.Y(:, L);
   elseif iX.isfun
      i_update_g_func;
      Ys = g_func(:, L);
   end
end

%adapt here the function'
if strncmp(fun, 'mov_',4) && (length(fun) > 4)
   fun1 = fun(5:end);
   fun = 'moving ';
end
if strncmp(fun, 'perc', 4)
   par1 = str2double(fun(5:end)) / 100;
   if isempty(par1)
      par1 = .5;
   end
   fun = 'perc';
end
switch fun
 case 'cumsum'
   res = cumsum(Ys,1);
 case 'min'
   res = min(Ys, [], 2);
 case 'max'
   res = max(Ys, [], 2);
 case 'mean'
   res = mean(Ys, 2);
 case 'totmax'
   res = nan(size(Ys));
   res(end, :) = max(Ys, [], 1);
 case 'totmedian'
   res = nan(size(Ys));
   res(end, :) = median(Ys, 1);
 case 'totmean'
   res = nan(size(Ys));
   res(end, :) = mean(Ys, 1);
 case 'totstd'
   res = nan(size(Ys));
   res(end, :) = std(Ys, 1);
 case 'totvar'
   res = nan(size(Ys));
   res(end, :) = var(Ys, 1);
 case 'totmin'
   res = nan(size(Ys));
   res(end, :) = min(Ys, [], 1);
 case 'tmaxdiff'
   res = nan(size(Ys));
   diffY = abs(diff(Ys));
   [i, ~]=find(repmat(max(diffY), size(Ys, 1) - 1, 1) == diffY);
   res(end, :) = transpose(g_t(i));
 case 'minmaxa'
   if size(Ys, 2) > 1
      [~, n] = max(sum(Ys));
      Ys = Ys(:, n);
   end
   index = [find( diff( sign( diff([0; Ys; 0]) ) ) < 0 ); find( diff( sign( diff([0; Ys; 0]) ) ) > 0 )];
   if length(index) > 2
      index = index(1:length(index) - 1);
   end
   if (length(index) > 2) && (index(1) == 1)
      index = index(2:length(index));
   end
   res = NaN * zeros(size(Ys));
   res(index) = Ys(index);
   
 case 'maxa'
   if size(Ys, 2) > 1
      [~, n] = max(sum(Ys));
      Ys = Ys(:, n);
   end
   index = find( diff( sign( diff([0; Ys; 0]) ) ) < 0 );
   if length(index) > 2
      index = index(1:length(index) - 1);
   end
   if (length(index) > 2) && (index(1) == 1)
      index = index(2:length(index));
   end
   res = NaN * zeros(size(Ys));
   res(index) = Ys(index);
 case 'mina'
   if size(Ys, 2) > 1
      [~, n] = max(sum(Ys));
      Ys = Ys(:, n);
   end
   index = find( diff( sign( diff([0; Ys; 0]) ) ) > 0 );
   if length(index) > 2
      index = index(1:length(index) - 1);
   end
   if (length(index) > 2) && (index(1) == 1)
      index = index(2:length(index));
   end
   res = NaN * zeros(size(Ys));
   res(index) = Ys(index);
 case 'obserror'
   if isempty(par1)
      par1 = 1;
   end
   res = Ys + randn(size(Ys)) * par1;
 case 'sum'
   res = sum(Ys, 2);
 case 'median'
   res = median(Ys, 2);
 case 'cv'
   res = mean(Ys, 2) ./ std(Ys, 0, 2);
 case 'std'
   res = std(Ys, 0, 2);
 case {'variance','var'}
   res = var(Ys, 0, 2);
 case 'perc'
   res = quantile(Ys', par1)';
 case {'cover','cover>'}
   if isempty(par1)
      par1 = 0.01;
   end
   res = sum(Ys > par1, 2);
 case 'realeigen'
   res = ones(size(g_Y));
   for i = 1:size(g_Y, 1)
      [~, eigenval] = i_eigen(isempty(g_grind.syms.Jacobian),g_grind.solver.iters,transpose(g_Y(i,:)));
      res(i, :) = real(transpose(eigenval(:)));
   end
 case 'imageigen'
   res = ones(size(g_Y));
   for i = 1:size(g_Y, 1)
      [~, eigenval] = i_eigen(isempty(g_grind.syms.Jacobian),g_grind.solver.iters,transpose(g_Y(i,:)));
      res(i, :) = imag(transpose(eigenval(:)));
   end
 case 'domeigen'
   res = ones(size(g_Y, 1), 1);
   for i = 1:size(g_Y, 1)
      res(i) = domeigen(transpose(g_Y(i, :)));
   end
 case 'smooth'
   res = ksmooth(g_t, Ys, par1);
 case 'cover<'
   if isempty(par1)
      par1 = 0.01;
   end
   res = sum(Ys < par1, 2);
 case 'moving '
   if isempty(par1)
      par1 = g_t(end) / 10;
   end
   if strcmp(fun1, 'autoregression') && (~isempty(par2) || (par2~=0))
      aTrend = ksmooth(g_t, Ys, par2);
      oldfig = gcf;
      hfig = i_figure(double(oldfig) + 10);
      oldhold = ishold;
      plot(g_t, Ys, 'b-');
      hold on;
      plot(g_t, aTrend, 'r-');
      if ~oldhold
         hold off;
      end
      i_plotdefaults(hfig)
      xlabel('Time');
      ylabel(avar);
      i_figure(oldfig);
      Ys = Ys - aTrend;
   end
   %   lags = (g_t(1):1:g_t(end))';
   %  Ys = interp1(g_t, Ys, lags);
   res = moving_window(g_t, Ys, par1, str2func(fun1));
  case 'gini'
   res = gini(Ys);  
 case 'shannon'
   if isempty(par1)
      par1 = 0.01;
   end
   res = shannon(Ys, 1, par1);
 otherwise
   error('GRIND:outf:UnknownFunction','error in outf: function "%s" not supported', fun);
end

function maxeig = domeigen(N0)
global g_grind;
[~, eigenval] = i_eigen(isempty(g_grind.syms.Jacobian),g_grind.solver.iters,N0);
maxeig = max(real(eigenval));

function res = ksmooth(x, y, bandwidth, method)
if nargin < 3
   bandwidth = [];
end
if nargin < 4
   method = 'normal';
end
n = length(x);
if isempty(bandwidth)
   % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
   hx = median(abs(x - median(x))) / 0.6745 * (4 / 3 / n)^0.2;
   hy = median(abs(y - median(y))) / 0.6745 * (4 / 3 / n)^0.2;
   h = sqrt(hy * hx);
   fprintf('bandwidth smoothing set to %g\n', h);
   if h < sqrt(eps) * n
      error('GRIND:outf:LittleVar','There is not enough variation in the data. Regression is meaningless.')
   end
else
   h = bandwidth;
end
switch lower(method(1))
 case 'n' %"Normal" = Gaussian kernel function
   h = h * 0.36055512754640; %variance of 0.13 to get quartiles at + /  -  0.25 * h (comparable with ksmooth (R))
   res = ones(size(y, 1), 1);
   for k = 1:n
      xx = abs(x(k) - x) / h;
      z = exp(-xx(xx < 4).^2 / 2); %x < 4 is more efficient for long time series (negligible effect)
      res(k) = sum(z .* y(xx < 4)) / sum(z);
   end
 case 'b' %"box" = moving average
   d = h / 2; % 0.25 quartiles
   res =  cell(1, n);
   for k = 1:n
      xx = abs(x(k) - x) / h;
      z = xx < d;
      res(k) = sum(z .* y) / sum(z);
   end
end


function res = moving_window(ts, A, lag, funct,varargin)
% if nargin <4
%    funct = lag;
%    lag = A;
%    A = ts;
%    ts = (1:size(A, 1))';
% end
res = A .* NaN;
i1 = 1;
i = 1;
if ts(end) - ts(1) > lag
   tlag = ts + lag;
   while i1 < length(ts)
      while (i < length(ts)) && (ts(i) < tlag(i1))
         i = i + 1;
      end
      if ts(i) >= tlag(i1)
         windowA = A(i1:i, :);
         %write at end/start of the window?
         if strcmp(func2str(funct), 'autoregression')
            res(i, :) = autoregression(windowA,1);
         else
            res(i, :) = feval(funct, windowA,varargin); %end
         end
         %res(i1,:)= feval(funct, windowA); %start
         %res(floor((i1+i)/2),:)= feval(funct, windowA); %centre
      end
      i1 = i1 + 1;
      i = i - 1;
   end
end

function [sfrom, sto] = permanentvarnos(avar)
global g_grind;
sfrom = [];
%fcolon = strfind(avar, ':');
fbrack = strfind(avar, '(');
if ~isempty(g_grind.permanent)
   if ~isempty(fbrack)
      avar1 = avar(1:fbrack(1) - 1);
   else
      avar1 = avar;
   end
   sto = 0;
   for i = 1:length(g_grind.permanent)
      sfrom = sto + 1;
      sto = sto + g_grind.permanent{i}.dims(1) * g_grind.permanent{i}.dims(2);
      if strcmp(g_grind.permanent{i}.name, avar1)
         if ~isempty(fbrack)
            l = i_findindices(avar, g_grind.permanent{i}.dims);
            sfrom = l - 1 + sfrom;
            sto = [];
         end
         break
      end
   end
end
%GINI - get gini index
% if x is a matrix the gini of each row is taken
% if x is a vector, just on gini index is calculated
%
%  res=gini(x);
%
function res=gini(x)
if size(x,2)>1&&size(x,1)>1
    p = repmat(mean(x, 2), 1, size(x, 2));
    YY = sort(x, 2);
    cumYY = cumsum(YY, 2);
    res = (sum(cumsum(p, 2), 2) - sum(cumYY, 2)) ./ sum(cumsum(p, 2), 2);
else
    x=x(:);
    p=repmat(mean(x),1,length(x));
    YY=sort(x);
    cumYY=cumsum(YY);
    res=(sum(cumsum(p))-sum(cumYY))./sum(cumsum(p));
end

function [sfrom, sto]  = funvarnos(avar)
global g_grind;
i_update_g_func;
fcolon = strfind(avar, ':');
fbrack = strfind(avar, '(');
if ~isempty(fbrack)
   avar1 = avar(1:fbrack(1) - 1);
else
   avar1 = avar;
end
sto = 0;
for i = 1:length(g_grind.funcnames.dims)
   sfrom = sto + 1;
   sto = sto + g_grind.funcnames.dims{i}.dim1 * g_grind.funcnames.dims{i}.dim2;
   if strcmp(g_grind.funcnames.names{i}, avar1)
      if ~isempty(fcolon) %#ok<STREMP>
         l = sub2ind2d([g_grind.funcnames.dims{i}.dim1, g_grind.funcnames.dims{i}.dim2], avar);
         sfrom = l - 1 + sfrom;
         sto = [];
      end
      break
   end
end
