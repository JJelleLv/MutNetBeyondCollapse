%CATFUN = analyse matrix by catagories, the matrix is assumed to be sorted on the
% category
%
% [MAT]=CATFUN('AFUN',MAT,CATCOL) analyse matrix MAT, category is in column CATCOL
% (or empty if there is only one category) default CATCOL=[].
% [MAT,NDXS]=CATFUN('AFUN',MAT,ACATS) analyse matrix MAT, the category column is in a
% separate array ACATS. NDXS gives the indexes of the results.
%
%
% AFUN = one of these functions:
%   unchanged = return 
%   rangexx = xx% of the data are within this range
%   sum = sum of the catagory
%   count = coefficient of variation
%   cv = coefficient of variation
%   var = variance
%   percxx = xx% percentile
%   max = maximum
%   min = minimum
%   mean = avarage
%   median = median
%   sd = standard deviation
%   minima = only local minima
%   maxima = only local maxima
%   minima+maxima = both local minima and maxima

function [A, ndxs] = catfun(afun, mat, acat, par1)
if nargin < 2
   error('GRIND:catfun:ArgEror','Too few arguments');
end
if nargin < 3
   acat = [];
else
   acat = i_checkstr(acat);
end
if nargin < 4
   par1 = [];
else
   par1 = i_checkstr(par1);
end
mat = i_checkstr(mat);
afun = lower(afun);
if size(acat, 1) == size(mat, 1)
   cats = acat;
   acat = [];
   firstcol = 1;
elseif isempty(acat)
   cats =  ones(size(mat, 1), 1);
   firstcol = 1;
elseif size(mat,2)==1
   if length(acat)==1
       cats=zeros(size(mat))+acat;
   else
      cats=acat;
   end
   firstcol=1;
else
   if ~isempty(mat)
      cats = mat(:, acat);
   else
      cats =mat;
   end
   if acat == 1
      firstcol = 2;
   else
      firstcol = 1;
   end
end
if strcmp(afun, 'unchanged')
   A = mat;
   ndxs = [];
   return;
end

if strncmp(afun, 'perc', 4) && (length(afun)>4)
   par1 = str2double(afun(5:end)) / 100;
   if isempty(par1)
      par1 = .5;
   end
   afun = 'perc';
end
if strncmp(afun, 'range', 5) && (length(afun)>5)
   par1 = str2double(afun(6:end)) / 100;
   afun = 'range';
end
if isempty(cats)
   icats=cats;
else
   icats=[0; find(diff(cats(:,1)) ~= 0); size(cats,1)];
end
noutputs=size(mat,2);
ncat=length(icats)-1;
A = zeros(ncat,noutputs);
ndxs  = zeros(ncat,1);
k=1;
for i = 2:length(icats)
   [a,xxs] = thefun(afun, mat(icats(i - 1) + 1:icats(i), :), par1, firstcol);
   if ~isempty(a)
      a(:, acat) = cats(icats(i));
%     c = zeros(size(a, 1), size(cats,2)) + cats(icats(i),:);
      A(k:k+size(a,1)-1,:) = a;
      %A= [A; a];
      %ndxs = [ndxs; xxs+icats(i-1)];
      ndxs(k:k+size(a,1)-1) = xxs+icats(i-1);
      k=k+size(a,1);
   end
end
function [a,xxs] = thefun(afun, m, par1,firstcol)
xxs=1;
switch afun
 case 'range'
   if isempty(par1) || isnan(par1)
      a = [min(m, [], 1); max(m, [], 1)];
   else
      a = [i_makepercentiles(m', (1-par1)/2)' ;i_makepercentiles(m', par1+(1-par1)/2)'];
   end
 case 'minima'
   if isempty(par1)
      par1 = firstcol;
   end
   xxs = find(diff(sign((diff(m(:, par1))))) > 0) + 1;
   if isempty(xxs)
      xxs = 1;
   end
   a = m(xxs, :);
 case 'maxima'
   if isempty(par1)
      par1 = firstcol;
   end
   xxs = find(diff(sign(diff(m(:, par1)))) < 0) + 1;
   if isempty(xxs)
      xxs = 1;
   end
   a = m(xxs, :);
 case {'minima+maxima' 'maxima+minima'}
   if isempty(par1)
      par1 = firstcol;
   end
   xxs = [find(diff(sign(diff(m(:, par1)))) > 0) + 1; find(diff(sign(diff(m(:, par1)))) < 0) + 1];
   if isempty(xxs)
      xxs = 1;
   end
   a = m(xxs, :);
case 'count'
   m1= ones(size(m));
   a=sum(m1,1);
 case 'sum'
   a = sum(m, 1);
 case 'cv'
   a = mean(m, 1) ./ std(m, 0, 1);
 case {'variance','var'}
   a = std(m, 0, 1).^2;
 case {'autoregr'}
   a = autoregression(m, 1);
 case 'perc'
   a = i_makepercentiles(m', par1)';
 case 'max'
   a = max(m, [], 1);
 case 'min'
   a = min(m, [], 1);
 case 'mean'
   a = mean(m, 1);
 case 'median'
   a = median(m, 1);
 case 'sd'
   a = std(m, 0, 1);
 otherwise
   a = [];
end
