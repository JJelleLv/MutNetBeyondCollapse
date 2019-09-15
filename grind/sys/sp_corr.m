function [I1] = sp_corr(A, lag, binsize)

% SP_CORR Function to calculate the spatial correlation of a matrix for
% increasing distance class of spatial lag d using the Moran's coefficient I given by:
%
% I(d)=[Sum Sum wij(d)(xi-xmean)(xj-xmean)/W(d)]/[Sum(xi-xmean)^2/n]
%
% Usage:
% [I] = sp_corr(A, lag, binsize);
% A        refers to the matrix for which the correlogram has to be estimated
% lag      is the distance up to which the Moran's coefficient will be calculated
%          by default the maximum width of the matrix)
% binsize  The size of the bins in the plot (default 1)
% I        the resulting spatial correlations are in this matrix (first column: lag; second column: correlation;
%          third column the number of points used. 
%
% References: Legendre P. & Fortin M.J.(1989). Spatial pattern and ecologival analysis. Vegetation 80: 107-138.

%function [lag1, lagthres, Imax, lagmax] = sp_corr(A, flag, lag, thres)

if nargin < 2
   lag = length(A);
end
if nargin <3
   binsize = 1; % the size of the bins in the distance
end
n = length(A(:));
s = size(A);
%optimize the memory usage
if lag > length(A)
   lengthM = n * (n - 1);
elseif lag <= 1
   lengthM = 4 * n;
elseif lag < max(s) / 2
   lengthM = n * (2 * round(lag * lag) - 1);
else
   lengthM = n * (round(lag * lag) - 1);
end
maxlen = 50^4;
if lengthM > maxlen
   M = zeros(maxlen, 1);
else
   M = zeros(lengthM, 1);
end
D=M;
%M=zeros(0,3);
xm = mean(A(:));
varA = var(A(:), 1);  	% estimate the variance of the matrix
rows = repmat(transpose(1:s(1)), 1, s(2));
cols = repmat(1:s(2), s(1), 1);
k = 1;
for i = 1:n - 1
   dists = transpose(sqrt((rows(i + 1:end) - rows(i)).^2 + (cols(i + 1:end) - cols(i)).^2));
   ndx=find(dists <= lag);
   L = length(ndx);
   xs = A(i) + zeros(L, 1);
   ys = transpose(A(i + 1:n));
   M(k:k + L - 1) = (xs - xm) .* (ys(ndx) - xm);
   D(k:k+L-1)= dists(ndx);
   k = k + L;
end
M = M(1:k - 1);
D =D(1:k-1);
[D,ndx]=sort(D);
M=M(ndx);
%M(:,1:2)=M(:,1:2)-xm; 					% calculate xi-xmean and xj-xmean
%MM=[M(:,1).*M(:,2), M(:,3)];        % calculate product xi-xmean and xj-xmean


m = length(M);
nbins = floor(D(end) / binsize);
I = zeros(nbins, 3);
k = 1;
previ = 1;
bin1 = D(1) + binsize;
for i =1:m
   if D(i) > bin1	% estimate breaks where clusters of same distances exist
      I(k, 1) = bin1 - binsize / 2;
      I(k, 2) = mean(M(previ:i));
      I(k, 3) = i - previ;
      previ = i + 1;
      bin1 = bin1 + binsize;
      k = k + 1;
   end
end
I = I(1:k-1, :);
I(k, 1) = bin1 - binsize / 2;
I(k, 2) = mean(M(previ:end));
I(k, 3) = m - previ;

%b=[0;bb';m];
%nn=length(b);  			% number of different distance lags
%I=zeros(nn-1,2);
%for i=1:nn-1
%   I(i,1)=mean(MM(b(i)+1:b(i+1)));
%   I(i,2)=MM(b(i+1),2);
%end
I(:, 2) = I(:, 2) / varA;

%lag1 = I(1, 1);
%Imax = max(I(:, 1));
lagmax = I(end, 1);
%a = find(I(:, 1) < thres);
%lagthres = I(a(1), );
%b=find(I(:,2)<0.01 & I(:,2)>-0.01);
%if ~isempty(b) %  lagzero=I(b(1),2);
%end


if nargout==0
   figure(1);
   plot(I(:, 1), I(:, 3), '.-');
   set(gca, 'xlim', [0 lagmax]);
   xlabel('distance');
   ylabel('n');
   title('Number of points used per bin');
   figure(2);
   plot(I(:, 1), I(:, 2), '.-');
   set(gca, 'xlim', [0 lagmax]);
   xlabel('distance');
   ylabel('spatial correlation');
   htitle=title('correlogram');
   set(htitle,'fontweight','normal');
else
   I1=I;
end



