%CCM - Convergent Cross Mapping of Sugihara et al., 2012
%A MATLAB implementation of convergent cross mapping to detect causation
%in connected time series. Basic function
%
%Citation:
%Sugihara, G., R. May, H. Ye, C. H. Hsieh, E. Deyle, M. Fogarty, and
%S. Munch. 2012. Detecting causality in complex ecosystems. Science 338:496-500.
%
%Usage:
%results = CCM(X, Ydata, E, L, tau, b)  calculates the CCM of X xmap Y (=the effect
%of Y on X, each column of Ydata) for embedding dimension E, length of library L,
%time lag (tau) for the lagged-vector construction, and number of
%point for the simplex averaging b (default E+1)
%
%results.N = the length of the predicted time series
%results.rho  = correlation coefficients between Y and Ypred
%results.mae = mean absolute errors
%results.rmse = root mean square errors
%
%Note that adding columns to Ydata takes little extra calculation time
%
%See also: CMM
function results = CCM(X, Ydata, E, L, tau, b)
if nargin < 5
   tau = 1;
end


if nargin < 6
   b = E + 1;
end


if L < E * tau
   error('CCM:Ltooshort','library L too short, must be larger than E*tau');
end

if L - E + 1 < b
   b = L - E + 1;
end

%calculate start and ends of the libraries
Lstart = 1:size(Ydata, 1) - L + 1;
Lend = Lstart + L - 1;

%generate time lags with length tau
lagsX = zeros(length(X), E) + NaN;
for i = 1:E
   ilag = (i - 1) * tau;
   lagsX(ilag + 1:length(X), i) = transpose(X(1:length(X) - ilag));
end


nY = size(Ydata, 2);
Ypred = zeros(size(Ydata,1), length(Lstart), nY) + NaN;
inm = zeros(1, b);
for i = (E-1) * tau + 1:size(Ydata, 1)
   Dtot = sqrt(sum(bsxfun(@minus, lagsX(i,:), lagsX).^2,2));
   for k = 1:length(Lend)
      %%%Why start at Lstart+E-1 and not with Lstart?
      %Lstart1=Lstart(k);
      Lstart1 = Lstart(k) + (E - 1) * tau;
      D = Dtot(Lstart1:Lend(k));
      Yobs1 =  Ydata(Lstart1:Lend(k), :);
      if (i >= Lstart1) && (i <= Lend(k))
         %leave-one-out cross validation?
         %D(i-Lstart(k)+1:i-Lstart(k)+E) = NaN;
         %else leave only the first out
         D(i - Lstart1 + 1) = NaN;
      end

      D1 = D;
      %find the b closest points
      for m = 1:min(b, length(D))
         [~, j] = min(D1);
         inm(m) = j;
         D1(j) = NaN;
      end

      Dmx = transpose(D(inm));
      %     for larger embedding dimensions b>50 is quicksort faster
      %     [~, inm] = sort(D);
      %     inm = inm(1:b);
      %     Dmx = D(inm)';
      %
      %weighted averaging of Y's of the closest points in embedded X
      if isnan(Dmx(1))
         u = nan(size(Dmx));
      elseif Dmx(1) > 0
         u = exp(-Dmx ./ Dmx(1));
         u(u < 1e-6) = 1e-6;
      else
         %perfect match at least once
         u = zeros(size(Dmx)) + 0.000001;
         u(Dmx == 0)=1;
      end

      Ypred(i, k,:) = u * Yobs1(inm,:) ./ sum(u); %summing the weights
      %    fprintf('%d\t%d\t%d\t%g\t%g  \n',i+k*1000,Lstart1,Lend(k),Ypred(i,k,1),Ydata(i,1))
   end

end
%Some statistics
rho =  zeros(length(Lend), nY);
mae = rho;
rmse = rho;
for k = 1:length(Lend)
   for j = 1:nY
      Ypred1 = Ypred(:, k, j);
      Yobs1 = Ydata(~isnan(Ypred1), j);
      Ypred1 = Ypred1(~isnan(Ypred1));
      Ypred1 = Ypred1(~isnan(Yobs1)); %remove NaNs from observations
      Yobs1 = Yobs1(~isnan(Yobs1));
      if ~isempty(Yobs1)
         rho(k,j) =  corr(Yobs1, Ypred1);
         mae(k, j) = mean(abs(Yobs1 - Ypred1));
         rmse(k, j) = sqrt(mean((Yobs1 - Ypred1).^2));
      else
         rho(k, j) =  NaN;
         mae(k, j)  = NaN;
         rmse(k, j) = NaN;
      end

   end

end

results.N = size(Ypred, 1);
results.rho  = median(rho);
results.mae = median(mae);
results.rmse = median(rmse);

