%autoregression ols (without intercept and with subtraction the mean)
%
function rc = autoregression(Ys, alag)
if nargin<2
    alag=1;
end

rc = zeros(length(alag), size(Ys, 2));
Ys=Ys-mean(Ys(2:end)); %subtracting the mean 
%                       (important if you do regression without intercept!)
for j = 1:length(alag)
   for i = 1:size(Ys, 2)
      res = linregres(Ys(1:end-alag, i), Ys(1+alag:end, i));
      rc(j, i) = res(1);
   end

end

function b=linregres(x,y) %without intercept
%b=(sum(x.*y)-1/length(x)*sum(x)*sum(y))/(sum(x.^2)-1/length(x)*sum(x).^2);
b=sum(x.*y)./(sum(x.^2)); %without intercept
%intercept
%a=mean(y)-b*mean(x);
