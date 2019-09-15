function A = randlogn(n1, n2, mean, sd)
%mean is the expected mean
%sd is expected standard dev
%randlogn draws from lognormal distribution
if sd==0 || mean==0
   A = ones(n1, n2) * mean;
else
   mu = log(mean^2 / sqrt(sd^2 + mean^2));
   sigma2 = log((sd / mean)^2 + 1);
   A = exp(mu + sigma2 * randn(n1, n2));
end
