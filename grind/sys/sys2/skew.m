%simple skewness function
%the unbiased estimator (Wikipedia)
%s=skew(x,dim)
%
function s=skew(x,dim)
if nargin==1
    dim=1;
end

if dim==2
    s=skew(x',1)';
elseif dim==1
   xdev=x-repmat(mean(x,1),size(x,1),1);
   n=size(x,1);
   s=sqrt(n.*(n-1))./(n-2).*sqrt(n).*sum(xdev.^3,1)./(sum(xdev.^2,1).^(3/2));
else
    error('grind:skew', 'dim should be 1 or 2')
end

