function r = stable_rnd(alpha,beta,gam,delta,sizeOut)
%STABLE_RND Random arrays from the stable distribution.
%
% stable_rnd(alpha,beta,gam,delta,sizeOut)
%    1,0,1,0 = Cauchy distribution
%    0.5 -1,1,0 = Levy distribution
%    2,0,1,0 = Normal distribution
%  alpha[0-2  = kurtosis parameter (very sensitive)
%  beta = skewness parameter
%  gam = scaling parameter (relates to sd)
%  delta = mean
% References:
%       [1] A. Weron and R. Weron (1995), "Computer Simulation of Levy
%       alpha-Stable Variables and Processes", Lecture Notes in Physics, 
%       Springer-Verlag.
%       [2] J.P. Nolan (2015), "Stable Distributions - Models for Heavy
%       Tailed Data", Birkhauser, Boston. In progress, Chapter 1 online at
%       academic2.american.edu/~jpnolan

if nargin<2
    beta=0;
end
if nargin<3
    gam=1;
end
if nargin<4
    delta=0;
end
if nargin<5
    sizeOut=1;
end
numelOut=prod(sizeOut);

if numelOut ~= 1
    if isscalar(alpha), alpha = repmat(alpha,sizeOut); end
    if isscalar(beta), beta = repmat(beta,sizeOut); end
    if isscalar(gam), gam = repmat(gam,sizeOut); end
    if isscalar(delta), delta = repmat(delta,sizeOut); end
end

% Return NaN for illegal parameter values.
alpha(alpha<=0 | alpha>2) = NaN;
beta(beta<-1 | beta>1) = NaN;
gam(gam <= 0) = NaN;

% Preallocation
r = zeros(sizeOut);

% Generate the STANDARD stable distribution, gam = 1 and delta = 0.
%   Use corresponding method for special cases (Normal and Cauchy), 
%   otherwise use general method.
Gau = alpha == 2;
if any(Gau(:)) % Gaussian distribution
    r(Gau) = sqrt(2) * randn([sum(Gau(:)),1]);
end

Cau = (alpha == 1) & (beta == 0);
if any(Cau(:))
    % Cauchy(0,1) is t-distribution with d.f. = 1
    r(Cau) = trnd(1,[sum(Cau(:)),1]);
end

Levy = (alpha == 0.5) & (abs(beta)==1);
if any(Levy(:)) % Levy distribution 
    % Levy(0,1) = S(0.5,1,1,1;0), convert it to S(0.5,1,1,0;0) by subtracting -1
    r(Levy) = 1./ randn([sum(Levy(:)),1]).^2 - 1;
    r(Levy) = beta(Levy).*r(Levy);
end

other = (~Cau) & (~Gau) & (~Levy);
if any(other(:))
    U = pi.*(rand(sizeOut)-0.5); % Uniform distribution on (-pi/2,pi/2)
    W = -log(rand(sizeOut)); % Exponential distribution with mean 1
   
    alphaE1 = (alpha == 1) & other;
    alphaN1 = (alpha ~= 1) & other;
    if any(alphaE1(:))
        r(alphaE1) = (2/pi) .* ( (pi/2 + beta(alphaE1) .*U(alphaE1)).*tan(U(alphaE1)) - ...
             beta(alphaE1).*log( (pi/2.*W(alphaE1).*cos(U(alphaE1)))./(pi/2 + beta(alphaE1).*U(alphaE1)) ));  
    end
    if any(alphaN1(:))
        Btan = beta(alphaN1).*tan(pi*alpha(alphaN1)/2);
        aUB = alpha(alphaN1).*(U(alphaN1) + atan(Btan)./alpha(alphaN1));
        r(alphaN1) = (1 + Btan.^2).^(1./(2*alpha(alphaN1))).*...
                    sin( aUB ) ./ (cos(U(alphaN1)).^(1./alpha(alphaN1))) .* ...
                    ( cos(U(alphaN1)-aUB) ./ W(alphaN1) ).^((1-alpha(alphaN1))./alpha(alphaN1));  
        % shift when alpha~=1
        r(alphaN1) = r(alphaN1) - Btan; 
    end
end

% Scale and shift from the STANDARD stable distribution
r = gam.*r + delta;

end