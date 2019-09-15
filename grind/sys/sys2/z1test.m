%Z1TEST   Gottwald-Melbourne 0-1 test for chaos
%  Simple and fast test to decide if a system is chaotic. Only one long time series 
%  is needed. The data should not be too densely sampled and not too spasely (ca 10 points
%  per peak works well).
%     
%  Citation:
%
%  Gottwald, G. A. and I. Melbourne. 2009. On the implementation of the 0-1 
%    test for chaos. SIAM Journal on Applied Dynamical Systems 8:129-145.
%    see  <a href="http://dx.doi.org/10.1137/080718851">http://dx.doi.org/10.1137/080718851</a> or
%    <a href="https://arxiv.org/pdf/0906.1418.pdf">https://arxiv.org/pdf/0906.1418.pdf</a>
%
%  Usage:
%   Z1TEST(X) is the result of the 0-1 test applied to the vector X.
%      Result is near to 0 for non-chaotic data and near 1 for chaotic data. By
%      default the current model run is used as time series. If X is a matrix
%      the column with the maximum average is used.
%  RES=Z1TEST - save the results of the analysis to the structure RES:
%    Fields:
%        N: number of used points
%       k: median expansion factor (1= chaotic 0 = non-chaotic
%      ks: k-values per iteration
%      cs: c-values per iteration (multiplication factor for sine and
%      cosine)
%       p, q: all p and q values for the transfomed dataseries.
%    disp: text if k<0.2 'Not chaotic' else if k>0.9 'Chaotic' else 'Indecisive'.
%
%  
%   See also stabil, lyapunov, nmaxima
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('z1test')">commands z1test</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function results=z1test(x)

% This function is adapted from code from Paul Matthews:
% Copyright Paul Matthews, July 2009. (adapted by E. van Nes)
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

global g_grind g_Y;
if nargin==0||numel(x)==1||ischar(x)
    %if no arguments are given (or x is an integer value denoting the number of days)
    %we use the last GRIND run (g_Y) for x. If the settings have been changed, we
    %update the run
    if nargin==1
       ndays=i_checkstr(x);
       if ~isempty(ndays)
           g_grind.ndays=ndays;
       end
    end
    if isnan(g_grind.tstep)
        g_grind.tstep=g_grind.ndays;
    end
    N0 = i_initvar;
    if i_settingschanged(N0)
        time('-s')
    end
    x=g_Y;
end
if size(x,1)==1 
    x=transpose(x); 
end
if size(x,2)>1
    %if x is a matrix we take the column with the maximum avalue)
    x=x(:,max(x)==max(max(x)));
end
N=length(x);
meanx=mean(x);
%check if there is a fixed point 
%the 0-1 test is often indicisive in fixed points
if meanx>1E-6
    cv=std(x)/meanx;
else
    cv=std(x);
end
if cv<1e-4
    disp('Constant (k=0)');
    if nargout>0
        results.disp='Constant';
        results.k=0;
        results.N=N;
    end
    return
end
%the test does not work if the points are too densely packed.
%here we find local maxima and leave on average at most 10 points between maxima 
%but we keep minimal 10% of the points.
%
%get the maxima;
nmaxima = sum(diff( sign( diff([0; x(:,1); 0]) ) ) < 0 );
n_between_max=10;
if nmaxima<N/n_between_max
    %remove data if they are too dense packed (maximum lag =10)
    lag=max(10,round(N/(nmaxima*n_between_max))+1);
    x=x(1:lag:end);
    s=sprintf('Removed data by taking a sampling period of %d\n',lag);
    if nargout==0
        disp(s);
    else
        results.warning=s;
    end
    N=length(x);
end
 
%%%% The 1-0 test for chaos starts here
niter=100;
j=transpose(1:N);
t=transpose(1:round(N/10));
M=zeros(round(N/10),1);
c=pi/5+rand(1,niter)*3*pi/5;  % 100 random c values in [pi/5,4pi/5]
meanx2=mean(x)^2; %precalculate for speed
p=zeros(N,niter);
q=zeros(N,niter);
kcorr=zeros(1,niter);
for its=1:niter
   p(:,its)=cumsum(x.*cos(j*c(its)));
   q(:,its)=cumsum(x.*sin(j*c(its)));
   for n=1:round(N/10)
      %improved convergence formula from Gottwald & Melbourne 2009
      M(n)=mean( (p(n+1:N)-p(1:N-n)).^2 + (q(n+1:N)-q(1:N-n)).^2 )- ...
           meanx2*(1-cos(n*c(its)))/(1-cos(c(its)));
   end
   kcorr(its)=corr(t,M);
end

kmed=median(kcorr);
if nargout==0
    if kmed>0.9
        fprintf('Chaotic (k=%g)\n',kmed);
    elseif kmed<0.1
        fprintf('Not chaotic (k=%g)\n',kmed);
    else 
        fprintf('Indecisive (k=%g)\n',kmed);
    end
else
    results.N=N;
    results.k=kmed;
    results.ks=kcorr;
    results.cs=c;
    results.p=p;
    results.q=q;
    results.disp='Indecisive';
    if kmed>0.9
         results.disp='Chaotic';
    elseif kmed<0.1
        results.disp='Not chaotic';
    end
end

        
