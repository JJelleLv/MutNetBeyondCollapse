%REDNOISE   Generate red noise
%   Function to generate white or red noise. To be used in differential
%   equations. Itself is a difference equation (Steele and Henderson):
%   T(t) = (1 -1/lambda) * (T(t-1) - T0) + T0 + beta*randn
%
%   In which:
%   T = variable with red noise
%   T0 = mean of red noise
%   lambda = 'period' of red noise. white: lambda=1 red: lambda>1.
%   beta = extend of noise
%   randn = normally distributed random number
%
%   Note: you can set the standard deviation of the resulting series to a certain value
%   with the following formula:(Ives et al., 2003):
%   beta=SD*sqrt(2/lambda-1/lambda^2) (in which SD is the resulting standard
%   deviation).
%
% 
%   Usage:
%   REDNOISE(t,T0,lambda,beta) for explanation of coefficients, see above (t=time).
%   REDNOISE(t,T0,lambda,beta,iset) - If you use several independent rednoise sets 
%   within one set of equations, number them using an integer starting with 1 (iset).
%   The iset can be a vector with numbers, rednoise will then return a vector.
%   REDNOISE(t,T0,lambda,beta,iset,deltat) - define the period deltat for which the
%   difference equation will be calculated (default=1)
%   REDNOISE -D - (from the command line) deactivates all rednoise sets, so the last generated
%   set will be used even if parameters change.
%   REDNOISE -A - (from the command line) activates rednoise again.
%   
%       
%   See also model, externvar, dwiener
%
%   Reference page in Help browser:
%      <a href="matlab:commands('rednoise')">commands rednoise</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function T = rednoise(tnow, T0, labda, beta, isets, deltat)
global g_noise t;
if nargin<=2
    if strncmpi(tnow,'-a',2)
        global g_grind;  %#ok<TLEV>
        if ~isfield(g_noise,'active')||~g_noise.active
           g_noise.active=1;
           g_noise.y={};
           g_grind.checks.lastsettings=[];
           disp('rednoise activated');
        else
           disp('rednoise was already active')
        end
     elseif strncmpi(tnow,'-d',2)
        g_noise.active=0;
        disp('rednoise deactivated, note that rednoise is even not updated if its parameters change');
    elseif strncmpi(tnow,'-u',2)
        if ~isfield(g_noise,'active')||g_noise.active
          g_noise.y={};
          g_noise.isets=[];
          g_noise.deltat=[];
          g_noise.active=1;
        end   
    end
    return;
end
if nargin < 5
   isets = 1;
end
if nargin < 6
   deltat = 1;
end
%update = min(isets) > 0;
if isets == 0
   isets = 1;
end
if numel(tnow)==1 %this gives a mistake if isets is not specified
   maxparlen=max([size(T0,1), size(labda,1),size( beta,1)]);
   if length(isets) < maxparlen
      if tnow==t
          disp('Warning rednoise: Size of isets is adapted as one of the parameters is a vector, this may give problems');
      end
      isets = transpose(isets(1):isets(1) + maxparlen - 1);
   end
end
if isempty(g_noise)
    g_noise.y={};
    g_noise.isets=[];
    g_noise.deltat=[];
end    
istart=find(g_noise.isets==isets(1),1);
if isempty(istart)
    istart=length(g_noise.isets)+1;
end
if length(isets)>1
   zer=zeros(length(isets),1);
   deltat=deltat + zer;
   labda=labda + zer;
   T0=T0+zer;
   beta=beta + zer;
end
if istart>length(g_noise.y)
    g_noise.deltat(istart:istart-1+length(isets))=deltat;
    g_noise.isets(istart:istart-1+length(isets))=isets;
    g_noise.y(istart:istart-1+length(isets))={[]};
end
T= zeros(length(tnow),length(isets));
for k = 1:length(tnow)
   if beta == 0
      g_noise.y{istart} = [];
      T=T0+T;
      if size(tnow,2)>1
         T=transpose(T);
      end
      return;
   end
   for i = 1:length(isets)
      iset = istart-1+i;
      if tnow(1)<t
          sdeltat=-1;
      else
          sdeltat=1;
      end
      nt=(round((tnow(k) - t) / g_noise.deltat(iset)*sdeltat) + 1);
      if (nt>length(g_noise.y{iset}))
           g_noise.y{iset}=[g_noise.y{iset}; zeros(1000,1)+NaN];
      end
      if isnan(g_noise.y{iset}(nt))
         ifirst=find(isnan(g_noise.y{iset}),1);
         if ifirst==1
             g_noise.y{iset}(1)=T0(i);
             ifirst=2;
         end
         for j=ifirst:nt
             g_noise.y{iset}(j) = (1 - 1 / labda(i)) * (g_noise.y{iset}(j-1) - T0(i)) + T0(i) + beta(i) * randn;
         end
      end
      if (tnow(k)*sdeltat) >= t
         T(k, i) = g_noise.y{iset}(nt);
      end
   end
end
if size(tnow,2)>1
    T=transpose(T);
end
