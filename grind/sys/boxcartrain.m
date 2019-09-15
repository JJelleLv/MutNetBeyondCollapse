%BOXCARTRAIN   Implementation of a boxcar train (Goudriaan)
%   The boxcar train concept is a way of modelling age structured models (Goudriaan, 1988). Three kinds 
%   of boxcar trains are implemented in GRIND:
%   fixed boxcar train - this kind of boxcar train has a fixed variation (dispersion) due to numerical
%   dispersion
%   escalator boxcar train - this kind of boxcar train has no variation (dispersion).
%   fractional boxcar train - this version has an adjustable variation
%   This function is used in the definition of models, the solver should be <a href="matlab:help euler">Euler</a>.
%
%
%   Usage:
%   RES=BOXCARTRAIN(VAR,DEVRATE,CV) - VAR is a vector state variable that contains the boxcars. DEVRATE 
%   is the development rate (per time unit). CV is the desired relative variation (dispersion): CV=-1 : fixed BT; 
%   CV=0 : Escalator BT and CV>0: fractional BT; The result RES contains the flow (RES.flow) and the outflow 
%   per time unit (RES.outflow)
%   RES=BOXCARTRAIN(VAR,DEVRATE) - if CV is omitted the escalator boxcar train is assumed (CV=0)
%
%
%   See also boxcarinflow, modelpanel, model, solver
%
%   Reference page in Help browser:
%      <a href="matlab:commands('boxcartrain')">commands boxcartrain</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function [res,Anew] = boxcartrain(boxnr,A, devrate, cv_devrate)
global g_grind;
if nargin < 4
   %if no sd at default a escalator boxtrain;
   cv_devrate = 0;
end
%initiations;
%
deltat = solver('step');
N = size(A,1);
N2= size(A,2);
%gamma = 1 / (N-0.5); %%%% should last a bit longer as the last class transfers gradually into the next stage
 gamma = 1 / N; %%%% origial code of Goudriaan!!
%flow calculation
%
if isempty(cv_devrate) || (cv_devrate < 0)
   %
   %fixed boxcar train
   %
   %flows from each boxcar
   res.flow = A ./ gamma .* devrate;
   %shift the flows to get the inflow to the next boxcar and keep outflow
   res.outflow = res.flow(N,:);
   res.flow = [0; res.flow(1:N - 1,:)]-res.flow; %nett flow: what goes in what out?
elseif cv_devrate == 0
   %
   %escalator boxcar train
   %     
   res.flow = zeros(N, N2);
   res.outflow = zeros(1,N2);
   if isnan(g_grind.boxcar.gcycl(boxnr))
      g_grind.boxcar.gcycl(boxnr) = 0.5 * gamma; % !! boxcar.gcycl should be reset to []
   else
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) + devrate * deltat;
   end 
%the last class empties gradually instead of all at once
% Goudriaan code
% if (gamma-boxcar.gcycl > 1E-30)
 %    res.outflow =  min(A(N)/deltat, A(N) / (gamma - boxcar.gcycl) * devrate);
%      res.flow(N) = -res.outflow;  
% end

   if g_grind.boxcar.gcycl(boxnr) >= gamma
      %flow is the difference between next and current
      res.outflow = res.outflow+A(N,:)/deltat; % all biomass is moved in one moment to the next stage
      res.flow = ([zeros(1,N2); A(1:N - 1,:)]-A) / deltat; %divide by deltat such that the whole A will be shifted   
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) - gamma;
   end  
else
   %
   %fractional boxcar train
   %
   res.flow = zeros(N, N2);
   if isnan(g_grind.boxcar.gcycl(boxnr))
      g_grind.boxcar.gcycl(boxnr) = 0.5 * gamma; % !! g_grind.boxcar.gcycl should be reset to NaN
   else
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) + devrate * deltat;
   end
   %
   %Fraction F to be moved (Goudriaan)
   F = 1 - (N * cv_devrate^2);
   if F<0
      error('GRIND:boxcartrain:CVtoobig','Boxcartrain: coefficent of variation (%g) too big or number of stages (%d) too big\nboxcar %s CV should at least be < %g (<sqrt(1/N))',cv_devrate,N,g_grind.boxcar.names{boxnr},sqrt(1/N));
   end
   if deltat > F * gamma / (devrate + 1e-30)
      error('GRIND:boxcartrain:IntSteptoobig','Boxcartrain: Integration step too large: deltat should be < %0.5g\nAlternatively decrease the number of classes of %s to %d\ncv = %g devr = %g',...
         F * gamma / (devrate + 1e-30),g_grind.boxcar.names{boxnr},floor(1/(deltat*devrate+cv_devrate^2)),cv_devrate,devrate);
   end
   if (gamma -g_grind.boxcar.gcycl(boxnr) > 1E-30)
      res.outflow = min(A(N,:)/deltat,A(N,:) / (gamma - g_grind.boxcar.gcycl(boxnr)) * devrate);
      res.flow(N,:) = -res.outflow;
   else
      res.outflow = zeros(1,N2);
   end
   if g_grind.boxcar.gcycl(boxnr) >= F * gamma
      %flow is the difference between next and current
      res.flow = ([zeros(1,N2); A(1:N - 1,:)]-A) .* F ./ deltat; %divide by deltat such that the whole A will be shifted
      g_grind.boxcar.gcycl(boxnr) = g_grind.boxcar.gcycl(boxnr) - F * gamma;
   end
end
%update boxcar train
Anew = A + res.flow * deltat;
%g_grind.boxcar.trains{boxnr} = boxcar;
