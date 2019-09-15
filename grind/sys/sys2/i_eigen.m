%i_EIGEN = calc Jacobian and eigenvalues of the current initial condition
%
function [Jacobian, eigenval, eigenvect] = i_eigen(donumerical,niters,N0)
global g_grind;
if nargin < 1
   donumerical = isempty(g_grind.syms.Jacobian);
end

if nargin < 3
   N0 = i_initvar;
end

if nargin < 2
   niters = g_grind.solver.iters;
end

if ~donumerical&&isempty(g_grind.syms.Jacobian)
    try
       enterjac('-sym')
    catch %#ok<CTCH>
        g_grind.syms.Jacobian=[];
    end
    if isempty(g_grind.syms.Jacobian)
        donumerical=1;
    end
end

if g_grind.solver.haslags
   %Delay differential equations are difficult
   J0 = i_calcjac(donumerical, niters, N0); %Jacobian for the non - lags
   Jtau = zeros(g_grind.statevars.dim, g_grind.statevars.dim, length(g_grind.dde.lags));
   tau = zeros(size(g_grind.dde.lags));
   for i = 1:length(g_grind.dde.lags)
      Jtau(:,:,i) = i_calcjac(donumerical,niters, N0, i); % Jacobian for the lags
      tau(i) = evalin('base', g_grind.dde.lags{i});
   end

   if g_grind.statevars.dim==1 && length(g_grind.dde.lags)==1 && i_hastoolbox('symbolic')
      %Eigenvalues in delay differential equations are complex
      %for the characteristic equation we need to solve
      %det(J0 + Jtau1*exp(-lambda*tau1)+ Jtau2*exp(-lambda*tau2) ... -lambda*I)=0
      %for one dimensional/one lag systems (lag tau):
      %J0+ Jtau*exp(-lambda*tau)-lambda=0
      %can only be solved using the lambert W function (or numerically)
      % lambda(k) = J0*tau + lambertw(k, (Jt*tau)/exp(J0*tau))/tau;
      % k=1 to Inf
      %
      eigenval = zeros(10, 1);
      for i = 1:10
         eigenval(i) = J0 * tau + lambertw(i - 6, (Jtau * tau) / exp(J0 * tau)) / tau;
      end

      Jacobian = J0 + Jtau *  exp(-eigenval(6) * tau);
      eigenvect = eigenval(6);
   elseif g_grind.statevars.dim==2 && length(g_grind.dde.lags)==1
      %characteristic equation for 2 variables 1 lag
      %solve for lambda: det(J0 + Jtau*exp(-lambda*tau) - lambda*I)=0
      %mupad:
      %> J0:=matrix(2,2,[[J0_11,J0_12],[J0_21,J0_22]]);
      %> Jtau:=matrix(2,2,[[Jtau_11,Jtau_12],[Jtau_21,Jtau_22]]);
      %> Jtot:=J0+Jtau*exp(-lambda*tau)-lambda*matrix(2,2,[[1,0],[0,1]])
      %> simplify(det(Jtot))
      %lambda^2 - J0_22*lambda - J0_11*lambda + J0_11*J0_22 - J0_12*J0_21 +
      %(J0_11*Jtau_22)/exp(lambda*tau) - (J0_12*Jtau_21)/exp(lambda*tau) -
      %(J0_21*Jtau_12)/exp(lambda*tau) + (J0_22*Jtau_11)/exp(lambda*tau) +
      %(Jtau_11*Jtau_22)/exp(2*lambda*tau) - (Jtau_12*Jtau_21)/exp(2*lambda*tau) -
      %(Jtau_11*lambda)/exp(lambda*tau) - (Jtau_22*lambda)/exp(lambda*tau)
      %
      %lambda^2 - F1*lambda + F2 + (F3 - F4*lambda) * exp(-lambda*tau) + F5*exp(-2*lambda*tau)
      F(1) = trace(J0); % = J0(2, 2) + J0(1, 1);
      F(2) = det(J0); % = J0(1, 1) * J0(2, 2) - J0(1, 2) * J0(2, 1);
      F(3) = J0(1, 1) * Jtau(2, 2) - J0(1, 2) * Jtau(2, 1) + J0(2, 2) * Jtau(1, 1) - J0(2, 1) * Jtau(1, 2); %kind of "cross determinant"
      F(4) = trace(Jtau); % = Jtau(1, 1) +  Jtau(2, 2);
      F(5) = Jtau(1, 1) * Jtau(2, 2)-Jtau(1, 2) * Jtau(2, 1); % = Jtau(1, 1) * Jtau(2, 2) - Jtau(1, 2) * Jtau(2, 1);
      F(6) = tau(1);
      eigs = zeros(100 * length(tau), 2);
      opt=struct('MaxIter',5000,'MaxFunEvals',5000,'TolX',1e-10,'TolFun',1e-10);
      for i = 1:100 * length(tau)
         lambdas = rand(1, 2) * 50;
         [aa1,xx] = fminsearch(@(lambdas)chardelay(lambdas,F),lambdas, opt);
         if xx < 1e-9
            eigs(i, :) = aa1;
         else
            eigs(i, :) = -9999;
         end

      end

      [~, ii] = unique(round(eigs(:, 1) * 10000) / 10000);
      eigs = eigs(ii, :);
      eigs=eigs(eigs(:, 1) ~= -9999, :);
      if size(eigs, 1) > 20
         eigs = eigs(end - 20:end, :);
      end

      eigenval = zeros(20, 1) + NaN;
      n = size(eigs, 1);
      eigenval(1:n) = eigs(:, 1) + eigs(:, 2) * 1i;
      Jacobian=J0 + Jtau *  exp(-eigenval(real(eigenval) == max(real(eigenval))) * F(6));
      [~, eigenvect] = eig(Jacobian);
      fprintf('Characteristic equation:\nlambda^2 - %g*lambda + %g + (%g-%g*lambda)*exp(-lambda*%g)+ %g*exp(-2*lambda*%g)=0\n', F(1),F(2),F(3),F(4),F(6),F(5),F(6));
      % error('grind:eigen','GRIND can currently not solve this equation');
   else
      %a bit slower but more versatile routine for delay differential
      %equaitons
      eigs = zeros(100 * length(tau), 2);
      opt=struct('MaxIter',5000,'MaxFunEvals',5000,'TolX',1e-10,'TolFun',1e-10);
      for i = 1:100 * length(tau)
         lambdas = rand(1, 2) * 50;
         [aa1,xx] = fminsearch(@(lambdas)chardelay1(lambdas,J0,Jtau,tau),lambdas, opt);
         if xx < 1e-9
            eigs(i, :) = aa1;
         else
            eigs(i, :) = -9999;
         end

      end

      [~, ii] = unique(round(eigs(:, 1) * 10000) / 10000);
      eigs = eigs(ii, :);
      eigs=eigs(eigs(:, 1) ~= -9999, :);
      d = g_grind.statevars.dim * 10;
      if size(eigs, 1) > d
         eigs = eigs(end - d + 1:end, :);
      end

      eigenval = zeros(d, 1) + NaN;
      n = size(eigs, 1);
      eigenval(1:n) = eigs(:, 1) + eigs(:, 2) * 1i;
      ndx=real(eigenval) == max(real(eigenval));
      Jacobian = J0;
      for i = 1:length(tau)
         Jacobian = Jacobian + Jtau(i) * exp(-eigenval(ndx) * tau(i));
      end

      [~, eigenvect] = eig(Jacobian);
   end

   %   else
   %      error('grind:eigen','Eigenvalues for multi-dimensional or multi-lag models not supported');
   
else
   %normal models
   Jacobian = i_calcjac(donumerical, niters,N0);
   if any(any(isnan(Jacobian)))||any(any(isinf(Jacobian)))
       warning('grind:eigen','Cannot determine eigenvalues as there are NaN values in the Jacobian')
       eigenval=nan(size(Jacobian,1),1);
       eigenvect=nan(size(Jacobian));
   else
       [eigenvect, eigenval] = eig(Jacobian);
       eigenval = diag(eigenval);
   end

end

%eigenval=eig(Jacobian);
function res = chardelay(lambdas, F)
lambda = lambdas(1) + lambdas(2) * 1i;
res = lambda^2 - F(1)*lambda + F(2) + (F(3) - F(4)*lambda) * exp(-lambda*F(6)) + F(5)*exp(-2*lambda*F(6));
res = abs(real(res)) + abs(imag(res));
function res = chardelay1(lambdas, J0, Jtau, tau)
lambda = lambdas(1) + lambdas(2) * 1i;
res = J0 - eye(size(J0)) * lambda;
for i = 1:length(tau)
   res = res + Jtau(:, :, i) * exp(-lambda * tau(i));
end

res = abs(det(res));

