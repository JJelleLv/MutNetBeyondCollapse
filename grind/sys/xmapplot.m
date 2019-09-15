%CAUSCORR - Convergent Cross Mapping of Sugihara et al., 2012
%A MATLAB implementation of convergent cross mapping to detect causation
%in connected time series. This function does the cross maps for a range of
%lengths of the library
%
%Citation:
%Sugihara, G., R. May, H. Ye, C. H. Hsieh, E. Deyle, M. Fogarty, and
%S. Munch. 2012. Detecting causality in complex ecosystems. Science 338:496-500.
%
%Usage:
%[corrs] = causcorr(X,Ydata, E, Ls, b)  calculates cross maps of X on all columns
%  of Ydata. (so the effect of each column of Y is evaluated on X.
%  Preferrable X and Ydata are datasets (see Statistics Toolbox). E =
%  embedding dimension, Ls is a vector with the lengths of the library,
%  b(optional) is the number of points used to map Y (by default b=E+1)
%
%
function xmapplot(tau, E, Ls, b)
global g_Y g_t t g_grind g_xmap;
evalin('base','global g_xmap;')
i_parcheck;
if nargin < 1
   tau = 10;
end
if nargin  < 2
   E = 3;
end
if nargin < 3
   Ls = 5:100;
end
if nargin < 4
   b = E + 1;
end
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(t, g_grind.ndays, N0, 1);
   g_xmap=[];
end
if xmapchanged(tau, E, Ls, b)
   g_xmap.tau = tau; g_xmap.E = E;  g_xmap.b = b;
   X = interp1(g_t, g_Y, g_t(1):tau:g_t(end));
   Y = [X, [X(2:end, :); nan(1, size(X, 2))]];
   %Y=X;
   Xnames = g_grind.statevars.names;
   Ynames = [g_grind.statevars.names g_grind.statevars.names];
   for i = length(Xnames) + 1:length(Ynames)
      Ynames{i} = sprintf('%s(t+%d)', Ynames{i}, round(tau));
   end
   %Ynames=Xnames;
   nY = length(Ynames);
   nX = length(Xnames);
   Ls = Ls(Ls < length(X));
   g_xmap.L = transpose(Ls);
   corrs = zeros(nX,  nY,length(Ls));
   i_waitbar(0, nX*length(Ls), 'xmapplot', 'Calculating',0.5)

   for j = 1:nX
      for i = 1:length(Ls)
         i_waitbar(1);
         res1 = CCM(X(:,j), Y, E, Ls(i), 1, b);
         corrs(j, :, i) = res1.rho;
      end
   end
   i_waitbar([]);
   names = cell(nX, nY,1);
   k = 1;
   for j = 1:length(Xnames)
      for i = 1:length(Ynames)
         names{j,i,1} = sprintf('%s xmap %s', Xnames{j}, Ynames{i});
         k = k + 1;
      end
   end
   g_xmap.corrs = transpose(reshape(corrs, nY * nX, length(Ls)));
   g_xmap.names = transpose(reshape(names, nY * nX, 1));
end

hfig=figure(998);
n=length(g_xmap.names);
plot(g_xmap.L, g_xmap.corrs(:, 1:n / 2))
i_plotdefaults(hfig);
xlabel('Length of library L')
ylabel('\rho')
legend(g_xmap.names(1:n / 2))

hfig=figure(999);
plot(Ls, g_xmap.corrs(:, n/2 + 1:n))
i_plotdefaults(hfig);
xlabel('Length of library L')
ylabel('\rho')
legend(g_xmap.names(n/ 2 + 1:n))

function res = xmapchanged(tau,E, Ls, b)
global g_xmap
res = isempty(g_xmap)||(g_xmap.tau~=tau)||(g_xmap.E~=E)||(g_xmap.b~=b);
if ~res
   for i = 1:length(Ls)
      if (i > length(g_xmap.L))||(Ls(i)~=g_xmap.L(i))
         res = 1;
         return;
      end
   end
end
