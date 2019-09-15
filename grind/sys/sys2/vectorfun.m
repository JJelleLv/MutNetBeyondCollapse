%i_vector internal function to create a matrix with vectors
function [Vect] = vectorfun(nlines,npoints, thefun)
global g_grind;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
Xaxis =  g_grind.xaxis.lim;
Yaxis =  g_grind.yaxis.lim;
if ~isempty(nlines) && (nlines > 1)
   iZ = i_getno(g_grind.zaxis.var);
   Zaxis =  g_grind.zaxis.lim;
   oldZ = i_getparvar(iZ);
   if iZ.isext
      oldactZ = g_grind.externvars{iZ.no}.options.active;
      g_grind.externvars{iZ.no}.options.active = 0;
   end
   N0 = i_initvar;
   try
      Vect = cell(1, nlines);
      Vect{1}  = zeros(npoints + 1, npoints + 1, nlines);
      for z1 = 1:nlines
         pz = (z1 - 1) * (Zaxis(2) - Zaxis(1)) / nlines + Zaxis(1);
         N0 = i_initvar;
         N0 =  i_setparvar(iZ, N0, pz);
         V1 = vectorfun1(Xaxis,Yaxis,iX,iY,npoints, thefun);
         %i_vector(npoints, iX, g_grind.xaxis.lim, iY, g_grind.yaxis.lim, N0);
         for i = 1:size(V1, 3)
            V = Vect{i};
            V(:, :, z1) = V1(:, :, i);
            Vect{i} = V;
         end

      end

      i_setparvar(iZ, N0, oldZ);
      if iZ.isext
         g_grind.externvars{iZ.no}.options.active = oldactZ;
      end
   catch err
      i_setparvar(iZ, N0, oldZ);
      if iZ.isext
         g_grind.externvars{iZ.no}.options.active = oldactZ;
      end
      rethrow(err);
   end

else
   Vect = vectorfun1(Xaxis,Yaxis,iX,iY,npoints, thefun);
end


function Vect = vectorfun1(Xaxis,Yaxis,iX,iY,npoints, thefun)
global g_grind;
N = i_initvar;
oldX = i_getparvar(iX);
oldY = i_getparvar(iY);
if iX.isext
   oldactX = g_grind.externvars{iX.no}.options.active;
   g_grind.externvars{iX.no}.options.active = 0;
end
if iY.isext
   oldactY = g_grind.externvars{iY.no}.options.active;
   g_grind.externvars{iY.no}.options.active = 0;
end
try
   %nvar = g_grind.statevars.dim;
   nvar = 1;
   Vect = zeros(npoints, npoints, nvar);
   incrY = (Yaxis(2) - Yaxis(1)) / npoints;
   incrX = (Xaxis(2) - Xaxis(1)) / npoints;
   for y1 = 1:npoints + 1
      py = (y1 - 1) * incrY + Yaxis(1);
      N = i_setparvar(iY, N, py);
      for x1 = 1:npoints + 1
         px = (x1 - 1) * incrX + Xaxis(1);
         N = i_setparvar(iX, N, px);
         i_keep(N);
         Nres = feval(thefun);
         Vect(y1, x1, :) = Nres;
         fprintf('(y=%g,x=%g) = %g\n',py,px,Nres);
      end

   end

   N = i_setparvar(iY, N, oldY);
   i_setparvar(iX, N, oldX);
   if iX.isext
      g_grind.externvars{iX.no}.options.active = oldactX;
   end
   if iY.isext
      g_grind.externvars{iY.no}.options.active = oldactY;
   end
catch err
   %  err=lasterror;
   N = i_setparvar(iY, N, oldY);
   i_setparvar(iX, N, oldX);
   if iX.isext
      g_grind.externvars{iX.no}.options.active = oldactX;
   end
   if iY.isext
      g_grind.externvars{iY.no}.options.active = oldactY;
   end
   rethrow(err);
end

