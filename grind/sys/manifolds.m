%MANIFOLDS   Trajectories in the direction of the eigenvectors
%   In the current phase plane the manifolds of the currently selected equilibrium
%   (see: <a href="matlab:help findeq">findeq</a>) are drawn. This is especially useful for saddle points, 
%   as the stable manifolds of these equilibria are the separatrix.
%
%   Usage:
%   
%   MANIFOLDS - draw all manifolds of the currently selected equilibrium
%   MANIFOLDS('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%  Order of normal arguments:
%     1: 'persize'
%     2: 'stable'
%     'persize' [number>0] - the size of  the perturbation (see <a href="matlab:help perturb">perturb</a>) (default=0.005)
%     'stable' [integer>=0 and integer<=2] - draw the stable manifolds (1) or unstable (0) or both (2) (default=2)
%   MANIFOLDS('-opt1','-opt2',...) - Valid command line options:
%     '-s' - stable: draw all stable manifolds (=separatrix).
%     '-u' - unstable: draw all unstable manifolds
%
%
%   See also findeq, perturb, ru, backw
%
%   Reference page in Help browser:
%      <a href="matlab:commands('manifolds')">commands manifolds</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res1=manifolds(varargin)
%(stable)
global g_grind g_Y;
fieldnams={'persize', 'n>0', 'the size of  the perturbation (see <a href="perturb.htm">perturb</a>) (default=0.005)',0.005;...
   'stable', 'i>=0&i<=2', 'draw the stable manifolds (1) or unstable (0) or both (2) (default=2)',2}';
[args,defaults]=i_parseargs(fieldnams,'persize,stable','-s,-u',varargin);
args=mergestructs(defaults,args);
if any(strcmp(args.opts,'-s'))
    args.stable =1;
elseif any(strcmp(args.opts,'-u'))
    args.stable = 0;
end
if g_grind.solver.haslags
    error('grind:manifolds:lags','The command "manifolds" cannot handle delay differential equations')
end
res=[];
[N1,found] = findeq('Display','off','maxiter',1,'tolfun', 1E-5,'tolx',1E-8);
if ~found
    error('grind:manifolds:noeq','First select an equilibrium (e.g. with <a href="matlab:findeq">findeq</a>) before using "manifolds"' )
end
[~, eigenvalues,  eigenvect] = i_eigen(1);
domeigen=find(eigenvalues == max(eigenvalues));
[isstable, issaddle]  = i_stability(eigenvalues, g_grind.solver.isdiffer);
if issaddle
   for i = 1:length(eigenvalues)
      stab = i_stability(eigenvalues(i), g_grind.solver.isdiffer);
      if stab && (args.stable~=0)
         i_keep(N1 + eigenvect(:, i) *  args.persize);
         backw;
         res.stable1=g_Y;
         i_keep(N1 - eigenvect(:, i) *  args.persize);
         backw;
         res.stable2=g_Y;
      elseif ~stab && (args.stable~=1)
         i_keep(N1 + eigenvect(:, i) *  args.persize);
         ru;
         res.unstable1=g_Y;
         i_keep(N1 - eigenvect(:, i) *  args.persize);
         ru;
         res.unstable2=g_Y;
      end
   end
else
   if isstable
      domeigen=find(eigenvalues == min(eigenvalues));
      i_keep(N1 + real(eigenvect(:, domeigen)) *  args.persize*10);
      backw;
      res.stable1=g_Y;
      i_keep(N1 - real(eigenvect(:, domeigen)) *  args.persize*10);
      backw;
      res.stable2=g_Y;
      disp('Can only show manifold of the lowest eigenvalue')
   else
      i_keep(N1 + real(eigenvect(:, domeigen)) *  args.persize);
      ru;
      res.unstable1=g_Y;
      i_keep(N1 - real(eigenvect(:, domeigen)) *  args.persize);
      ru;
      res.unstable2=g_Y;
      disp('Can only show manifold of the dominant eigenvalue')
   end
end
i_keep(N1)
if nargout>0
    res1=res;
end


