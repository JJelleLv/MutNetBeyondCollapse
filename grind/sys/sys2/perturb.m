%PERTURB   Perturb saddle point
%   perturb a saddle point in the direction of the attracting eigenvector
%   use <a href="matlab:help backw">backw</a> after perturbation of the equilibrium to approximate one part of
%   the separatrix
%
%   To find the whole separatrix:
%   perturb 1
%   backw 30
%   perturb -1
%   backw 30
%  
%   Usage:
%   PERTURB 1 - perturb along one unstable eigenvector
%   PERTURB -1 PERSIZE - perturb in the opposite direction of the first unstable eigenvector. 
%        PERSIZE is the size of the perturbation
%   PERTURB('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'dir' [integer>=-2 and integer<=2] - The direction of the perturbation: values 1/-1 = stable manifolds; 2/-2 = unstable manifolds
%     'persize' [number>0] - the size of  the perturbation (default=0.005)
%      
%
%   Usage:
%   PERTURB - Perturb the equilibrium in positive direction with a size of
%   0.05
%   PERTURB -1 SIZE - direction can be 1 or -1. SIZE should be a small
%   number.
%
%See also backw, findeq, optimlib/perturb
%
%   Reference page in Help browser:
%      <a href="matlab:commands('perturb')">commands perturb</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [pert,x] = perturb(varargin)
%(dir, persize,u,del)
if (nargin>2)||(nargout>0)&&length(which('perturb','-all'))>1
   d=which('perturb','-all');
   ndx = strcontains(d, 'toolbox');
   d = d(~ndx);
   if ~isempty(d)
      cur = pwd;
      cd(fileparts(d{1}));
      [pert, x] = perturb(varargin{:});
      cd(cur);
      return;
   end

end
global g_grind;
fieldnams={'dir', 'i>=-2&i<=2', 'The direction of the perturbation: values 1/-1 = stable manifolds; 2/-2 = unstable manifolds',1;...
   'persize', 'n>0', 'the size of  the perturbation (default=0.005)',0.005}';
args=i_parseargs(fieldnams,'dir,persize','',varargin);
i_parcheck;
if ~isfield(args,'dir')
   args.dir = 1;
end

if ~isfield(args,'persize')
   args.persize = 0.005;
end

[N1] = findeq('Display','off');
[~, eigenvalues, eigenvect] = i_eigen(1);
if abs(args.dir) > 10
   vec = abs(args.dir) - 10;
   args.dir = args.dir / abs(args.dir);
else
   if i_stability(eigenvalues(1), g_grind.solver.isdiffer)
      vec = 1;
   elseif i_stability(eigenvalues(2), g_grind.solver.isdiffer)
      vec = 2;
   else
      vec = 1;
   end

end

if abs(args.dir) == 2
   args.dir = args.dir / 2;
   if vec == 1
      vec = 2;
   else
      vec = 1;
   end

end

% len = length(real(eigenvect(:, vec)));
% if len < 0.1
%    len = 0.1;
% end

N0 = N1 + args.dir * real(eigenvect(:, vec)) / norm(real(eigenvect(:, vec))) * args.persize;
i_keep(N0);

