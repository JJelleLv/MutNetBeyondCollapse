%CONNECTMAT - matrix with connections between neighbors
%  This function returns a sparse matrix defining the connections between a matrix 
%  the matrix gives the connections between the elements of the space matrix
%  if A= the space matrix and M the connections then
%  d*(M*A(:)) defines the diffusion between 4 neighbors.
%
% Usage:
%   M=CONNECTMAT([SIZX,SIZY],ISROUND) SIZX is the size in X direction, SIZY is the size in the Y direction, ISROUND = 1 for 
%   periodic 0 for bouncing boundaries;
%   CONNECTMAT('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'border' [periodic | neumann | dirichlet or integer>0] - the kind of border in space
%     'matrix' [number] - matrix for defining the size
%     'setdiag' [logical] - set the diagonal to negative sum of neighbors
%     'size' [integer>1] - size of the matrix.
% 
%  See also neighborcells, leftcells, rightcells, upcells, downcells     
%
%   Reference page in Help browser:
%      <a href="matlab:commands('connectmat')">commands connectmat</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function [M]=connectmat(varargin)
%(s1,s2,isround)
fieldnams={'size','i>1','size of the matrix.',[];...
   'matrix','n','matrix for defining the size',[];...
   'border','e[periodic|neumann|dirichlet]#i>0','the kind of border in space','periodic';...
   'setdiag','l','set the diagonal to negative sum of neighbors',[]}';
args=i_parseargs(fieldnams,'matrix,border','',varargin);
if isfield(args,'matrix')
    if ischar(args.matrix)
        args.matrix=evalin('base',args.matrix);
    end
    if numel(args.matrix)>2
       args.size=size(args.matrix);
    else
       args.size=args.matrix;
    end
end
if ~isfield(args,'size')
    error('GRIND:connectmat:UnknownSize','What is the size of the connections matrix?');
elseif numel(args.size)==1
    args.size=[args.size args.size];
end
if ~isfield(args,'border')
    args.border= 'neumann';
end
M=mat2connect(args.size,@neighborcells,4,args.border);
psiz=prod(args.size);
if ~isfield(args,'setdiag')||args.setdiag
   M(1:psiz+1:psiz*psiz) = M(1:psiz+1:psiz*psiz)-sum(M);
end

  
