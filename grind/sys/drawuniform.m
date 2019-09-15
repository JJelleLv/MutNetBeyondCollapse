%DRAWUNIFORM  Draw from an uniform distribution
%  Convenient alternative for rand.
% 
%  Usage:
%  A=DRAWUNIFORM(10,0.5,0.7) - draw a 10x10 matrix with random numbers between 0.5 and 0.7
%  A=DRAWUNIFORM(10,20,0,4) - draw a 10x20 matrix between 0 and 4
% 
%  See also <a href="matlab:help rand">rand</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('drawuniform')">commands drawuniform</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function A=drawuniform(dim1,dim2,minA,maxA)
if nargin==3
   maxA=minA;
   minA=dim2;
   dim2=dim1;
end
if isempty(minA)||isempty(maxA)||isnan(minA)||isnan(maxA)
    error('GRIND:drawuniform:ArgError','Error drawuniform: range [%g,%g] is not correct\n',minA,maxA);
end
A=rand(dim1,dim2)*(maxA-minA)+minA;
