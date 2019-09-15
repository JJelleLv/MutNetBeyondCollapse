%BOXCARINFLOW   Inflow into a boxcartrain
%   This simple function is used to add biomass to a boxcartrain in a model definition. It returns a vector 
%   with all zero’s except for the first element, which is filled with the inflow.
%
%
%   Usage:
%   res=BOXCARINFLOW(VAR, INFLOW) - returns a zero vector of sizeof(VAR) with inflow in the first element.
%
%   See also boxcartrain, model, insim
%
%   Reference page in Help browser:
%      <a href="matlab:commands('boxcarinflow')">commands boxcarinflow</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function res=boxcarinflow(A,inflow)
res=zeros(size(A));
if numel(inflow)>1
    warning('grind:boxcarinflow','Multiple inflows to boxcar are summed');
    inflow=sum(inflow(:));
end
res(1,:)=inflow;

