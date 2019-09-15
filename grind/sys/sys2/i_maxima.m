function [maxima,ndx] = i_maxima(A, iX)
ndx=diff(sign([0;diff(A(:,iX));0]))<0;
maxima = A(ndx,iX);

%maxima = ones(size(A, 1), size(A, 2));
%imax = 0;
%siz = size(A, 1);
%for i = 2:siz - 1
%   if (A(i - 1, iX) < A(i, iX)) && (A(i, iX) > A(i + 1, iX));
%      imax = imax + 1;
%      maxima(imax, :) = A(i, :);
%   end

%end

%maxima = maxima(1:imax, :);
