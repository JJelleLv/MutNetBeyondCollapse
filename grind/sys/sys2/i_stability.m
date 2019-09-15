function [isstable, issaddle, isspiral, domcomplex] =  i_stability(eigenval, isdiffer, haslags)
if nargin < 3
   haslags = false;
end

if isdiffer
   isstable = max(abs(eigenval), [], 1) < 1;
   issaddle = (min(abs(eigenval), [], 1) < 1) & (max(abs(eigenval), [], 1) > 1); %No convert
else
   %    if size(eigenval, 1) == 1
   %       isstable = real(eigenval) < 0;
   %       issaddle = (real(eigenval) > 0) & (real(eigenval) < 0); %No convert
   %    else
   isstable =  max(real(eigenval), [], 1) < 0;
   issaddle =  (max(real(eigenval), [], 1) > 0) & (min(real(eigenval), [], 1) < 0); %No convert
   %   end

end

if size(eigenval, 1) == 1
   isspiral = abs(imag(eigenval)) > 0;
else
   if haslags
      ismax = real(eigenval) < repmat(max(real(eigenval)), size(eigenval, 1), 1); %only the maximum eigenvalues are important
      eigenval1 = eigenval;
      eigenval1(ismax) = -1;
      isspiral = max(imag(eigenval1)) > 0;
   else
      isspiral = max(imag(eigenval)) > 0;
   end

end
domcomplex = zeros(size(issaddle)) > 0;
if any(isspiral)
   siz = size(eigenval);
   if min(siz) == 1
      domeig = eigenval;
   else
      [~, ii] =  max(real(eigenval));
      domeig = eigenval(sub2ind(siz, ii, 1:siz(2)));
   end

   domcomplex = imag(domeig) > 0;
end

