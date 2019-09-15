% roundconv2 = convolution with periodical boundaries
%
% C=roundconv2(A,B) B is smaller matrix
% see also conv2
%
function C1=roundconv2(A,B)
sizA=size(A);
sizB=size(B);
A1=roundexp(A,sizA+2*sizB);
C=conv2(A1,B,'same');
C1=C(sizB(1)+1:sizB(1)+sizA(1),sizB(2)+1:sizB(2)+sizA(2));
