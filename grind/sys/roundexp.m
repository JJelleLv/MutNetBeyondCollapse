% Expand a matrix periodically
% new size should be larger than old size
% at most doubled to all sides   
function A1=roundexp(A,newsize)
newsize=newsize-size(A);
n1=floor(newsize(1)/2);
n2=newsize(1)-n1;
A1=[A(end-n2+1:end,:);A;A(1:n1,:)];
n1=floor(newsize(2)/2);
n2=newsize(2)-n1;
A1=[A1(:,end-n2+1:end),A1,A1(:,1:n2)];