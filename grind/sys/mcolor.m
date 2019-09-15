%function mcolor Checkerboard plot of matrix
%the same as PCOLOR but with all rows and columns
%a nan row is added to the right and bottom of the matrix
%See PCOLOR
function h=mcolor(X,Y,C,p)
if ishandle(X)
    hax=X;
    if nargin==2
       C=Y;
       [X,Y]=meshgrid(1:size(C,2),1:size(C,1));
    elseif nargin==4
        X=Y;
        Y=C;
        C=p;
    end
else
    hax=gca;
end
if nargin==1
    C=X;
    [X,Y]=meshgrid(1:size(C,2),1:size(C,1));
end
C1=nan(size(C)+1);
C1(1:end-1,1:end-1)=C;
difx1=mean(diff(X(1,:)));
difx2=mean(diff(X(:,1)));
X1= [X X(:,end)+difx1];
X1= [X1;X1(end,:)+difx2];
X1=X1-max(difx1,difx2)/2;
dify1=mean(diff(Y(1,:)));
dify2=mean(diff(Y(:,1)));
Y1= [Y Y(:,end)+dify1];
Y1= [Y1;Y1(end,:)+dify2];
Y1=Y1-max(dify1,dify2)/2;
h1=pcolor(hax,X1,Y1,C1);
if nargout>0
    h=h1;
end