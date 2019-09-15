function plotcomplexfnc(fnc,rangeReal,rangeIm,npoints)
if nargin<1
    fnc=@(x)lambertw(x);
end

if nargin<2
    rangeReal=[-pi pi];
end

if nargin<3
    rangeIm=[-pi pi];
end

if nargin<4
    npoints=20;
end

close all;
figure(1);
cs=zeros(npoints,npoints,3);
c=[0 0 0];
for i=1:npoints-1
    for j=1:npoints-1
        x=([i i i+1 i+1]-1).*(rangeReal(2)-rangeReal(1))/(npoints+1)+rangeReal(1);
        y=([j j+1 j+1 j]-1).*(rangeIm(2)-rangeIm(1))/(npoints+1)+rangeIm(1);
        if x(1)>0&&y(1)>0
            c(1)=(i-npoints/2)/(npoints/2);
            c(2)=(j-npoints/2)/(npoints/2);
            c(3)=0.5;
        elseif x(1)<0&&y(1)>0
            c(1)=1;
            c(2)=i/(npoints/2);
            c(3)=(j-npoints/2)/(npoints/2);
        elseif x(1)>0&&y(1)<0
            c(1)=(i-npoints/2)/(npoints/2);
            c(2)=1;
            c(3)=(j)/(npoints/2);
        else 
            c(1)=(i)/(npoints/2);
            c(2)=(j)/(npoints/2);
            c(3)=1;
        end

        c(c>1)=1;
        c(c<0)=0;
        cs(i,j,:)=c;
        patch(x,y,cs(i,j,:));
       % set(h,'Edgecolor','none')
    end

end

hold on;
figure(2);
for i=1:npoints-1
    for j=1:npoints-1
        x=([i i i+1 i+1]-1).*(rangeReal(2)-rangeReal(1))/(npoints+1)+rangeReal(1);
        y=([j j+1 j+1 j]-1).*(rangeIm(2)-rangeIm(1))/(npoints+1)+rangeIm(1);
        c=x+y*1i;
        c1=fnc(c);
        patch(real(c1),imag(c1),cs(i,j,:));
    %    set(h,'Edgecolor','none')
    end

end



