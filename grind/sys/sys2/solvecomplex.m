function solvecomplex(fn,rangereal,rangecomplex,npoints)
if nargin<2
    rangereal=[-10,10];
end

if nargin<3
    rangecomplex=[-100,100];
end

if nargin<4
    npoints=100;
end

if nargin==0
    fn=@(x)(x+exp(-x*1));
end

[reX,imX] = meshgrid(linspace(rangereal(1),rangereal(2),npoints), linspace(rangecomplex(1),rangecomplex(2),npoints));
x=reX+1i*imX;
y=fn(x);
figure;
contour(reX, imX, real(y), [0, 0],'b:');
hold on
contour(reX, imX, imag(y), [0, 0],'r:');
xlabel('real part of \lambda')
ylabel('complex part of \lambda')
legend('Re(\lambda)==0','Imag(\lambda)==0')

