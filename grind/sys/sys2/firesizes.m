%firesizes - special function for the Bak Forest fire model
% Evaluate the sizes of fires after a run.
%
function [firesizs,timestart,P] = firesizes(Fires, siz)
global g_grind g_t t;
N0 = i_initvar;
if i_settingschanged(N0, g_grind.ndays)
    disp('running');
   i_ru(t, g_grind.ndays, N0, 1);
end

if nargin < 1
   Fires=outfun(g_grind.statevars.vectnames{1}) == 2;
   siz = [g_grind.statevars.dims{1}.dim1, g_grind.statevars.dims{1}.dim2];
end

prevFNo = zeros(siz);
FireNos=zeros(size(Fires));
timestart1=[];
fNo = 1;
mask=ones(3);
mask(2,2)=0;
i_waitbar(0, size(Fires,1), 'firesizes', 'Analyzing',0.5)        
for i = 1:size(Fires, 1)
   i_waitbar(1);
   u=unique(prevFNo);
   prevFNo=roundexp(prevFNo,siz+6); %expand the FireNumbers for periodic boundaries
   %expand the previous fire numbers
   expFNo=zeros(size(prevFNo));
   for j=1:length(u)
       if u(j)~=0
         aa=conv2(double(prevFNo==u(j)),mask,'same')>0;
  %       aa=roundconv2(double(prevFNo==u(j)),mask)>0;
         expFNo(aa)=u(j); % if fires are overlapping the largest fireno (=youngest) is taken
       end

   end

   expFNo=expFNo(3+1:3+siz(1),3+1:3+siz(2));
   F = reshape(Fires(i, :), siz);
   FNo = zeros(siz);
   FNo(F>0)=expFNo(F>0);
   newfires=find(F>0&expFNo==0);
   FNo(newfires)=fNo:fNo+length(newfires)-1;
   timestart1(fNo:fNo+length(newfires))=g_t(i);
   fNo=fNo+length(newfires);
   FireNos(i,:)=FNo(:);
   prevFNo=FNo;
end

timestart1=transpose(timestart1);
firesizs1=histc(FireNos(:),1:fNo);
hfig=figure;
i_waitbar([]);
ndx=(timestart1>g_t(1))&(firesizs1>0);
firesizs1=firesizs1(ndx);
timestart1=timestart1(ndx);
plot(timestart1,firesizs1,'o');
i_plotdefaults(hfig);
xlabel('time start of fire');
ylabel('size of fire');
set(gca,'Yscale','log');
figure
P=plotpowerlaw(firesizs1,20,'fire size');
if nargout>0
    firesizs=firesizs1;
    timestart=timestart1;
end
function P=plotpowerlaw(x,nbins,nam)
if nargin<2
    nbins=20;
end

if nargin<3
    nam='x';
end

minx=min(x);
bins=transpose(logspace(log10(minx),log10(max(x)+minx),nbins))-0.5*minx;
diffbins=diff(bins);
freq=histc(x,bins);
freq=freq(1:end-1)./diffbins;
bincenter=transpose(logspace(log10(minx),log10(max(x)),nbins));
bincenter=bincenter(1:end-1);
loglog(bincenter,freq,'o');
hold on;
logbins=log10(bincenter);
freq=freq(freq>0);
logbins=logbins(freq>0);
P=polyfit(logbins,log10(freq),1);
v=polyval(P,logbins);
loglog(10.^(logbins),10.^(v),'k:');
xlim([min(x),max(x)])
htitle=title(sprintf('slope = %g',P(1)));
set(htitle,'fontweight','normal');
ylabel(sprintf('P(%s)',nam),'FontSize',16);
xlabel(nam,'FontSize',16);
