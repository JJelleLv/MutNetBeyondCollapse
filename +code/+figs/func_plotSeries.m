function func_plotSeries(noiseseries_MEAN,noiseseries_LOW,noiseseries_HIGH,Xaxis,windowSizeFrac,CritRange_Mmin,CritRange_Mmax,CritRange_FOUND,COLs)
%%(DATA,NRSTEPS,LOWER,HIGHER,NRDATA,COLs)

hold on

%% set resolution
set(gcf,'units','pixel','position',[500,500,800,400],'papersize',[800,400],'Color',[1 1 1]);

%% set colors
%set(0,'DefaultAxesColorOrder',COLs)

%% data length and number
DATAlength=length(noiseseries_MEAN(:,1));
NRdata=length(noiseseries_MEAN(1,:));

%% set axis
Ylim(1,1)=0;
Ylim(1,2)=(ceil((1.1.*(max(max(noiseseries_HIGH))))*2))/2;
Xlim=[-0.01, 1.01];
set(gca,'Xlim',Xlim,'Ylim',Ylim);

%% plot critical region
if CritRange_FOUND==1
    patch([CritRange_Mmin, CritRange_Mmin, CritRange_Mmax, CritRange_Mmax],[0, Ylim(1,2), Ylim(1,2), 0],[0.85 0.85 0.85],'FaceAlpha',1,'EdgeColor','none')
end

%% make noise useing patch
X_patch([1:DATAlength],:)=Xaxis;
X_patch([DATAlength+1:2*DATAlength],:)=flipud(Xaxis);
for dataNR=1:NRdata
    Y_patch([1:DATAlength],:)=noiseseries_LOW(:,dataNR);
    Y_patch([DATAlength+1:2*DATAlength],:)=flipud(noiseseries_HIGH(:,dataNR));
    patch(X_patch,Y_patch,COLs(dataNR,:),'EdgeColor','none','FaceAlpha',0.3)
    
    %% plot data
    plot(Xaxis,noiseseries_MEAN(:,dataNR),'Color',COLs(dataNR,:),'LineWidth',2.5);
end

%% plot windowsize
deltaY=(Ylim(1,2)-Ylim(1,1));
windowSizeFrac_Ypos=Ylim(1,1)+0.04.*deltaY;
line([0 windowSizeFrac],[windowSizeFrac_Ypos windowSizeFrac_Ypos],'Color',[0 0 0],'LineWidth',1.5)
line([0 0],[windowSizeFrac_Ypos-0.01*deltaY windowSizeFrac_Ypos+0.01*deltaY],'Color',[0 0 0],'LineWidth',1.5)
line([windowSizeFrac windowSizeFrac],[windowSizeFrac_Ypos-0.01*deltaY windowSizeFrac_Ypos+0.01*deltaY],'Color',[0 0 0],'LineWidth',1.5)

%% set x- y-ticks
%line([0.781,0.781],[0 6.5])
%set(gca,'Xtick',[0:0.2:1]);
set(gca,'Xtick',[]);
Yticks=[0:0.5:max(Ylim)];
set(gca,'Ytick',Yticks);

%% set labels
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35])
hget=get(gca,'title'); set(hget,'FontSize',15,'Color',[0.25 0.25 0.25]); 
hget=get(gca,'xlabel'); set(hget,'FontSize',15); set(hget, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0],'Color',[0.25 0.25 0.25]);
hget=get(gca,'ylabel'); set(hget,'FontSize',15); set(hget, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'Color',[0.25 0.25 0.25]);

%axis square
set(gca,'Box','off', 'TickDir', 'out');
axesPosition = get(gca,'Position');          %# Get the current axes position
axes('Position',axesPosition, 'Color','none','YAxisLocation','left','Ytick',[],'XAxisLocation','bottom','XTick',[],'Box','off');
axes('Position',axesPosition, 'Color','none','YAxisLocation','right','Ytick',[],'XAxisLocation','top','XTick',[],'Box','off');
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35])
%axis square

%% prevents white lines from being saved as black
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'auto')

hold off