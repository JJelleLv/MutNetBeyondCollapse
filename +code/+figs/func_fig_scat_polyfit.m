function func_fig_scat_polyfit(scatDATA,Xaxis,CritRange_Mmin,CritRange_Mmax,CritRange_FOUND,COLs)

NRDATA=length(scatDATA(:,1));
NRsteps=length(scatDATA(1,:));

hold on

%% set resolution
set(gcf,'units','pixel','position',[50,50,800,400],'papersize',[800,400],'Color',[1 1 1]);

%% determine Ylim and Xlim
minScatDATA=min(min(scatDATA));
maxScatDATA=max(max(scatDATA));
Ylim=[(floor(minScatDATA*20))/20 (ceil(maxScatDATA*20))/20];
Xlim=[-0.01, 1.01];
set(gca,'Xlim',Xlim,'Ylim',Ylim);

%% set x- y-ticks
%set(gca,'Xtick',[0:0.2:1]);
set(gca,'Xtick',[]);
Yticks=[(round(minScatDATA*10))/10:0.1:maxScatDATA];
if (maxScatDATA-minScatDATA)>2 && (maxScatDATA-minScatDATA)<=20
    Yticks=[(round(minScatDATA)):1:maxScatDATA];
elseif (maxScatDATA-minScatDATA)>20
    Yticks=[(round(minScatDATA/10))*10:10:maxScatDATA];
end
set(gca,'Ytick',Yticks);

%% plot critical region
if CritRange_FOUND==1
    patch([CritRange_Mmin, CritRange_Mmin, CritRange_Mmax, CritRange_Mmax],[Ylim(1,1), Ylim(1,2), Ylim(1,2), Ylim(1,1)],[0.85 0.85 0.85],'FaceAlpha',1,'EdgeColor','none')
end

%% put lines in y direction
NRYticks=length(Yticks);
for YtickNR=1:NRYticks
    %line([0 1],[Yticks(1,YtickNR) Yticks(1,YtickNR)],'Color',[0.5 0.5 0.5],'LineWidth',0.5);
end
line([-0.01 1.01],[0 0],'Color',[0 0 0],'LineWidth',2.5)

%% polyfit
% for SpecNR=1:NRDATA
%     maxrange=max(find(isnan(scatDATA(SpecNR,[1:NRsteps]))==0));
%     COEFplf=polyfit(Xaxis([1:maxrange],1)',scatDATA(SpecNR,[1:maxrange]),6);
%     Ypoly=polyval(COEFplf,Xaxis([1:maxrange],1));
%     plot(Xaxis([1:maxrange],1),Ypoly,'Color',COLs(SpecNR,:),'LineWidth',2.5);
% end

%% take mean over range
TRNDrange=3;
TRND_BEGIN=[1:1:(NRsteps-TRNDrange)];
TRND_END=[TRNDrange:1:NRsteps];
NRTRNDdata=length(TRND_BEGIN);
TRNDXaxis=nan(NRTRNDdata,1);
TRNDscatDATA=nan(NRDATA,NRTRNDdata);
for TRNDdataNR=1:NRTRNDdata
    TRND_BEGIN(1,TRNDdataNR);
    TRND_END(1,TRNDdataNR);
    TRNDXaxis(TRNDdataNR,1)=mean(Xaxis([TRND_BEGIN(1,TRNDdataNR):TRND_END(1,TRNDdataNR)],1));
    TRNDscatDATA(:,TRNDdataNR)=mean(scatDATA(:,[TRND_BEGIN(1,TRNDdataNR):TRND_END(1,TRNDdataNR)]),2);
end

%% plot trend
for SpecNR=1:NRDATA
    plot(TRNDXaxis,TRNDscatDATA(SpecNR,:),'Color',COLs(SpecNR,:),'LineWidth',2.5);
end

%% make scatter plot
%for SpecNR=1:NRDATA
%    scatter(Xaxis,scatDATA(SpecNR,:),40,[1 1 1],'filled');
%    scatter(Xaxis,scatDATA(SpecNR,:),20,COLs(SpecNR,:),'filled');
%end

%% make scatter plot with patch
AX_Xrange=(Xlim(1,2)-Xlim(1,1))/2;
AX_Yrange=Ylim(1,2)-Ylim(1,1);
C_angle=0:359;

C_size=0.0085;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;
for SpecNR=1:NRDATA
    for StepNR=1:NRsteps
        px=Xaxis(StepNR,1);
        py=scatDATA(SpecNR,StepNR);
        COL=COLs(SpecNR,:);
        patch(px+C_xcoor,py+C_ycoor,[1 1 1],'FaceAlpha',1,'EdgeColor','none')
    end
end

C_size=0.007;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;
for SpecNR=1:NRDATA
    for StepNR=1:NRsteps
        px=Xaxis(StepNR,1);
        py=scatDATA(SpecNR,StepNR);
        COL=COLs(SpecNR,:);
        patch(px+C_xcoor,py+C_ycoor,COL,'FaceAlpha',0.75,'EdgeColor','none')
    end
end
        
%% put arrows
%annotation('textarrow',[0.925 0.925],[0.9-0.6.*(abs(Ylim(1,2))./sum(abs(Ylim))) 0.9])
%annotation('textarrow',[0.925 0.925],[0.175+0.6.*(abs(Ylim(1,1))./sum(abs(Ylim))) 0.175])

% %% plot windowsize
% deltaY=(Ylim(1,2)-Ylim(1,1));
% windowSizeFrac_Ypos=Ylim(1,1)+0.9.*deltaY;
% line([0 windowSizeFrac],[windowSizeFrac_Ypos windowSizeFrac_Ypos],'Color',[0 0 0],'LineWidth',1.5)

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
