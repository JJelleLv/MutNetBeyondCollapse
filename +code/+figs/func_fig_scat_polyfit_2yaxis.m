function func_fig_scat_polyfit_2yaxis(scatDATA_1,scatDATA_2,Xaxis_1,Xaxis_2,CritRange_Mmin,CritRange_Mmax,CritRange_FOUND,COLs_1,COLs_2)

NRDATA_1=length(scatDATA_1(:,1));
NRsteps_1=length(scatDATA_1(1,:));
NRDATA_2=length(scatDATA_2(:,1));
NRsteps_2=length(scatDATA_2(1,:));

hold on

%% set resolution
set(gcf,'units','pixel','position',[50,50,800,400],'papersize',[800,400],'Color',[1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lines and other properties %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% determine Ylim and Xlim, Y and X ticks
minscatDATA_1=min(min(scatDATA_1));
maxscatDATA_1=max(max(scatDATA_1));
Xlim_1=[-0.01, 1.01];
Ylim_1=[(floor(minscatDATA_1*20))/20-0.025 (ceil(maxscatDATA_1*20))/20+0.025];
Yticks_1=[0:0.05:Ylim_1(1,2)-0.025];
if (maxscatDATA_1-minscatDATA_1)>2 && (maxscatDATA_1-minscatDATA_1)<=20
    Ylim_1=[(floor(minscatDATA_1)) (ceil(maxscatDATA_1))];
    Yticks_1=[Ylim_1(1,1):1:Ylim_1(1,2)];
elseif (maxscatDATA_1-minscatDATA_1)>20
    Ylim_1=[(floor(minscatDATA_1/10)*10) (ceil(maxscatDATA_1/10)*10)];
    Yticks_1=[Ylim_1(1,1):10:Ylim_1(1,2)];
end
set(gca,'Xtick',[0:0.2:1]);
set(gca,'Ytick',Yticks_1);
set(gca,'Xlim',Xlim_1,'Ylim',Ylim_1);

%% plot critical region
if CritRange_FOUND==1
    patch([CritRange_Mmin, CritRange_Mmin, CritRange_Mmax, CritRange_Mmax],[Ylim_1(1,1), Ylim_1(1,2), Ylim_1(1,2), Ylim_1(1,1)],[0.85 0.85 0.85],'FaceAlpha',1,'EdgeColor','none')
end

%% put lines in y direction
NRYticks_1=length(Yticks_1);
for YtickNR=1:NRYticks_1
    %line([0 1],[Yticks_1(1,YtickNR) Yticks_1(1,YtickNR)],'Color',[0.5 0.5 0.5],'LineWidth',0.5);
end
%line([0 1],[0 0],'Color',[0 0 0],'LineWidth',2.5) %% no thick line on 0

%% formatting
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35], 'TickDir', 'out')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% data for right axis %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create second Y axes on the right.
set(gca,'Box','off')
axes('YAxisLocation', 'Right', 'color', 'none','XaxisLocation','top','Xtick',[])
set(gca,'Box','off', 'TickDir', 'out');
hold on

%% determine Ylim and Xlim, Y and X ticks
minscatDATA_2=min(min(scatDATA_2));
maxscatDATA_2=max(max(scatDATA_2));
Xlim_2=[-0.01, 1.01];
Ylim_2=[(floor(minscatDATA_2*50))/50 (ceil(maxscatDATA_2*50))/50];
Yticks_2=[Ylim_2(1,1):0.01:Ylim_2(1,2)];
if (maxscatDATA_2-minscatDATA_2)>0.14 && (maxscatDATA_2-minscatDATA_2)<=2
    Ylim_2=[(floor(minscatDATA_2*20))/20-0.025 (ceil(maxscatDATA_2*20))/20];
    Yticks_2=[0:0.05:Ylim_2(1,2)];
elseif (maxscatDATA_2-minscatDATA_2)>2 && (maxscatDATA_2-minscatDATA_2)<=10
    Ylim_2=[(floor(minscatDATA_2)) (ceil(maxscatDATA_2))];
    dYtick_2=ceil(Ylim_2(1,2)/5);
    Yticks_2=[Ylim_2(1,1):dYtick_2:Ylim_2(1,2)];
elseif (maxscatDATA_2-minscatDATA_2)>10 && (maxscatDATA_2-minscatDATA_2)<=20
    Ylim_2=[(floor(minscatDATA_2)) (ceil(maxscatDATA_2))];
    dYtick_2=ceil(Ylim_2(1,2)/5);
    Yticks_2=[Ylim_2(1,1):dYtick_2:Ylim_2(1,2)];
elseif (maxscatDATA_2-minscatDATA_2)>20
    Ylim_2=[(floor(minscatDATA_2/10)*10) (ceil(maxscatDATA_2/10)*10)];
    dYtick_2=5;
    Yticks_2=[Ylim_2(1,1):dYtick_2:Ylim_2(1,2)];
end

%% Yticks_2 relative to Yticks_1
%dYtick_2=(Ylim_2(1,2)-Ylim_2(1,1))./(length(Yticks_1)-1);
%Yticks_2=round([Ylim_2(1,1):dYtick_2:Ylim_2(1,2)]);
set(gca,'Ytick',Yticks_2);
set(gca,'Xlim',Xlim_2,'Ylim',Ylim_2);

%% take mean over range
TRNDrange=3;
TRND_BEGIN=[1:1:(NRsteps_2-TRNDrange)];
TRND_END=[TRNDrange:1:NRsteps_2];
NRTRNDdata=length(TRND_BEGIN);
TRNDXaxis_2=nan(NRTRNDdata,1);
TRNDscatDATA_2=nan(NRDATA_2,NRTRNDdata);
for TRNDdataNR=1:NRTRNDdata
    TRND_BEGIN(1,TRNDdataNR);
    TRND_END(1,TRNDdataNR);
    TRNDXaxis_2(TRNDdataNR,1)=mean(Xaxis_2([TRND_BEGIN(1,TRNDdataNR):TRND_END(1,TRNDdataNR)],1));
    TRNDscatDATA_2(:,TRNDdataNR)=mean(scatDATA_2(:,[TRND_BEGIN(1,TRNDdataNR):TRND_END(1,TRNDdataNR)]),2);
end

%% plot trend
for SpecNR=1:NRDATA_2
    plot(TRNDXaxis_2,TRNDscatDATA_2(SpecNR,:),'Color',COLs_2(SpecNR,:),'LineWidth',2.5);
end

%% make scatter plot with patch
AX_Xrange=(Xlim_2(1,2)-Xlim_2(1,1))/2;
AX_Yrange=Ylim_2(1,2)-Ylim_2(1,1);
C_angle=0:359;

C_size=0.0085;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;
for SpecNR=1:NRDATA_2
    for StepNR=1:NRsteps_2
        px=Xaxis_2(StepNR,1);
        py=scatDATA_2(SpecNR,StepNR);
        COL=COLs_2(SpecNR,:);
        patch(px+C_xcoor,py+C_ycoor,[1 1 1],'FaceAlpha',1,'EdgeColor','none')
    end
end

C_size=0.007;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;
for SpecNR=1:NRDATA_2
    for StepNR=1:NRsteps_2
        px=Xaxis_2(StepNR,1);
        py=scatDATA_2(SpecNR,StepNR);
        COL=COLs_2(SpecNR,:);
        patch(px+C_xcoor,py+C_ycoor,COL,'FaceAlpha',0.75,'EdgeColor','none')
    end
end

%% formatting
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35], 'TickDir', 'out')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% data for left axis %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create second Y axes on the right.
set(gca,'Box','off');
axes('YAxisLocation', 'Left', 'color', 'none','XaxisLocation','bottom','Ytick',[])
set(gca,'Box','off', 'TickDir', 'out');
set(gca,'Xlim',Xlim_1,'Ylim',Ylim_1);
set(gca,'Xtick',[0:0.2:1]);
hold on

%% take mean over range
TRNDrange=3;
TRND_BEGIN=[1:1:(NRsteps_1-TRNDrange)];
TRND_END=[TRNDrange:1:NRsteps_1];
NRTRNDdata=length(TRND_BEGIN);
TRNDXaxis_1=nan(NRTRNDdata,1);
TRNDscatDATA_1=nan(NRDATA_1,NRTRNDdata);
for TRNDdataNR=1:NRTRNDdata
    TRND_BEGIN(1,TRNDdataNR);
    TRND_END(1,TRNDdataNR);
    TRNDXaxis_1(TRNDdataNR,1)=mean(Xaxis_1([TRND_BEGIN(1,TRNDdataNR):TRND_END(1,TRNDdataNR)],1));
    TRNDscatDATA_1(:,TRNDdataNR)=mean(scatDATA_1(:,[TRND_BEGIN(1,TRNDdataNR):TRND_END(1,TRNDdataNR)]),2);
end

%% plot trend
for SpecNR=1:NRDATA_1
    plot(TRNDXaxis_1,TRNDscatDATA_1(SpecNR,:),'Color',COLs_1(SpecNR,:),'LineWidth',2.5);
end

%% make scatter plot with patch
AX_Xrange=(Xlim_1(1,2)-Xlim_1(1,1))/2;
AX_Yrange=Ylim_1(1,2)-Ylim_1(1,1);
C_angle=0:359;

C_size=0.0085;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;
for SpecNR=1:NRDATA_1
    for StepNR=1:NRsteps_1
        px=Xaxis_1(StepNR,1);
        py=scatDATA_1(SpecNR,StepNR);
        COL=COLs_1(SpecNR,:);
        patch(px+C_xcoor,py+C_ycoor,[1 1 1],'FaceAlpha',1,'EdgeColor','none')
    end
end

C_size=0.007;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;
for SpecNR=1:NRDATA_1
    for StepNR=1:NRsteps_1
        px=Xaxis_1(StepNR,1);
        py=scatDATA_1(SpecNR,StepNR);
        COL=COLs_1(SpecNR,:);
        patch(px+C_xcoor,py+C_ycoor,COL,'FaceAlpha',0.75,'EdgeColor','none')
    end
end

%% formatting
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35])

%% set labels
hget=get(gca,'title'); set(hget,'FontSize',15,'Color',[0.25 0.25 0.25]); 
hget=get(gca,'xlabel'); set(hget,'FontSize',15); set(hget, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0],'Color',[0.25 0.25 0.25]);
hget=get(gca,'ylabel'); set(hget,'FontSize',15); set(hget, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'Color',[0.25 0.25 0.25]);

%% prevents white lines from being saved as black
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'auto')

hold off
