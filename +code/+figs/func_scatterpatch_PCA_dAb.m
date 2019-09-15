function [Xlimit,Ylimit,AX_Xrange,AX_Yrange]=func_scatterpatch_PCA_dAb(PC1_ALL_vect,DELTAab_tpp,NRs_collapsed,regParOrigin,COL_all)

Xdata=PC1_ALL_vect;
Ydata=DELTAab_tpp;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make scatterplot %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on

%% set resolution
set(gcf,'units','pixel','position',[50,50,400,400],'papersize',[700,700],'Color',[1 1 1]);

%% determine axis limits
Xlim=max(abs(Xdata));
Ylim=regParOrigin.*Xlim;
if (1.05*Ylim)<=max(abs(Ydata))
   Ylim=max(abs(Ydata));
end

AX_Xrange=2.2*Xlim;
AX_Yrange=2.2.*Ylim;
AX_ratio=AX_Yrange./AX_Xrange;

Xlimit=[-1.1.*Xlim,1.1.*Xlim];
Ylimit=[-1.1.*Ylim,1.1.*Ylim];

%% plot lines through origin
line(Xlimit,[0 0],'Color',[0.35 0.35 0.35],'LineWidth',1)
line([0 0],Ylimit,'Color',[0.35 0.35 0.35],'LineWidth',1)

%% plot regression line through origin
line(Xlimit,Xlimit.*regParOrigin,'Color',[0.5 0.5 0.5], 'LineWidth',2,'LineStyle','--')

%% make scatter plot with patch
C_angle=0:359;
C_size=0.018;
C_xcoor=cosd(C_angle).*C_size.*AX_Xrange;
C_ycoor=sind(C_angle).*C_size.*AX_Yrange;

%% cross
Crs_length=1.1.*AX_Xrange;%1.2.*AX_Xrange;
Crs_width=0.55.*AX_Xrange; %0.6.*AX_Xrange;
Crs_xcoor=[-Crs_width; -Crs_length; -Crs_length+Crs_width; 0; Crs_length-Crs_width; Crs_length; Crs_width; Crs_length; Crs_length-Crs_width; 0; -Crs_length+Crs_width; -Crs_length].*C_size;
Crs_ycoor=[0; Crs_length-Crs_width; Crs_length; Crs_width; Crs_length; Crs_length-Crs_width; 0; -Crs_length+Crs_width; -Crs_length; -Crs_width; -Crs_length; -Crs_length+Crs_width].*C_size.*AX_ratio;

crossLineLength=0.016;
dCrossLineX=crossLineLength.*AX_Xrange;
dCrossLineY=crossLineLength.*AX_Yrange;

for DATANR=1:length(Xdata(:,1));
    
    %% scatterpoint - if not collapsed
    if (sum(NRs_collapsed==DATANR))==0
        px=Xdata(DATANR,1);
        py=Ydata(DATANR,1);
        COL=COL_all(DATANR,:);
        patch(px+C_xcoor,py+C_ycoor,COL,'FaceAlpha',0.9,'EdgeColor','none')
    end
    
    %% cross - if collapsed
    if (sum(NRs_collapsed==DATANR))==1
        px=Xdata(DATANR,1);
        py=Ydata(DATANR,1);
        COL=COL_all(DATANR,:);
        patch(px+Crs_xcoor,py+Crs_ycoor,COL,'FaceAlpha',0.9,'EdgeColor','none')
    end
    
end

%% set axis limits
set(gca,'Xlim',Xlimit,'Ylim',Ylimit,'Ytick',[-10:1:10],'Xtick',[-1:0.1:1])

%% set labels
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35])
hget=get(gca,'title'); set(hget,'FontSize',15,'Color',[0.25 0.25 0.25]); 
hget=get(gca,'xlabel'); set(hget,'FontSize',15); set(hget, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0],'Color',[0.25 0.25 0.25]);
hget=get(gca,'ylabel'); set(hget,'FontSize',15); set(hget, 'Units', 'Normalized', 'Position', [-0.075, 0.5, 0],'Color',[0.25 0.25 0.25]);

%% make square and other stuff...
axis square
set(gca,'Box','off', 'TickDir', 'out');
axesPosition = get(gca,'Position');          %# Get the current axes position
axes('Position',axesPosition, 'Color','none','YAxisLocation','right','Ytick',[],'XAxisLocation','top','XTick',[],'Box','off');
set(gca,'FontSize',15,'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35])
axis square

%% prevents white lines from being saved as black                        
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'auto', 'Color',[1 1 1])

hold off
