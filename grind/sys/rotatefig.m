function rotatefig(hax,aviname)
if nargin<1
    hax=gca;
end
if nargin<2
    answer=inputdlg({'Enter filename for rotate function (*.avi):'},'rotate figure',1,{'rotatefig.avi'});
    aviname=answer{1};
end
hfig=get(gca,'Parent');
g_v=get(hax,'view');
winsize=get(hfig,'position');
oldset=get(hax);
set(hax,'YLimMode', 'manual');
set(hax,'ZLimMode', 'manual');
set(hax,'XLimMode', 'manual');
winsize(1:2)=0;
%g_F(1:365)=struct( 'cdata',[],'colormap',[]);
%axis tight;
set(hax,'NextPlot','replacechildren')
if ~isoctave&&verLessThan('matlab','8.6')
    set(hfig,'renderer','ZBuffer'); %in fact this is solving a window bug
end
k=1;
%cm=get(hfig,'colormap');
for g_i=0:1:364
    set(hax,'view',[g_v(1)+g_i,g_v(2)]);
   % pause(0.1)
    drawnow;
    g_F(k) = getframe(hfig,winsize);
    k=k+1;
end
g_F=g_F(1:k-1);
aviobj = VideoWriter(aviname);
open(aviobj);
for g_i=1:length(g_F)
    writeVideo(aviobj,g_F(g_i));
end
close(aviobj);
%set all changed settings back to the original
set(hax,'XLimMode', oldset.XLimMode);
set(hax,'YLimMode', oldset.YLimMode);
set(hax,'ZLimMode', oldset.ZLimMode);
set(hax,'NextPlot',oldset.NextPlot)
set(hax,'View',g_v);
fprintf('written to %s\n',aviname);

