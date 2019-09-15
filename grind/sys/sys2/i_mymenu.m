function i_mymenu(ud)
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

%  i_mymenu(mess,struct('ndx',ndx,'N0',N0,'ix',ix,'iy',iy,'p1',p1,'p2',p2,'oldNO',oldNO));  %dialog box for selecting an action

if nargin<1
    ud=get(gca,'userdata');
end
if length(ud.initpt)==2
    mess=sprintf('Mouse cursor at:\n%s%s%s',dispval(ud.meta.xname,ud.initpt(1)),dispval(ud.meta.yname,ud.initpt(2)));
else
    mess=sprintf('Mouse cursor at:\n%s%s%s',dispval(ud.meta.xname,ud.initpt(1)),dispval(ud.meta.yname,ud.initpt(2)),dispval(ud.meta.zname,ud.initpt(3)));
end
%load i_mymenu
hfig=i_figno('dialog');
if ishandle(hfig)
   close(hfig);
end

hfig = figure(hfig);
set(hfig,'userdata',ud);
set(hfig, 'KeypressFcn',@(hobj,ev)i_callb('keypressed',hobj));
set(hfig,'Units','points', ...
	'Color',[0.914 0.914 0.914], ...
	'Colormap',[], ...
	'MenuBar','none', ...
	'Name','State variables', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[200 170 160 240] ,...
	'Resize','off', ...
	'Tag','Fig1', ...
   'ToolBar','none',...
    'CreateFcn',@(h,evnt)movegui(h, 'center'));

uicontrol('Parent',hfig, ...
	'Units','points', ...
	'BackgroundColor',[0.914 0.914 0.914], ...
	'ListboxTop',0, ...
	'Position',[10 170 140 60], ...
	'String',mess, ...
	'Style','text', ...
	'Tag','StaticText1');
uicontrol('Parent',hfig, ...
	'Units','points', ...
	'Callback',@(hObject,callbackdata)runclick(hObject,callbackdata,1), ...
	'ListboxTop',0, ...
	'Position',[10 160 140 25], ...
	'String','Run forward', ...
    'TooltipString','Set initial conditions to the mouse cursor and simulate from that point (ru)', ...
	'Tag','Pushbutton6');
uicontrol('Parent',hfig, ...
	'Units','points', ...
	'Callback',@(hObject,callbackdata)runclick(hObject,callbackdata,2), ...
	'ListboxTop',0, ...
	'Position',[10 130 140 25], ...
	'String','Find nearest equilibrium', ...
  	'TooltipString','Set initial conditions to the mouse cursor and try to find equilibrium (findeq)', ...
	'Tag','Pushbutton5');
uicontrol('Parent',hfig, ...
	'Units','points', ...
	'Callback',@(hObject,callbackdata)runclick(hObject,callbackdata,3), ...
	'ListboxTop',0, ...
	'Position',[10 100 140 25], ...
	'String','Set initial values only', ...
  	'TooltipString','Set initial conditions to the values at the mouse cursor', ...
	'Tag','Pushbutton4');
uicontrol('Parent',hfig, ...
	'Units','points', ...
	'Callback',@(hObject,callbackdata)runclick(hObject,callbackdata,4), ...
	'ListboxTop',0, ...
	'Position',[10 70 140 25], ...
	'String','Enable plot editing', ...
  	'TooltipString','Neglect mouse click and enable plot editing', ...
	'Tag','Pushbutton3');
uicontrol('Parent',hfig, ...
	'Units','points', ...
	'Callback',@(hObject,callbackdata)runclick(hObject,callbackdata,5), ...
	'ListboxTop',0, ...
	'Position',[10 40 140 25], ...
	'String','Cancel', ...
   	'TooltipString','Neglect mouse click', ...
	'Tag','Pushbutton2');
uicontrol('Parent',hfig, ...
	'Units','points', ...
	'Callback','commands phasclick', ...
	'ListboxTop',0, ...
	'Position',[10 10 140 25], ...
	'String','Help', ...
   	'TooltipString','Open help text', ...
	'Tag','Pushbutton1');
%set(hfig, 'Visible', 'on' );

%------------------------------------------------------------------------
% Wait for choice to be made (i.e UserData must be assigned)...
%------------------------------------------------------------------------
% waitfor(hfig,'userdata')
% if ishandle(hfig)
%    result = get(hfig,'userdata');
%    delete(hfig);
% else
%    result=5;
% end


function runclick(hObject,~,no) 
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
close(hfig);
drawnow;
i_callb('mymenucallback',gcf,no,ud);

function s=dispval(avar,aval)
if ~isempty(avar)&&~any(avar=='''')&&~strcontains(avar,'(t+')
    s=sprintf('%s = %g\n',avar,aval);
else
    s='';
end




