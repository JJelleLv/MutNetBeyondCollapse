function result = i_myerrordlg(ErrorString, DlgName)
if nargin==0
   ErrorString='Error';
end

if nargin<2
   DlgName='Error';
end

BoxHandle = errordlg(ErrorString, DlgName, 'modal');
OKHandle=findobj( BoxHandle,'Tag','OKButton');
OKPos = get(OKHandle, 'Position');
OKPos(3)=floor(OKPos(3)*1.5);
OKPos(1) = OKPos(1) - floor(0.6 * OKPos(3));
set(OKHandle, 'Position', OKPos);
set(OKHandle,'CallBack',@(h,eventdata)set(get(h,'parent'),'Userdata',1),...
   'String','Edit model');
OKPos(1) = OKPos(1) + floor(1.2 * OKPos(3));
Font.FontUnits = 'points';
Font.FontSize = get(0, 'FactoryUIControlFontSize');
uicontrol(BoxHandle, Font, ...
   'Style', 'pushbutton', ...
   'Units', 'points', ...
   'Position', OKPos, ...
   'CallBack', @(h,eventdata)set(get(h,'parent'),'Userdata',2), ...
   'String', 'Ignore error', ...
   'HorizontalAlignment', 'center', ...
   'Tag', 'IgnoreButton');
waitfor(BoxHandle , 'userdata')
if ishandle(BoxHandle )
   result = get(BoxHandle , 'userdata');
   delete(BoxHandle );
else
   result = 1;
end


