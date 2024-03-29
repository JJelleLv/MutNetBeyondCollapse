function [vidObj,frames1]=i_videodlg(frames)
if nargin==0
   frames=[];
end

profs=VideoWriter.getProfiles();
profilenames=cell(size(profs));
for i=1:length(profs)
    Fs=fields(profs(i));
    profilenames{i}=profs(i).Name;
    prof.Name=profs(i).Name;
    prof.FileExtension=profs(i).FileExtensions{1};
    if any(strcmp(Fs,'Quality'))
        prof.Quality=profs(i).Quality;
    else
        prof.Quality=nan;
    end

    if any(strcmp(Fs,'FrameRate'))
        prof.FrameRate=profs(i).FrameRate;
    else
        prof.FrameRate=nan;
    end

    if any(strcmp(Fs,'CompressionRatio'))
        prof.CompressionRatio=profs(i).CompressionRatio;
    else
        prof.CompressionRatio=nan;
    end

    if any(strcmp(Fs,'LosslessCompression'))
        prof.LosslessCompression=profs(i).LosslessCompression;
    else
        prof.LosslessCompression=nan;
    end

    ud.profs(i)=prof;
end

ud.No=find(strcmp(profilenames,'Motion JPEG AVI'));
if isempty(ud.No)
    ud.No=1;
end

ud.frames=frames;
ud.Ok=0;
ud.vidObj=[];
ud.curr=ud.profs(ud.No);

hfig = figure('Units','characters',...
    'PaperUnits',get(0,'defaultfigurePaperUnits'),...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
    'IntegerHandle','off',...
    'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
    'MenuBar','none',...
    'Name','Create video file',...
    'NumberTitle','off',...
    'PaperPosition',get(0,'defaultfigurePaperPosition'),...
    'PaperSize',get(0,'defaultfigurePaperSize'),...
    'PaperType',get(0,'defaultfigurePaperType'),...
    'Position',[103.8 35.7692307692308 94 25.6923076923077],...
    'Resize','off',...
    'HandleVisibility','callback',...
    'UserData',[],...
    'Tag','figure1',...
    'Visible','on' );

uicontrol('Parent',hfig,...
    'Units','characters',...
    'Callback',@OK_callback,...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Position',[36.6 2.69230769230769 13.8 1.69230769230769],...
    'String','OK',...
    'Tag','OKbutton' );

uicontrol('Parent',hfig,...
    'Units','characters',...
    'Callback',@Cancel_callback,...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Position',[53.2 2.69230769230769 13.8 1.69230769230769],...
    'String','Cancel',...
    'Tag','Cancelbutton');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'Callback',@profilepopup_callback,...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'Position',[31 22.6923076923077 40.2 1.53846153846154],...
    'String',profilenames,...
    'Style','popupmenu',...
    'Value',ud.No,...
    'Tag','Profilepopup');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[3.4 22.6923076923077 17.8 1.07692307692308],...
    'String','Video profile:',...
    'Style','text',...
    'Tag','text1');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[3.4 19.8461538461538 17.8 1.07692307692308],...
    'String','File Name:',...
    'Style','text',...
    'Tag','text2' );

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[31 19.7692307692308 36.2 1.69230769230769],...
    'String',{  'video.avi' },...
    'Style','edit',...
    'Tag','FileNameEdit');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'Callback' ,@fileedit_callback,...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Position',[67.2 19.7692307692308 4.4 1.69230769230769],...
    'String','...',...
    'Tag','FileEditButton' );

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[3.4 14.1538461538462 20 1.07692307692308],...
    'String','Compression ratio',...
    'Style','text',...
    'Tag','text3');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[3.4 17 17.8 1.07692307692308],...
    'String','Frame rate:',...
    'Style','text',...
    'Tag','text4');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[3.4 11.3076923076923 17.8 1.07692307692308],...
    'String','Quality:',...
    'Style','text',...
    'Tag','text5' );

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[3.4 8.46153846153846 29.8 1.07692307692308],...
    'String','Lossless compression',...
    'Style','text',...
    'Tag','text6');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[31 16.8461538461538 20.2 1.69230769230769],...
    'String','',...
    'Style','edit',...
    'Tag','FrameRateEdit');


uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[31 13.9230769230769 20.2 1.69230769230769],...
    'String','',...
    'Style','edit',...
    'Tag','CompressionRatioEdit');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'Callback',@saveframe_callback,...
    'FontSize',10.6666666666667,...
    'Position',[20 2.69230769230769 13.8 1.69230769230769],...
    'String','Save Frames',...
    'Tag','SaveFramesbutton');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'Callback',@loadframe_callback,...
    'FontSize',10.6666666666667,...
    'Position',[3.4 2.69230769230769 13.8 1.69230769230769],...
    'String','Load Frames',...
    'Tag','LoadFramesbutton');

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'Callback',@help_callback,...
    'FontSize',10.6666666666667,...
    'Position',[70.2 2.69230769230769 13.8 1.69230769230769],...
    'String','Help',...
    'Tag','pushbutton7' );

uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[31 11 20.2 1.69230769230769],...
    'String','',...
    'Style','edit',...
    'Tag','QualityEdit');
uicontrol('Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'Position',[31 8.23076923076923 19.8 1.53846153846154],...
    'String',{  'Yes'; 'No'; 'N/A'},...
    'Style','popupmenu',...
    'Value',1,...
    'Tag','Losslesspopup');
set(hfig,'userdata',ud);
movegui(hfig,'center');
updateprofile(hfig,ud.curr);
uiwait(hfig);
if ishandle(hfig)
    ud=get(hfig,'userdata');
    frames1=ud.frames;
    vidObj=ud.vidObj;
    close(hfig);
end

function loadframe_callback(hObject, ~)
hfig=get(hObject,'parent');
ud=get(hfig,'userdata');
 filter={'*.mat','MAT-files (*.mat)'; '*.*',  'All Files (*.*)'};
[filename,pathname]= uigetfile(filter,'Select file name to load frames');
load([pathname filename],'frames')
ud.frames=frames;
set(hfig,'userdata',ud);

function saveframe_callback(hObject, ~)
hfig=get(hObject,'parent');
ud=get(hfig,'userdata');
h=findobj(hfig,'tag','FileNameEdit');
filename=get(h,'string');
if iscell(filename)
    filename=filename{:};
end

[pathname,filename]=fileparts(char(filename));
if isempty(pathname)
    filename=[filename '.mat'];
else
    filename=[pathname filesep filename '.mat'];
end

 filter={'*.mat','MAT-files (*.mat)'; '*.*',  'All Files (*.*)'};
[filename,pathname]= uiputfile(filter,'Select file name to write frames',filename);
frames=ud.frames; 
save([pathname filename],'frames')


function help_callback(~, ~)
doc('videoWriter')
function OK_callback(hObject, ~)
hfig=get(hObject,'parent');
ud=get(hfig,'userdata');
ud.curr=getprofile(hfig, ud.curr);
ud.profs(ud.No)=ud.curr;
if ~isempty(ud.frames)
    h=findobj(hfig,'Tag','FileNameEdit');
    filename=get(h,'string');
    if isempty(filename)
        filename='video.avi';
    end

    ud.vidObj=VideoWriter(char(filename),ud.curr.Name);
    if ~isnan(ud.curr.Quality)
        ud.vidObj.Quality=ud.curr.Quality;
    end

    if ~isnan(ud.curr.FrameRate)
        ud.vidObj.FrameRate=ud.curr.FrameRate;
    end

    if ~isnan(ud.curr.CompressionRatio)
        ud.vidObj.CompressionRatio=ud.curr.CompressionRatio;
    end

    if ~isnan(ud.curr.LosslessCompression)
        ud.vidObj.LosslessCompression=ud.curr.LosslessCompression;
    end

end

set(hfig,'userdata',ud);
uiresume;
function Cancel_callback(hObject, ~)
hfig=get(hObject,'parent');
ud=get(hfig,'userdata');
ud.vidObj=[];
set(hfig,'userdata',ud);
uiresume;
function fileedit_callback(hObject, ~)
hfig=get(hObject,'parent');
h=findobj(hfig,'tag','FileNameEdit');
filename=get(h,'string');
[filename,pathname]=uiputfile(filename);
if ischar(filename)
  set(h,'string',[pathname filename]);
end

function profilepopup_callback(hObject, ~)
hfig=get(hObject,'parent');
v=get(hObject,'value');
ud=get(hfig,'userdata');
ud.curr=getprofile(hfig, ud.curr);
ud.profs(ud.No)=ud.curr;
ud.No=v;
ud.curr=ud.profs(ud.No);
set(hfig,'userdata',ud);
updateprofile(hfig,ud.curr);
h=findobj(hfig,'tag','FileNameEdit');
filename=get(h,'string');
if iscell(filename)
    filename=filename{1};
end

[pathname,filename]=fileparts(char(filename));
if isempty(pathname)
    filename=[filename ud.curr.FileExtension];
else
    filename=[pathname filesep filename ud.curr.FileExtension];
end

set(h,'string',filename);

function updateprofile(hdlg,prof)
h=findobj(hdlg,'Tag','QualityEdit');
if ~isnan(prof.Quality)
    set(h, 'String',sprintf('%g',prof.Quality));
    set(h, 'Enable', 'on');
else
    set(h, 'String', 'N/A');
    set(h, 'Enable', 'off');
end

h=findobj(hdlg,'Tag','FrameRateEdit');
if ~isnan(prof.FrameRate)
    set(h, 'String',sprintf('%g',prof.FrameRate));
    set(h, 'Enable', 'on');
else
    set(h, 'String', 'N/A');
    set(h, 'Enable', 'off');
end

h=findobj(hdlg,'Tag','CompressionRatioEdit');
if ~isnan(prof.CompressionRatio)
    set(h, 'String',sprintf('%g',prof.CompressionRatio));
    set(h, 'Enable', 'on');
else
    set(h, 'String', 'N/A');
    set(h, 'Enable', 'off');
end

h=findobj(hdlg,'Tag','Losslesspopup');
if ~isnan(prof.LosslessCompression)
    if prof.LosslessCompression
        set(h, 'Value',1);
    else
        set(h, 'Value',2);
    end

    set(h, 'Enable', 'on');
else
    set(h, 'Value', 3);
    set(h, 'Enable', 'off');
end

function prof=getprofile(hdlg,prof)
h=findobj(hdlg,'Tag','QualityEdit');
v=str2num(get(h,'string')); 
prof.Quality=v;
h=findobj(hdlg,'Tag','FrameRateEdit');
v=str2num(get(h,'string')); 
prof.FrameRate=v;
h=findobj(hdlg,'Tag','CompressionRatioEdit');
v=str2num(get(h,'string')); 
prof.CompressionRatio=v;
h=findobj(hdlg,'Tag','Losslesspopup');
v=get(h,'value');
prof.LosslessCompression=v==2;
if v==3
    prof.LosslessCompression=NaN;
end



