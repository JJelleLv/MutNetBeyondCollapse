function [m] = i_table(A, title)
if ~exist('uitable','builtin')
    varcopy(A);
    m = A;
    disp('a table is copied to the clipboard');
else
    if nargin<2
        title='Table';
    end
    if nargin<1
        A=rand(30);
    end

    hfig = figure('MenuBar','none', ...
        'Color',[0.914 0.914 0.914], ...
        'NumberTitle','off', ...
        'Name',title,...
        'CreateFcn',@(h,evnt)movegui(h, 'center'));
    uitable(hfig,'Units','normalized','Position',...
        [0 0.1 1 0.9], 'Data', A, ...
        'ColumnEditable', true, ...
        'Tag','Table');
    uicontrol('Parent',hfig, ...
        'Units','normalized', ...
        'Callback',@c_cancel, ...
        'ListboxTop',0, ...
        'Position',[0.4,0.02,0.09,0.06], ...
        'String','Cancel', ...
        'Tag','Cancelbutton');
    uicontrol('Parent',hfig, ...
        'Units','normalized', ...
        'Callback',@c_ok, ...
        'ListboxTop',0, ...
        'Position',[0.5,0.02,0.09,0.06], ...
        'String','OK', ...
        'Tag','OKbutton');
    uicontrol('Parent',hfig, ...
        'Units','normalized', ...
        'Callback',@c_clipboard, ...
        'ListboxTop',0, ...
        'Position',[0.6,0.02,0.09,0.06], ...
        'String','Clipboard', ...
        'Tag','Clipboardbutton');
    uicontrol('Parent',hfig, ...
        'Units','normalized', ...
        'Callback',@c_paste, ...
        'ListboxTop',0, ...
        'Position',[0.7,0.02,0.09,0.06], ...
        'String','Paste', ...
        'Tag','Pastebutton');
    uicontrol('Parent',hfig, ...
        'Units','normalized', ...
        'Callback',@c_addcol, ...
        'ListboxTop',0, ...
        'Position',[0.8,0.02,0.09,0.06], ...
        'String','Add column', ...
        'Tag','AddColbutton');
    uicontrol('Parent',hfig, ...
        'Units','normalized', ...
        'Callback',@c_addrow, ...
        'ListboxTop',0, ...
        'Position',[0.9,0.02,0.09,0.06], ...
        'String','Add row', ...
        'Tag','AddRowbutton');
    if nargout>0
        set(hfig,'userdata','uiwait');
        uiwait;
        ud = get(hfig, 'userdata');
        if ~isempty(ud) && strcmp(ud, 'OK')
            m=get(findobj(hfig,'Tag','Table'),'data');
        else
            m = [];
        end

        close(hfig);
    end

end

function c_addcol(hobject,~)
hfig=getparentfig(hobject);
dat=get(findobj(hfig,'Tag','Table'),'data');
if iscell(dat)
    dat{1,end+1}=[];
else
    dat(1,end+1)=0;
end

set(findobj(hfig,'Tag','Table'),'data',dat)
function c_addrow(hobject,~)
hfig=getparentfig(hobject);
dat=get(findobj(hfig,'Tag','Table'),'data');
if iscell(dat)
    dat{end+1,1}=[];
else
    dat(end+1,1)=0;
end

set(findobj(hfig,'Tag','Table'),'data',dat)
function c_clipboard(hobject,~)
hfig=getparentfig(hobject);
dat=get(findobj(hfig,'Tag','Table'),'data');
varcopy(dat);
disp('Copied table to clipboard')
function c_paste(hobject,~)
hfig=getparentfig(hobject);
set(findobj(hfig,'Tag','Table'),'data',varpaste);
function c_ok(hobject,~)
hfig=getparentfig(hobject);
ud=get(hfig,'userdata');
set(hfig,'userdata','OK');
if strcmp(ud,'uiwait')
    uiresume;
else
    close;
end

function c_cancel(hobject,~)
hfig=getparentfig(hobject);
ud=get(hfig,'userdata');
set(hfig,'userdata','cancel');
if strcmp(ud,'uiwait')
    uiresume;
else
    close;
end


