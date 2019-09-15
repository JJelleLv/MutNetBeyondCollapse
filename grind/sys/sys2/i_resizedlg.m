function newpos=i_resizedlg(hfig,tags,anchors,origsizes,figorigpos)
%if you run with one or no arguments you can create the anchors
%interactively
newpos=[];
if nargin<5
    if nargin>1&&ischar(tags)
        switch tags
            case 'vertdistrib'
                tags=anchors;
                h1=findobj(hfig,'tag',tags{1});
                if isempty(h1)
                    return;
                end
                %the first tag should be the top object
                pos1=get(h1,'position');
                h2=findobj(hfig,'tag',tags{2});
                %the second tag should be the bottom object
                pos2=get(h2,'position');
                newpos=[0 pos1(4)];
                if pos1(4)~=pos2(4)
                    commonheight = (pos1(4) + pos2(4)) / 2;        
                    %move the top of the parameters panel
                    movetop = commonheight - pos2(4);
                    pos1(2) = pos1(2) + movetop;
                    pos1(4) = commonheight;
                    pos2(4) = commonheight;
                    newpos=[movetop commonheight];
                    set(h1,'position',pos1);
                    set(h2,'position',pos2);
                end
           case 'horizdistrib'
                tags=anchors;
                h1=findobj(hfig,'tag',tags{1});
                if isempty(h1)
                    return;
                end
                %the first tag should be the left object
                pos1=get(h1,'position');
                h2=findobj(hfig,'tag',tags{2});
                %the second tag should be the right object
                pos2=get(h2,'position');
                newpos=[0 pos1(3)];
                if pos1(3)~=pos2(3)
                    commonwidth = (pos1(3) + pos2(3)) / 2;        
                    %move the top of the parameters panel
                    moveleft = pos2(3)- commonwidth;
                    pos2(1) = pos2(1) + moveleft;
                    pos1(3) = commonwidth;
                    pos2(3) = commonwidth;
                    newpos=[moveleft commonwidth];
                    set(h1,'position',pos1);
                    set(h2,'position',pos2);
                end
            case 'horizcenter'
                tags=anchors;
                posfig=get(hfig,'position');
                pos1=get(findobj(hfig,'tag',tags{1}),'position');
                if isempty(pos1)
                    return;
                end
                pos2=get(findobj(hfig,'tag',tags{end}),'position');
                width=pos2(1)+pos2(3)-pos1(1);
                leftadapt=(posfig(3)-width)/2-pos1(1);
                for i=1:length(tags)
                    h=findobj(hfig,'tag',tags{i});
                    pos=get(h,'position');
                    pos(1)=pos(1)+leftadapt;
                    set(h,'position',pos);
                end

        end
        return;
    end
                % posmod=newpos{strcmp(data(:,1),'Model')};
% posparlab=newpos{strcmp(data(:,1),'parameterslabel')};
% pospar=newpos{strcmp(data(:,1),'Parameters')};
% 
% if ~isempty(posmod) && (posmod(4) ~= pospar(4))
%     % even distribution between model and parameters
%     %average height
%     commonheight = (posmod(4) + pospar(4)) / 2;
%     
%     %move the top of the parameters panel
%     movetop = commonheight - pospar(4);
%     posmod(2) = posmod(2) + movetop;
%     %move the label too
%     posparlab(2) = posparlab(2) + movetop;
% 
%     pospar(4) = commonheight;
%     posmod(4) = commonheight;
%     set(findobj(hfig,'tag','Model'),'position',posmod);
%     set(findobj(hfig,'tag','Parameters'),'position',pospar);
%     set(findobj(hfig,'tag','parameterslabel'),'position',posparlab);
% end
    if nargin==0
        hfig=gcf;
    end
    if nargin>2
        %you an supply the old anchors if the positions need to be updated
        oldanchors=anchors;
        oldtags=tags;
    else
        oldanchors=[];
    end
    ch=flipud(get(hfig,'children'));
    typ=get(ch,'type');
    %list here the kind of objects that need to be adapted
    if isempty(typ)
        disp('No objects in the figure');
        return;
    end
    ndx=ismember(typ, {'uicontrol','uitable','uipanel','uibuttongroup'});
    ch=ch(ndx);
    set(hfig,'units','points');
    set(ch,'units','points');
    tags=get(ch,'tag');
    utags=unique(tags);
    if length(utags)<length(tags)
        fprintf('%s\n',tags{:})
        error('resize:dlg','We need unique tags')
    end
    figorigpos=get(hfig,'position');
    anchors=repmat({'[1 0 1 0]'},length(tags),1);
    origsizes=get(ch,'pos');
    for i=1:size(origsizes,1)
        origsizes{i}=mat2str(origsizes{i});
    end
    if ~isempty(oldanchors)
        for i=1:size(oldanchors)
            ndx=strcmp(oldtags{i},tags);
            if any(ndx)
                anchors{ndx}=mat2str(oldanchors{ndx});
            end
        end
    end
    resizetagdlg(hfig,[tags,anchors,origsizes],figorigpos);
    return
end
% positionobj(f, 'Helpbutton', [0 1 0 1], [339 78 83 25], [figwidth figheight], p(3:4));
%figorigpos,tags,anchors,origsizes
%function pos = positionobj(hfig, tag, anchors, origsizes, figorigpos)
% adj = left right top bottom
% position = left bottom width height
if length(figorigpos)==4
    figorigpos=figorigpos(3:4);
end
set(hfig,'units','points');
p = get(hfig, 'position');
newfigsize=p(3:4);
newpos=cell(size(tags));
for i=1:length(tags)
    newpos{i}=positionobj(hfig, tags{i}, anchors{i}, origsizes{i}, figorigpos, newfigsize);
end

function pos = positionobj(f, tag, adj, defsize, oldsize, newsize)
%posmod = positionobj(f, 'Model', [1 1 1 1], [11 143 317 116], [figwidth figheight], p(3:4));
% adj = left right top bottom
% position = left bottom width height
h = findobj(f, 'tag', tag);
set(h,'units','points');
if ~isempty(h)
    pos = defsize;
    if adj(2) && ~adj(1)
        pos(1) = newsize(1) - (oldsize(1) - pos(1));
    end
    if adj(2) && adj(1)
        pos(3) = pos(3) - (oldsize(1) - newsize(1));
        if pos(3) < 1
            pos(3) = 1;
        end
    end
    if adj(3) && ~adj(4)
        pos(2) = newsize(2) - (oldsize(2) - pos(2));
    end

    if adj(4) && adj(3)
        pos(4) = pos(4) - (oldsize(2) - newsize(2));
        if pos(4) < 1
            pos(4) = 1;
        end

    end
    set(h, 'position', pos);
else
    pos = [];
end


function resizetagdlg(fig,data,figorigpos)
hfig = figure(...
'Units','characters',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','Resize set anchors',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',get(0,'defaultfigurePaperSize'),...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[103.8 29.3076923076923 112 32.1538461538462],...
'Resize','on',...
'HandleVisibility','on',...
'UserData',struct('figorigpos',figorigpos),...
'Tag','figure1',...
'ResizeFcn',@resizefn,...
'CreateFcn', @(h,evnt)movegui(h,'center') ,...
'Visible','on');

uitable(...
'Parent',hfig,...
'Units','characters',...
'ColumnFormat',{  [] [] [] },...
'ColumnEditable',[false true false],...
'ColumnName',{  'tag'; 'anchors [left right top bott]'; 'origposition' },...
'ColumnWidth',{  166 166 166 },...
'Data',data,...
'Position',[2.4 5.23076923076923 105.4 26.1538461538462],...
'RearrangeableColumns','on',...
'RowName',blanks(0),...
'UserData',[],...
'Tag','datatable');

uicontrol(...
'Parent',hfig,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObj,ev)testdlg(hObj,ev,fig),...
'FontSize',10.6666666666667,...
'Position',[25 1.46153846153846 13.8 1.69230769230769],...
'String','Test dlg',...
'Tag','pushbutton1');

uicontrol(...
'Parent',hfig,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@okbtn,...
'FontSize',10.6666666666667,...
'Position',[48.8 1.46153846153846 13.8 1.69230769230769],...
'String','Ok',...
'Tag','pushbutton2');

uicontrol(...
'Parent',hfig,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObj,ev)delete(hfig),...
'FontSize',10.6666666666667,...
'Position',[73.4 1.46153846153846 13.8 1.69230769230769],...
'String','Cancel',...
'Tag','pushbutton3');

function okbtn(hObj,~)
hfig=getparentfig(hObj);
ud=get(hfig,'userdata');
h=findobj(hfig,'Tag','datatable');
data=get(h,'data');
disp('function resize(hobj,~)');
disp('hfig=getparentfig(hobj);');
disp('data={...');
for i=1:size(data,1)
    fprintf('''%s'',%s,%s;...\n',data{i,:});
end
fprintf('};\nnewpos=i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),%s);\nend %%maybe need to delete this\n',mat2str(ud.figorigpos));

function resizefn(hobj,~)
hfig=getparentfig(hobj);
data={'datatable',[1 1 1 1],[9 51 395.25 255];...
'pushbutton1',[1 0  0 1],[93.75 14.25 51.75 16.5];...
'pushbutton2',[1 0  0 1],[183 14.25 51.75 16.5];...
'pushbutton3',[1 0  0 1],[275.25 14.25 51.75 16.5]};
i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[301.5 122.25 420 313.5]);


function testdlg(hObj,~,fig)
hfig=getparentfig(hObj);
ud=get(hfig,'userdata');
h=findobj(hfig,'Tag','datatable');
data=get(h,'data');
for i=1:size(data,1)
    data{i,2}=str2num(data{i,2}); %#ok<ST2NM>
    data{i,3}=str2num(data{i,3}); %#ok<ST2NM>
end
set(fig,'ResizeFcn',@(~,~)i_resizedlg(fig,data(:,1),data(:,2),data(:,3),ud.figorigpos));
