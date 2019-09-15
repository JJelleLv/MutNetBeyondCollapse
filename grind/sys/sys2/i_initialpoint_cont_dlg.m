function  [pointlist,aval,afilter]=i_initialpoint_cont_dlg(obj, pointlist,aval,codim,afilter)
% --- Executes on button press in pushbutton10.
%function pushbutton10_Callback(hObject, eventdata, handles)
% eventdata  reserved - to be defined in a future version of MATLAB
% hObject    handle to pushbutton10 (see GCBO)
% handles    structure with handles and user data (see GUIDATA)
if nargin<1
    obj=evalin('base','g_cont');
end

if nargin<2
    pointlist= obj.show(1);
end

if nargin<3
    aval=1;
end
if nargin<4
    ud.codim=1;
else
    ud.codim=codim;
end
[~,~,filter]=obj.show;

if nargin<5
    if ~isempty(filter)
        afilter=filter(1);
    else
        afilter=[];
    end
end
strngs={'no filter'};
ud.filter=afilter;
filterval=1;
if ~isempty(filter)
    ud1.filter={[]};
    for j=1:length(filter)
        ud1.filter{end+1}=filter(j);
        if struccmp(filter(j),afilter)
            filterval=length(ud1.filter);
        end
        f=fieldnames(filter(j));
        if length(f)==1
            strngs{end+1}=sprintf('%s==%g',f{1},filter(j).(f{1}));
        else
            vals=struct2cell(filter(j));
            parvals=sprintf('%.4g,',vals{:});
            pars=sprintf('%s,',f{:});
            strngs{end+1}=sprintf('[%s]==[%s]',pars(1:end-1),parvals(1:end-1));
        end
    end
    vis='on';
else
    ud1.filter={};
    vis='off';
end
if isempty(aval)
    aval=1;
end
if length(aval)>1
    aval=aval(1);
end
if aval>length(pointlist)
    aval=length(pointlist);
end
ud=updateud(ud,pointlist,aval);

if ~isempty(aval)&&~isempty(ud.ids)
    data1=[transpose(obj.settings.derived.allpars),num2cell(obj.points(ud.pointndx(aval(1))).p0)];
    data2=[obj.settings.derived.statevars,num2cell(obj.points(ud.pointndx(aval(1))).x0)];
else
    %   val=[];
    data1={'',[]};
    data2={'',[]};
end
hfig = figure(...
    'PaperUnits','inches',...
    'Units','characters',...
    'Position',[135.8 19.9230769230769 123.2 31.2307692307692],...
    'Visible',get(0,'defaultfigureVisible'),...
    'Color',[0.914 0.914 0.914],...
    'IntegerHandle','off',...
    'MenuBar','none',...
    'Name','Initial conditions',...
    'NumberTitle','off',...
    'Tag','figure1',...
    'Resize','on',...
    'ResizeFcn',@resize,...
    'PaperPosition',get(0,'defaultfigurePaperPosition'),...
    'PaperSize',[8.5 11],...
    'PaperType','usletter',...
    'HandleVisibility','on',...
    'CreateFcn',@(h,evnt)movegui(h, 'center') );

uitable(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuitableFontUnits'),...
    'Units','characters',...
    'BackgroundColor',get(0,'defaultuitableBackgroundColor'),...
    'ColumnName',{  'Symbol'; 'Value' },...
    'ColumnWidth',{  100 100 },...
    'RowName','',...
    'Position',[5 5.85 43.6 18.2307692307692],...
    'ColumnEditable',[false true],...
    'ColumnFormat',{  [] [] },...
    'Data',data2,...
    'RearrangeableColumns',get(0,'defaultuitableRearrangeableColumns'),...
    'RowStriping',get(0,'defaultuitableRowStriping'),...
    'CellEditCallback',@cellEditvar,...
    'CellSelectionCallback',blanks(0),...
    'Children',[],...
    'ForegroundColor',get(0,'defaultuitableForegroundColor'),...
    'Enable',get(0,'defaultuitableEnable'),...
    'TooltipString','Edit values of state variables (only for EP or P points)',...
    'Visible',get(0,'defaultuitableVisible'),...
    'ButtonDownFcn',blanks(0),...
    'Tag','vartable');
uitable(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuitableFontUnits'),...
    'Units','characters',...
    'BackgroundColor',get(0,'defaultuitableBackgroundColor'),...
    'ColumnName',{  'Symbol'; 'Value' },...
    'ColumnWidth',{  100 100 },...
    'RowName','',...
    'Position',[51.4 5.85 43.6 18.2307692307692],...
    'ColumnEditable',[false true],...
    'ColumnFormat',{  [] [] },...
    'Data',data1,...
    'RearrangeableColumns',get(0,'defaultuitableRearrangeableColumns'),...
    'RowStriping',get(0,'defaultuitableRowStriping'),...
    'CellEditCallback',@cellEditpar,...
    'CellSelectionCallback',blanks(0),...
    'Children',[],...
    'ForegroundColor',get(0,'defaultuitableForegroundColor'),...
    'Enable',get(0,'defaultuitableEnable'),...
    'TooltipString','Edit parameter values (only for EP or P points)',...
    'Visible',get(0,'defaultuitableVisible'),...
    'ButtonDownFcn',blanks(0),...
    'DeleteFcn',blanks(0),...
    'Tag','partable',...
    'UserData',[],...
    'KeyPressFcn',blanks(0),...
    'HandleVisibility',get(0,'defaultuitableHandleVisibility'),...
    'FontSize',get(0,'defaultuitableFontSize'),...
    'FontName',get(0,'defaultuitableFontName'),...
    'FontAngle',get(0,'defaultuitableFontAngle'),...
    'FontWeight',get(0,'defaultuitableFontWeight'));


uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'BackgroundColor',[1 1 1],...
    'String',pointlist,...
    'Style','popupmenu',...
    'Value',aval,...
    'TooltipString','Select initial point for continuation',...
    'Callback',@popup_Callback,...
    'Position',[35 26.8 40 3],...
    'Children',[],...
    'Tag','PopupEqlist');


uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left',...
    'Units','characters',...
    'String','Filter:',...
    'Style','text',...
    'Position',[100 28 20 3],...
    'Visible',vis,...
    'Tag','FilterText');

uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'BackgroundColor',[1 1 1],...
    'Callback',@filter_Callback,...
    'String',strngs,...
    'Style','popupmenu',...
    'Value',filterval,...
    'Position',[100 26.8 20 3],...
    'Tooltipstring','Filter points with the same parameter setting',...
    'Visible',vis,...
    'Userdata',ud1,...
    'Tag','PopupFilter');

uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'Callback',@ok_Callback,...
    'String','OK',...
    'Style',get(0,'defaultuicontrolStyle'),...
    'Position',[30.2 2 19 2.23076923076923],...
    'Children',[],...
    'ButtonDownFcn',blanks(0),...
    'DeleteFcn',blanks(0),...
    'Tooltipstring','Save changes and close',...
    'Tag','OKbutton',...
    'KeyPressFcn',blanks(0));

uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'HorizontalAlignment','left',...
    'String','Select point',...
    'Style','text',...
    'Position',[5 27.8461538461538 22.6 1.07692307692308],...
    'Children',[],...
    'Tag','text2');


uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'String','Cancel',...
    'Style',get(0,'defaultuicontrolStyle'),...
    'Position',[56 2 19 2.23076923076923],...
    'Children',[],...
    'Callback',@cancel_Callback,...
    'Tooltipstring','Cancel, ignore changes',...
    'ButtonDownFcn',blanks(0),...
    'Tag','Cancelbutton',...
    'KeyPressFcn',blanks(0));


%         h7 = uicontrol(...
%             'Parent',hfig,...
%             'FontUnits',get(0,'defaultuicontrolFontUnits'),...
%             'Units','characters',...
%             'String','Search',...
%             'Style',get(0,'defaultuicontrolStyle'),...
%             'Position',[100.6 28.2307692307692 19 2.23076923076923],...
%             'Children',[],...
%             'Tag','pushbutton8',...
%             'KeyPressFcn',blanks(0));
dat=load('iconrefresh');
uicontrol('Parent',hfig, ...
    'Units','characters', ...
    'Callback',@refresh_Callback, ...
    'ListboxTop',0, ...
    'Position',[35+42 26.8+1 4 2], ...
    'cdata',dat.iconrefresh, ...
    'Tag','RefreshButton', ...
    'TooltipString','Refresh list of equilibria (list may be incomplete)');

uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'String','New IC point',...
    'Style',get(0,'defaultuicontrolStyle'),...
    'Position',[100.6 22 19 2.23],...
    'Callback',@newIC_Callback,...
    'TooltipString','Create new initial point (only P points)',...
    'Tag','NewButton',...
    'KeyPressFcn',blanks(0));
uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'Units','characters',...
    'String','Delete',...
    'Style',get(0,'defaultuicontrolStyle'),...
    'Position',[100.6 19 19 2.23],...
    'TooltipString','Delete current initial point',...
    'Callback',@delete_Callback,...
    'Tag','DeleteButton',...
    'KeyPressFcn',blanks(0));
uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'Units','characters',...
    'HorizontalAlignment','left',...
    'String','State variables:',...
    'Style','text',...
    'Position',[5 24.2 40.2 1.07692307692308],...
    'Children',[],...
    'ButtonDownFcn',blanks(0),...
    'DeleteFcn',blanks(0),...
    'Tag','text3');


uicontrol(...
    'Parent',hfig,...
    'FontUnits',get(0,'defaultuicontrolFontUnits'),...
    'BackgroundColor',[0.914 0.914 0.914], ...
    'Units','characters',...
    'HorizontalAlignment','left',...
    'String','Parameters:',...
    'Style','text',...
    'Position',[51 24.2 40.2 1.08],...
    'Children',[],...
    'ButtonDownFcn',blanks(0),...
    'DeleteFcn',blanks(0),...
    'Tag','partext');
set(hfig,'userdata',ud);
uiwait(hfig);
if ishandle(hfig)
    ud=get(hfig,'userdata');
    aval=ud.value;
    afilter=ud.filter;
    pointlist=ud.pointlist;
    close(hfig);
else
    pointlist=[];
    aval=[];
end

    function ud=updateud(ud,list,aval)
        if nargin>1&&(isempty(list)||all(strcmp(list,'None available')))
            ud.pointlist={'None available'};
            ud.ids={};
            ud.pointndx=0;
            ud.value=1;
            return;
        end
        if nargin>1
            ud.pointlist=list;
        end
        if nargin>2
            ud.value=aval;
        end
        if ud.value>length(ud.pointlist)
            ud.value=length(ud.pointlist);
        end
        ud.ids=regexp(ud.pointlist,'[A-Za-z+0-9]*','match','once');
        ud.pointndx=zeros(size(ud.ids));
        for i2=1:length(ud.ids)
            ud.pointndx(i2)=find(obj.getndx('points','id',ud.ids{i2}),1);
        end
    end
    function refresh_Callback(hObject, ~, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'Userdata');
        obj.add_points(findeqs('-a'));
        ud=updateud(ud, obj.show(ud.codim,ud.filter));
        
        h=findobj(hfig,'tag','PopupEqlist');
        set(h, 'string', ud.pointlist);
        set(h, 'Value', ud.value);
        
        set(hfig,'Userdata',ud)
    end

    function cellEditvar(hObject, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'Userdata');
        data=get(hObject,'data');
        vect=transpose([data{:,2}]);
        obj.points(ud.pointndx(ud.value)).x0=vect;
        s2 = sprintf('%g ', round(vect * 10000) / 10000);
        if strcmp(obj.pointprops(obj.points(ud.pointndx(ud.value)).propndx).label,'P')
            obj.points(ud.pointndx(ud.value)).data.msg=sprintf('initial point: %s',s2);
            h=findobj(hfig,'Tag','PopupEqlist');
            i_keep(vect);
            str=get(h,'string');
            str{ud.value}=sprintf('%s - %s',obj.points(ud.pointndx(ud.value)).id,obj.points(ud.pointndx(ud.value)).data.msg);
            set(h,'string',str);
        end
        %ud.points(ud.idndx).N0=vect;
        %set(hfig,'Userdata',ud);
    end

    function cellEditpar(hObject, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'Userdata');
        data=get(hObject,'data');
        vect=transpose(data{:,2});
        obj.points(ud.pointndx(ud.value)).p0=vect;
        
    end

    function cancel_Callback(hObject, ~, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'userdata');
        ud.pointlist=[];
        ud.value=[];
        set(hfig,'userdata',ud)
        uiresume;
    end

    function ok_Callback(hObject, ~, ~)
        hfig=get(hObject,'parent');
        ud = get(hfig, 'userdata');
        h=findobj(hfig,'Tag','PopupEqlist');
        ud.pointlist=get(h,'string');
        ud.value=get(h,'value');
        set(hfig,'userdata',ud)
        uiresume;
    end

    function popup_Callback(hObject, ~, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'Userdata');
        ud.value=get(hObject,'Value');
        selectpoint(hfig,ud);
        set(hfig,'Userdata',ud);
    end

    function filter_Callback(hObject,~,~)
        hfig=get(hObject,'parent');
        h=findobj(hfig,'Tag','PopupFilter');
        ud1=get(h,'Userdata');
        aval=get(h,'Value');
        
        
        ud=get(hfig,'Userdata');
        ud.filter=ud1.filter{aval};
        oldpoint='';
        if isfield(ud,'value')
            oldpoint=ud.pointlist{ud.value};
        end
        ud=updateud(ud, obj.show(ud.codim,ud.filter));
        f=find(strcmp(oldpoint,ud.pointlist));
        if isempty(f)
            ud.value=1;
        else
            ud.value=f;
        end
        selectpoint(hfig,ud);
        h=findobj(hfig,'tag','PopupEqlist');
        set(h, 'string', ud.pointlist);
        set(h, 'Value', ud.value);
        set(hfig,'Userdata',ud)
        
    end

    function delete_Callback(hObject, ~, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'Userdata');
        %id0=ud.points(ud.idndx).id;
        ud.pointlist(ud.value)=[];
        obj.points(ud.pointndx(ud.value))=[];
        ud=updateud(ud);
        h=findobj(hfig,'Tag','PopupEqlist');
        set(h,'Value',ud.value);
        set(h,'string',ud.pointlist);
        selectpoint(hfig,ud);
        set(hfig,'Userdata',ud);
    end

    function newIC_Callback(hObject, ~, ~)
        hfig=get(hObject,'parent');
        ud=get(hfig,'Userdata');
        
        if ~isempty(ud.ids)
            i=1;
            newEP=ud.ids{1};
            while any(strcmp(ud.ids,newEP))
                newEP=sprintf('P%d',i);
                i=i+1;
            end
        else
            newEP='P1';
        end
        N0=i_initvar;
        newpoint=struct('id',newEP,'label', 'P','msg', 'User-defined','data', struct('eigenval',[]),'N0', N0);
        
        %newpoint=struct('ndx',[],'labnr', 'user-defined','labname', 'IC','runid', '','id', newEP,'p0', p0,'N0', N0,'dim', 1,'used_conteq', false,'used_contbif',false,'xtra', []);
        obj.add_points(newpoint);
        ud.pointlist=obj.show(ud.codim);
        ud.value=numel(ud.pointlist);
        ud=updateud(ud);
        h=findobj(hfig,'Tag','PopupEqlist');
        set(h,'string',ud.pointlist);
        set(h,'Value',ud.value);
        selectpoint(hfig,ud);
        set(hfig,'Userdata',ud);
    end
    function resize(hobj,~)
        hfig=getparentfig(hobj);
        data={...
            'vartable',[1 1 1 1],[18.75 57.0375 163.5 177.75];...
            'partable',[0 1 1 1],[192.75 57.0375 163.5 177.75];...
            'PopupEqlist',[1 1 1 0],[131.25 261.3 150 29.25];...
            'FilterText',[0 1 1 0],[375 273 75 29.25];...
            'PopupFilter',[0 1 1 0],[375 261.3 75 29.25];...
            'OKbutton',[1 0 0 1 ],[113.25 19.5 71.25 21.75];...
            'text2',[1 0 1 0],[18.75 271.5 84.75 10.5];...
            'Cancelbutton',[1 0 0 1 ],[210 19.5 71.25 21.75];...
            'RefreshButton',[0 1 1 0],[288.75 271.05 15 19.5];...
            'NewButton',[0 1 1 0],[377.25 214.5 71.25 21.7425];...
            'DeleteButton',[0 1 1 0],[377.25 185.25 71.25 21.7425];...
            'text3',[1 0 1 0],[18.75 235.95 150.75 10.5];...
            'partext',[0 1 1 0],[191.25 235.95 150.75 10.53];...
            };
        i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[248.25 126.75 462 304.5]);
        res=i_resizedlg(hfig,'horizdistrib',{'vartable','partable'});
        h=findobj(hfig,'tag','partext');
        if ~isempty(h)
            pos=get(h,'position');
            pos(1)=pos(1)+res(1);
            set(h,'position',pos);
            i_resizedlg(hfig,'horizcenter',{'OKbutton','Cancelbutton'});
        end
    end

    function selectpoint(hfig,ud)
        
        if ~isempty(ud.value)
            if isempty(ud.ids)
                data1={'',[]};
                data2={'',[]};
            else
                id0=regexp(ud.ids(ud.value),'R[1-4]|[A-Za-z]*','match','once');
                data1=[transpose(obj.settings.derived.allpars),num2cell(obj.points(ud.pointndx(ud.value)).p0)];
                data2=[obj.settings.derived.statevars,num2cell(obj.points(ud.pointndx(ud.value)).x0)];
            end
            h=findobj(hfig,'tag','partable');
            if ~isempty(ud.ids)
                set(h, 'ColumnEditable',[false any(strcmp({'EP','P'},id0))]);
            end
            set(h,'data',data1);
            h=findobj(hfig,'tag','vartable');
            if ~isempty(ud.ids)
                set(h, 'ColumnEditable',[false any(strcmp({'EP','P'},id0))]);
            end
            set(h,'data',data2);
        end
    end
end
